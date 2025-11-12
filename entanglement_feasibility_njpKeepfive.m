function entanglement_feasibility_njp
% ================================================================
% Prospective feasibility simulation + robust fitting + UQ
% with Filter-function panel and CHSH feasibility
% ================================================================
rng(7,'twister');                       % reproducible

% ---- Global constants ------------------------------------------
G     = 6.67430e-11;                    % m^3 kg^-1 s^-2
hbar  = 1.054e-34;                      % J*s

% ---- Frequency / filter-function settings ----------------------
omega_g = 2*pi*1;                       % rad/s (sets lambda for common-mode)
T_filter = 1.0;                         % s (integration time used in filter panel)

% ---- PSD model for proper-time noise eps(t) --------------------
S_eps.type = 'white';                   % 'white' or 'pink'
S_eps.S0   = 1e-28;                     % ~best-clock order (Hz^-1) for white
S_eps.k    = 1e-34;                     % (for pink) one-sided k/f, illustrative
S_eps.flo  = 1e-4; S_eps.fhi = 1e2;     % numerical band for integral

% ---- Acquisition / noise settings ------------------------------
theta_K      = linspace(0,2*pi,12);     % analyzer phases
Nmean        = 150;                     % trials per (theta,t)
phase_jit_sd = 0.08;                    % rad (laser/phase)
contrast_jit_sd = 0.05;                 % slow drift on C
B_boot       = 500;                     % bootstrap replicates
SNR_min      = 1.5;                     % keep only points with V/se >= SNR_min

% ---- Platforms (toggle "heavy" for Option B) -------------------
include_heavy = false;                  % <- set true to add QGEM_heavy projection

Pspec = { ...
  {'QGEM',        1e-14, 2.5e-7, 2.0e-4, 3e-2, 0.85, 0.0}, ...
  {'MAQRO',       1e-17, 1.0e-7, 3.0e-4, 8e-4, 0.60, 0.0}, ...
  {'Levitated',   1e-17, 8.0e-8, 8.0e-4, 2.5e-1,0.55, 0.0} ...
};
if include_heavy
  Pspec{end+1} = {'QGEM_heavy', 1.0, 5e-5, 2.0e-2, 1e-3, 0.90, 0.0};
end

P = repmat(struct('name','','m',[],'dx',[],'r',[],'gamma_loc',[], ...
                  'C0',[],'phi0',[],'tset',[]), numel(Pspec), 1);
for i=1:numel(Pspec)
    P(i).name      = Pspec{i}{1};
    P(i).m         = Pspec{i}{2};
    P(i).dx        = Pspec{i}{3};
    P(i).r         = Pspec{i}{4};
    P(i).gamma_loc = Pspec{i}{5};
    P(i).C0        = Pspec{i}{6};
    P(i).phi0      = Pspec{i}{7};
end

% ---- Results collectors ----------------------------------------
budget_names = strings(0);
gamma_com_v  = [];  gamma_com_se = [];
gamma_loc_v  = [];  gamma_loc_se = [];
chsh_S       = [];  chsh_topt    = [];

% ==================== Figure A: Filter function =================
figure('Name','Figure A - Filter function & PSD','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% left: PSD
ff = logspace(log10(S_eps.flo), log10(S_eps.fhi), 2000);
Splot = S_of_eps(ff, S_eps);
nexttile; hold on; box on;
loglog(ff, Splot, 'LineWidth', 1.8);
xlabel('f (Hz)','Interpreter','tex'); ylabel('S_\epsilon(f) (Hz^{-1})','Interpreter','tex');
title('Proper-time PSD model','Interpreter','tex');

% right: |H(f)|^2 and integral value
H2 = (sin(pi*ff*T_filter)./(pi*ff)).^2;
Int = trapz(ff, Splot.*H2);
lambda = omega_g/2;
Vcm = exp(-(lambda^2)*Int);
gamma_com_panel = max(-log(max(Vcm,eps))/T_filter, 0);

nexttile; hold on; box on;
loglog(ff, H2./max(H2), 'LineWidth', 1.8);
xlabel('f (Hz)','Interpreter','tex'); ylabel('|H(f)|^2 (arb)','Interpreter','tex');
title(sprintf('Filter @ T=%.2gs;  \\gamma_{com}~%.2e s^{-1}', T_filter, gamma_com_panel), 'Interpreter','tex');

% ==================== Figure 1 & 2: Visibility + Budget =========
figure('Name','Figures 1-2 - Visibility & Dephasing budget','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% -- left two tiles: per-platform visibility fits ----------------
% Precompute common-mode rate (from filter map) used in sim+budget
gamma_com_global = gamma_from_S(@(f)S_of_eps(f,S_eps), omega_g, 1.0);

% storage for CHSH feasibility
S_cap = 2*sqrt(2); % Tsirelson bound

for u = 1:numel(P)
    % chi scaling ~ G m^2 dx^2 /(hbar r^3)
    chi_true = G*(P(u).m^2)*(P(u).dx^2)/(hbar*P(u).r^3);   % s^-1

    % time grid focusing around first fringe (which is huge for small chi)
    t_lim      = 0.6 * (pi/max(chi_true,1e-16)); % avoid inf when chi tiny
    base_grid  = linspace(0.05*t_lim, 1.2*t_lim, 6);
    dense_peak = linspace(0.4*t_lim, 0.9*t_lim, 6);
    tset       = unique(round([base_grid dense_peak],3));
    tset(~isfinite(tset)) = [];
    tset(tset<=0) = [];
    if isempty(tset), tset = linspace(0.5,2.0,6); end
    P(u).tset  = tset;

    % simulate counts
    Gamma_true = P(u).gamma_loc + gamma_com_global;
    [Theta,Tgrid] = ndgrid(theta_K, tset);
    K  = numel(theta_K);  J = numel(tset);
    Nshots = poissrnd(Nmean, K, J);
    C_t  = max(min(P(u).C0 + contrast_jit_sd*randn(1,J),0.999),0.05);
    phi_jit = phase_jit_sd*randn(K,J);
    p = 0.5*(1 + C_t(ones(K,1),:).*exp(-2*Gamma_true*Tgrid).* ...
              cos(2*chi_true*Tgrid + P(u).phi0 + Theta + phi_jit));
    p = min(max(p,1e-6),1-1e-6);
    nPlus = binornd(Nshots, p);

    % per-time complex visibility
    [Vhat, Vse, theta0_hat] = per_time_complex_visibility(nPlus, Nshots, theta_K);
    snr = Vhat./max(Vse,eps);
    keep = snr >= SNR_min;

    % unwrap phase & init for MLE
    [chi0, phi00] = fit_phase_vs_time(theta0_hat(keep), tset(keep), Vse(keep));
    C0  = median(Vhat(keep));
    G0  = max(0, robust_decay_init(Vhat(keep), tset(keep)));
    x0 = [logit(C0), phi00, log(max(G0,1e-6)), log(max(chi0,1e-10))];

    % MLE
    [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat] = fit_fringe_MLE(theta_K, tset, Nshots, nPlus, x0, true);

    % visibility band
    tt = linspace(min(tset), max(tset), 120);
    Vfit = C_hat*exp(-2*Gamma_hat*tt);
    band = visibility_band(tt, C_hat, Gamma_hat, Cov_nat([1 3],[1 3]));

    % bootstrap bands at measured t
    [Vlow, Vhigh] = bootstrap_band(B_boot, theta_K, tset, Nshots, nPlus, ...
                                   @(nk,Nk) per_time_complex_visibility(nk,Nk,theta_K), 0.16);

    % AIC (decay vs const)
    [nll1,k1] = model_nll(theta_K,tset,Nshots,nPlus,[C_hat,phi_hat,Gamma_hat,chi_hat]);
    [C0c,ph0c,chi_c] = constant_visibility_init(Vhat(keep),theta0_hat(keep),tset(keep));
    [nll0,k0] = model_nll(theta_K,tset,Nshots,nPlus,[C0c,ph0c,0,chi_c]);
    AIC1 = 2*k1 + 2*nll1; AIC0 = 2*k0 + 2*nll0; dAIC = AIC1 - AIC0;

    % ---- Plot per-platform panel (left two tiles combined) -----
    nexttile( (u<=2) + 1 ); % tile 1 and 2 for first two; 3rd platform overplots on tile 2
    hold on; box on;
    % param band
    fill([tt, fliplr(tt)],[band.low, fliplr(band.high)], [0.8 0.9 1.0], ...
        'EdgeColor','none','FaceAlpha',0.35);
    % bootstrap band at discrete t
    fill([tset, fliplr(tset)],[Vlow, fliplr(Vhigh)], [0.6 0.8 1.0], ...
        'EdgeColor','none','FaceAlpha',0.35);
    % fit curve
    plot(tt, Vfit, 'LineWidth', 2);
    % data
    errorbar(tset(keep), Vhat(keep), Vse(keep), 'o', 'MarkerSize', 5, ...
             'MarkerFaceColor',[0.9 0.4 0.1], 'Color',[0.9 0.4 0.1], ...
             'CapSize',0, 'LineStyle','none');
    scatter(tset(~keep), Vhat(~keep), 18, 'filled', 'MarkerFaceColor',[0.7 0.7 0.7]);

    xlabel('t (s)','Interpreter','tex'); ylabel('Visibility','Interpreter','tex');
    ylim([0, 1]);
    title(sprintf('%s',P(u).name),'Interpreter','tex');

    txt = sprintf('C=%.2f\\pm%.2f, \\Gamma=%.2e\\pm%.2e s^{-1}\\n\\chi=%.2e\\pm%.2e s^{-1}\\n\\DeltaAIC=%.1f',...
        C_hat, SE(1), Gamma_hat, SE(3), chi_hat, SE(4), dAIC);
    ax = gca; xP = ax.XLim(1) + 0.02*range(ax.XLim); yP = 0.08;
    text(xP, yP, txt, 'Units','normalized','VerticalAlignment','bottom','FontSize',9,'Interpreter','tex');

    % ---- Budget numbers ----------------------------------------
    budget_names(end+1,1) = string(P(u).name);
    gamma_com_v(end+1,1)  = gamma_com_global;
    gamma_com_se(end+1,1) = 0.05*gamma_com_global;  % nominal
    gamma_loc_v(end+1,1)  = max(Gamma_hat - gamma_com_global, 0);
    gamma_loc_se(end+1,1) = SE(3);

    % ---- CHSH feasibility (best-time projection) ----------------
    [Sproj, t_opt] = chsh_projection(C_hat, Gamma_hat, chi_hat, tt);
    chsh_S(end+1,1)    = Sproj;
    chsh_topt(end+1,1) = t_opt;

    fprintf('[%s] gamma_com %.3e s^-1, gamma_loc %.3e +/- %.2e s^-1; chi %.3e +/- %.2e s^-1; C %.2f +/- %.2f; Sproj=%.2f @ t=%.2fs\n',...
        P(u).name, gamma_com_global, gamma_loc_v(end), gamma_loc_se(end), chi_hat, SE(4), C_hat, SE(1), Sproj, t_opt);
end

% -- right tile: dephasing budget (stacked) ----------------------
nexttile; hold on; box on;
bh = barh(categorical(budget_names), [gamma_com_v gamma_loc_v], 'stacked');
legend({'\gamma_{com}','\gamma_{loc}'},'Location','southeast');
for i=1:numel(budget_names)
    x = gamma_com_v(i) + gamma_loc_v(i);
    errorbar(x, i, gamma_loc_se(i), 'horizontal', 'k', 'LineStyle','none', 'CapSize',8, 'LineWidth',1.2);
end
xlabel('rate (s^{-1})','Interpreter','tex');
title('Dephasing budget (68% error bars)','Interpreter','tex');

% ==================== Figure 3: Total decoherence rate ==========
figure('Name','Figure 3 - Total decoherence rate','Color','w');
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
tot = gamma_com_v + gamma_loc_v;
err = gamma_loc_se; % dominant
bar(categorical(budget_names), tot); hold on; box on;
for i=1:numel(tot)
    errorbar(i, tot(i), err(i), 'k','LineStyle','none','CapSize',8,'LineWidth',1.2);
end
ylabel('\Gamma_\phi (s^{-1})','Interpreter','tex'); title('Total dephasing rate \Gamma_\phi','Interpreter','tex');

% ==================== Figure 4: Scaling-law check ===============
figure('Name','Figure 4 - Scaling law check','Color','w');
scaling_check(G, hbar); % re-uses internal plotting

% ==================== Figure 5: CHSH feasibility ================
figure('Name','Figure 5 - CHSH feasibility','Color','w');
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
bar(categorical(budget_names), chsh_S); hold on; box on;
yline(2,'r--','LineWidth',1.5); yline(2*sqrt(2),'k:','LineWidth',1.2);
ylabel('S (CHSH)','Interpreter','tex'); ylim([0, 2.9]);
title('CHSH feasibility (best-time projected S)','Interpreter','tex');
text(0.02,0.92,'Red dashed: S=2 threshold; black dotted: 2\surd2','Units','normalized','FontSize',9,'Interpreter','tex');

end % ============================ main ===========================


% ======================= Helper functions ===============================

function S = S_of_eps(f, Smodel)
    switch lower(Smodel.type)
        case 'white'
            S = Smodel.S0 + 0*f;
        case 'pink'
            ff = max(min(f, Smodel.fhi), Smodel.flo);
            S = Smodel.k ./ max(ff, 1e-12);
        otherwise
            error('Unknown PSD type.');
    end
end

function gamma_com = gamma_from_S(Sfun, omega_g, t)
    ff  = logspace(-4, 3, 6000);
    H2  = (sin(pi*ff*t)./(pi*ff)).^2;
    Int = trapz(ff, Sfun(ff).*H2);
    lambda = omega_g/2;
    V = exp(-lambda^2 * Int);
    gamma_com = max(-log(max(V,eps))/t, 0);
end

function [Vhat, Vse, theta0] = per_time_complex_visibility(nPlus, Nshots, thetas)
    K = numel(thetas); J = size(nPlus,2);
    Vhat = zeros(1,J); theta0 = zeros(1,J); Vse = zeros(1,J);
    w = ones(K,1); w = w/sum(w);
    for j=1:J
        p_hat = nPlus(:,j)./max(Nshots(:,j),1);
        y = 2*p_hat - 1;
        z = sum(w .* y .* exp(-1i*thetas(:)));
        Vhat(j)   = abs(z);
        theta0(j) = angle(z);
        var_p = max(p_hat.*(1-p_hat)./max(Nshots(:,j),1), 1e-9);
        c = cos(thetas(:) - theta0(j));
        var_y = 4*sum((w.^2).*var_p.*(c.^2));
        Vse(j) = sqrt(var_y);
    end
    Vhat = max(min(Vhat,0.999), 0.001);
end

function [chi, phi0] = fit_phase_vs_time(phi_t, t, Vse)
    if isempty(phi_t)
        chi = 1e-10; phi0 = 0; return;
    end
    phi_u = unwrap(phi_t(:));
    W = diag(1./max(Vse(:),1e-3));
    X = [2*t(:), ones(numel(t),1)];
    b = (X'*W*X)\(X'*W*phi_u);
    chi  = max(b(1), 1e-10);
    phi0 = b(2);
end

function G0 = robust_decay_init(V, t)
    if isempty(V) || isempty(t)
        G0 = 1e-3; return;
    end
    V = min(max(V,1e-3),0.99);
    y = -log(V);
    X = [2*t(:), ones(numel(t),1)];
    b = X\y(:);
    G0 = max(b(1), 0);
end

function [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat] = fit_fringe_MLE(theta, t, Nshots, nPlus, x0, wantCov)
    if nargin<6, wantCov=false; end
    options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5e4,'MaxIter',5e3);
    nll = @(x) fringe_nll(x, theta, t, Nshots, nPlus);
    xhat = fminsearch(nll, x0, options);

    [C_hat, phi_hat, Gamma_hat, chi_hat] = decode_params(xhat);

    if nargout>=5
        H = num_hessian(nll, xhat);
        Cov = pinv(H);
        J = diag([C_hat*(1-C_hat), 1, Gamma_hat, chi_hat]); % d[natural]/d[x]
        Cov_nat = J * Cov * J';
        SE = sqrt(diag(Cov_nat));               % [SE_C, SE_phi, SE_Gamma, SE_chi]
    else
        SE=[]; Cov_nat=[];
    end
end

function val = fringe_nll(x, theta, t, Nshots, nPlus)
    [C, phi0, Gamma, chi] = decode_params(x);
    [Theta,T] = ndgrid(theta, t);
    p = 0.5*(1 + C*exp(-2*Gamma*T).*cos(2*chi*T + phi0 + Theta));
    p = min(max(p,1e-9),1-1e-9);
    val = -sum( nPlus(:).*log(p(:)) + (Nshots(:)-nPlus(:)).*log(1-p(:)) );
end

function H = num_hessian(fun, x)
    h = 1e-5;
    n = numel(x);
    H = zeros(n);
    for i=1:n
        for j=i:n
            ei = zeros(n,1); ej = zeros(n,1);
            ei(i)=1; ej(j)=1;
            fpp = fun(x + h*ei + h*ej);
            fpm = fun(x + h*ei - h*ej);
            fmp = fun(x - h*ei + h*ej);
            fmm = fun(x - h*ei - h*ej);
            H(i,j) = (fpp - fpm - fmp + fmm)/(4*h*h);
            H(j,i) = H(i,j);
        end
    end
    H = (H+H')/2;
end

function [C, phi0, Gamma, chi] = decode_params(x)
    C     = invlogit(x(1));
    phi0  = x(2);
    Gamma = exp(x(3));
    chi   = exp(x(4));
end

function y = invlogit(z), y = 1./(1+exp(-z)); end
function z = logit(p),   p = min(max(p,1e-9),1-1e-9); z = log(p./(1-p)); end

function band = visibility_band(tt, C, Gamma, Cov_CG)
    band.low = zeros(size(tt)); band.high = band.low;
    for k=1:numel(tt)
        t = tt(k);
        V = C*exp(-2*Gamma*t);
        g = [exp(-2*Gamma*t), -2*t*C*exp(-2*Gamma*t)];
        varV = max(g * Cov_CG * g.', 0);
        seV  = sqrt(varV);
        band.low(k)  = max(V - seV, 0);
        band.high(k) = min(V + seV, 1);
    end
end

function [Vlow, Vhigh] = bootstrap_band(B, theta, t, Nshots, nPlus, perTimeFn, q)
    if nargin<7, q=0.16; end
    J = numel(t);
    Vmat = zeros(B,J);
    for b=1:B
        phat = nPlus./max(Nshots,1);
        nb   = binornd(Nshots, min(max(phat,1e-6),1-1e-6));
        [Vb, ~] = perTimeFn(nb, Nshots);
        Vmat(b,:) = Vb;
    end
    Vlow  = quantile(Vmat, q, 1);
    Vhigh = quantile(Vmat, 1-q, 1);
end

function [nll,k] = model_nll(theta,t,Nshots,nPlus,pars)
    C=pars(1); phi0=pars(2); Gamma=pars(3); chi=pars(4);
    [Theta,T] = ndgrid(theta, t);
    p = 0.5*(1 + C*exp(-2*Gamma*T).*cos(2*chi*T + phi0 + Theta));
    p = min(max(p,1e-9),1-1e-9);
    nll = -sum( nPlus(:).*log(p(:)) + (Nshots(:)-nPlus(:)).*log(1-p(:)) );
    k   = 4 - (Gamma==0);
end

function [C0,phi0,chi0] = constant_visibility_init(Vhat,theta0,t)
    if isempty(Vhat)
        C0=0.5; phi0=0; chi0=1e-10; return;
    end
    C0   = max(min(median(Vhat),0.99), 0.05);
    phi0 = theta0(1);
    chi0 = max(polyfit(t(:), unwrap(theta0(:)), 1), 1e-10)/2;
end

function scaling_check(G,hbar)
    m  = logspace(-32,-14,40);
    dx = logspace(-16,-7,40);
    r  = logspace(-6,-3,40);
    [M,DX,R] = ndgrid(m, dx, r);
    chi = G.*(M.^2).*(DX.^2)./(hbar.*(R.^3));
    y = log10(chi(:));
    lm = log10(M(:)); ldx = log10(DX(:)); lr = log10(R(:));
    X = [ones(size(y)), lm, ldx, -lr];
    beta = X\y;
    b0=beta(1); a=beta(2); b=beta(3); c=beta(4);
    resid = y - X*beta;
    s2 = sum(resid.^2)/(numel(y)-size(X,2));
    Cov = s2 * inv(X.'*X);
    se  = sqrt(diag(Cov)); CI  = 1.96 * se(2:4);

    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; box on;
    scatter(lm, y, 6, 'filled');
    xx = sort(lm); yy = b0 + a*xx + b*mean(ldx) - c*mean(lr);
    plot(xx, yy, 'r-', 'LineWidth', 2);
    xlabel('log(m/kg)','Interpreter','tex'); ylabel('log(\chi/s^{-1})','Interpreter','tex');
    title(sprintf('a=%.2f \\pm %.2f (exp 2.0)', a, CI(1)),'Interpreter','tex');

    nexttile; hold on; box on;
    scatter(ldx, y, 6, 'filled');
    xx = sort(ldx); yy = b0 + a*mean(lm) + b*xx - c*mean(lr);
    plot(xx, yy, 'r-', 'LineWidth', 2);
    xlabel('log(\Delta x/m)','Interpreter','tex'); ylabel('log(\chi/s^{-1})','Interpreter','tex');
    title(sprintf('b=%.2f \\pm %.2f (exp 2.0)', b, CI(2)),'Interpreter','tex');

    nexttile; hold on; box on;
    scatter(lr, y, 6, 'filled');
    xx = sort(lr); yy = b0 + a*mean(lm) + b*mean(ldx) - c*xx;
    plot(xx, yy, 'r-', 'LineWidth', 2);
    xlabel('log(r/m)','Interpreter','tex'); ylabel('log(\chi/s^{-1})','Interpreter','tex');
    title(sprintf('c=%.2f \\pm %.2f (exp 3.0)', c, CI(3)),'Interpreter','tex');
end

function [Sproj, t_opt] = chsh_projection(C, Gamma, chi, tt)
    % Best-time CHSH under dephased-Bell family: S = 2*sqrt(2)*V(t),
    % choose t at nearest fringe maximum (2*chi*t = pi/2) if available.
    if chi <= 0
        t_opt = tt(end);
    else
        t_star = (pi/4)/chi;
        t_opt  = min(max(tt(1), t_star), tt(end));
    end
    V = C*exp(-2*Gamma*t_opt);
    Sproj = 2*sqrt(2)*max(min(V,1),0);
end
