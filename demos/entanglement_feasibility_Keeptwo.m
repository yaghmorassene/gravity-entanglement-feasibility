function entanglement_feasibility_demo8
% Prospective feasibility simulation with robust fitting + UQ + reporting
rng(7,'twister');                                % reproducible seed
run_id = datestr(now,'yyyymmdd_HHMMSS');         % run tag for logs/files
fprintf('== entanglement_feasibility_demo8 | run %s ==\n', run_id);

% ------------- Presentation / Export -------------
use_log_budget = true;
export_figs    = true;
export_tables  = true;
outdir         = 'entanglement_outputs';
scenario      = 'B';      % 'A' comparable gamma_com, 'B' subdominant/paper

set(groot,'defaultAxesFontName','Times', ...
          'defaultTextFontName','Times', ...
          'defaultAxesFontSize',11, ...
          'defaultTextInterpreter','latex', ...
          'defaultLegendInterpreter','latex', ...
          'defaultLineLineWidth',1.6);

% ------------- Physics / Platforms --------------
G      = 6.67430e-11;
hbar   = 1.054e-34;
omega_g = 2*pi*1;               % 1 Hz
t_filter = 1.0;                 % filter-function representative time

P = repmat(struct('name',"", 'm',0, 'dx',0, 'r',0, ...
                  'gamma_loc',0, 'C0',0, 'phi0',0), 1, 3);
P(1) = struct('name',"QGEM",     'm',1e-14, 'dx',2.5e-7, 'r',2e-4,  'gamma_loc',3e-2,  'C0',0.85, 'phi0',0.0);
P(2) = struct('name',"MAQRO",    'm',1e-17, 'dx',1.0e-7, 'r',3.0e-4,'gamma_loc',8e-4, 'C0',0.60, 'phi0',0.0);
P(3) = struct('name',"Levitated",'m',1e-17, 'dx',8.0e-8, 'r',8.0e-4,'gamma_loc',2.5e-1,'C0',0.55, 'phi0',0.0);

% ------------- Acquisition -----------------------
theta_K      = linspace(0,2*pi,12);
Nmean        = 150;
phase_jit_sd = 0.08;
contrast_jit_sd = 0.05;
B_boot       = 500;
SNR_min      = 1.5;

% ------------- PSD + Scenario toggle -------------
S_eps.type = 'white';
S_eps.S0   = 1e-32;          % baseline "good clock"
S_eps.k    = 1e-34;
S_eps.flo  = 1e-3; S_eps.fhi = 1e2;

if scenario=='A'
    gamma_com_target = 5e-1;                      % 0.5 s^-1 stress test
    S_eps.type = 'white';
    S_eps.S0   = 8*gamma_com_target/(omega_g^2);
elseif scenario=='B'
    % keep "good clock"
else
    error('scenario must be ''A'' or ''B''');
end
gamma_com_expected = (omega_g^2/8)* (strcmpi(S_eps.type,'white')*S_eps.S0 + 0);
fprintf('Scenario %s | S0=%.3e -> expected gamma_com=%.3e s^-1\n',...
        scenario, S_eps.S0, gamma_com_expected);

% ------------- Collectors for Table III ----------
names = strings(0);  C_hat_v=[]; SE_C=[]; Gamma_hat_v=[]; SE_Gamma=[]; chi_hat_v=[]; SE_chi=[];
gamma_com_v=[]; gamma_loc_v=[]; gamma_tot_v=[];
Smax_v = []; t_star_v=[]; violatesCHSH = [];

% ------------- FILTER FUNCTION PLOT (Figure A) ---
fig_filter = filter_function_plot(S_eps, omega_g, t_filter, outdir);

% ------------- VISIBILITY + FITS per platform ----
fig_vis = figure('Name',['Figure 1: Visibility by time (' scenario ')'],'Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for u=1:numel(P)
    % gravitational chi
    chi_true = G*(P(u).m^2)*(P(u).dx^2)/(hbar*P(u).r^3);

    % gamma_com from PSD (white analytic, pink numeric)
    gamma_com = gamma_from_psd(S_eps, omega_g, t_filter);

    % time grid (around first fringe)
    t_lim      = 0.6 * (pi/chi_true);
    base_grid  = linspace(0.05*t_lim, 1.2*t_lim, 6);
    dense_peak = linspace(0.4*t_lim, 0.9*t_lim, 6);
    tset       = unique(round([base_grid dense_peak],3));

    % simulate counts
    Gamma_true = P(u).gamma_loc + gamma_com;
    [Theta,Tgrid] = ndgrid(theta_K, tset); %#ok<NASGU>
    K=numel(theta_K); J=numel(tset);
    Nshots = poissrnd(Nmean, K, J);
    C_t  = max(min(P(u).C0 + contrast_jit_sd*randn(1,J),0.999),0.05);
    phi_jit = phase_jit_sd*randn(K,J);
    p = 0.5*(1 + C_t(ones(K,1),:).*exp(-2*Gamma_true*Tgrid).* ...
              cos(2*chi_true*Tgrid + P(u).phi0 + theta_K(:) + phi_jit));
    p = min(max(p,1e-6),1-1e-6);
    nPlus = binornd(Nshots, p);

    % per-time visibility/phase
    [Vhat,Vse,theta0_hat] = per_time_complex_visibility(nPlus, Nshots, theta_K);
    snr  = Vhat./max(Vse,eps);
    keep = snr >= SNR_min;
    if sum(keep)<3, [~,ord]=sort(snr,'descend'); keep=false(size(snr)); keep(ord(1:3))=true; end

    % robust inits
    [chi0, phi00] = fit_phase_vs_time(theta0_hat(keep), tset(keep), Vse(keep));
    C0  = median(Vhat(keep));
    G0  = max(0, robust_decay_init(Vhat(keep), tset(keep)));

    % MLE & Hessian/parametric UQ
    x0 = [logit(C0), phi00, log(max(G0,1e-6)), log(max(chi0,1e-8))];
    [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat, okH] = ...
        fit_fringe_MLE_robust(theta_K, tset, Nshots, nPlus, x0);

    badSE = any(~isfinite(SE)) || SE(1)>1 || SE(3)>10*max(Gamma_hat,1e-12) || SE(4)>10*max(chi_hat,1e-12);
    if ~okH || badSE
        Pse = parametric_SEs(300, theta_K, tset, Nshots, C_hat, phi_hat, Gamma_hat, chi_hat);
        SE = Pse; Cov_nat = diag(Pse.^2);
    end
    SE = max(SE, [0.02; 1e-3; 0.20*max(Gamma_hat,1e-12); 0.20*max(chi_hat,1e-12)]);

    % predicted band + bootstrap
    tt = linspace(min(tset), max(tset), 120);
    Vfit = C_hat*exp(-2*Gamma_hat*tt);
    band = visibility_band(tt, C_hat, Gamma_hat, Cov_nat([1 3],[1 3]));
    [Vlow,Vhigh] = bootstrap_band(B_boot, theta_K, tset, Nshots, nPlus, ...
                                  @(nk,Nk) per_time_complex_visibility(nk,Nk,theta_K), 0.16);

    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    fill([tt, fliplr(tt)],[band.low, fliplr(band.high)], [0.8 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.35);
    fill([tset, fliplr(tset)],[Vlow, fliplr(Vhigh)], [0.6 0.8 1.0], 'EdgeColor','none','FaceAlpha',0.35);
    plot(tt, Vfit, 'k-');
    errorbar(tset(keep), Vhat(keep), Vse(keep), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor',[0.9 0.4 0.1], 'Color',[0.9 0.4 0.1], ...
        'CapSize',0, 'LineStyle','none');
    scatter(tset(~keep), Vhat(~keep), 18, 'filled', 'MarkerFaceColor',[0.7 0.7 0.7]);
    xlabel('t (s)'); ylabel('Visibility'); ylim([0,1]); xlim([min(tset)*0.9, max(tset)*1.05]);
    title(P(u).name);

    txt = sprintf( ...
        ['$C=%.2f\\,\\pm\\,%.2f,\\; \\Gamma=%.2e\\,\\pm\\,%.2e\\;\\mathrm{s}^{-1}$\n' ...
         '$\\chi=%.2e\\,\\pm\\,%.2e\\;\\mathrm{s}^{-1}$'], ...
        C_hat, SE(1), Gamma_hat, SE(3), chi_hat, SE(4));
    text(0.05, 0.08, txt, 'Units','normalized', ...
        'VerticalAlignment','bottom','FontSize',9);

    % budget numbers
    names(end+1,1) = string(P(u).name);
    gamma_com_v(end+1,1)  = gamma_com;
    gamma_loc_v(end+1,1)  = max(Gamma_hat - gamma_com, 0);
    gamma_tot_v(end+1,1)  = gamma_com_v(end) + gamma_loc_v(end);
    C_hat_v(end+1,1)      = C_hat;   SE_C(end+1,1) = SE(1);
    Gamma_hat_v(end+1,1)  = Gamma_hat;  SE_Gamma(end+1,1)=SE(3);
    chi_hat_v(end+1,1)    = chi_hat; SE_chi(end+1,1)=SE(4);

    % -------- CHSH feasibility (conservative bound) --------
    % For sinusoidal correlations with visibility V_corr <= fitted envelope,
    % S_max = 2*sqrt(2)* V_corr. Violation requires V_corr > 1/sqrt(2).
    Vcorr_t = C_hat*exp(-2*Gamma_hat*tt);             % upper bound on 2-qubit visibility
    Smax_t  = 2*sqrt(2)*Vcorr_t;
    [Smax, idx] = max(Smax_t); t_star = tt(idx);
    Smax_v(end+1,1) = Smax; t_star_v(end+1,1)=t_star; violatesCHSH(end+1,1)= Smax>2;

    % Optional: small inset marker at t* where Smax occurs
    plot(t_star, Vfit(idx), 'ks','MarkerFaceColor','y','MarkerSize',6);
end

% ------------- Dephasing budget (Fig 2) ------------
fig_budget = figure('Name',['Figure 2: Dephasing budget (68% CI) (' scenario ')'],'Color','w');
cats = categorical(names);
barh(cats, [gamma_com_v, gamma_loc_v], 'stacked'); hold on; box on;
set(gca,'TickLabelInterpreter','latex');
legend({'$\gamma_{\mathrm{com}}$','$\gamma_{\mathrm{loc}}$'},'Location','southeast');
xlabel('rate (s$^{-1}$)'); title('Posterior dephasing budget (68\% error bars)');
x_end = gamma_com_v + gamma_loc_v;
for i=1:numel(names)
    % conservative SE for loc (from Hessian already captured via SE_Gamma)
    err = SE_Gamma(i);
    errorbar(x_end(i), i, err, 'horizontal','k','LineStyle','none','CapSize',8);
end
if use_log_budget
    set(gca,'XScale','log'); xmin = max(min([gamma_com_v; x_end]/2), 1e-12);
    xmax = max(x_end + max(SE_Gamma,1e-18)); xlim([xmin, 10^(ceil(log10(xmax)))]);
end

% ------------- Total dephasing (Fig 3) -------------
fig_total = figure('Name',['Figure 3: Total decoherence rate (' scenario ')'],'Color','w');
barh(cats, gamma_tot_v, 'FaceColor',[0.8 0.45 0.2]); hold on; box on;
set(gca,'TickLabelInterpreter','latex');
% independent-propagation (com vs loc)
gamma_tot_se = sqrt( (0.05*gamma_com_v).^2 + (SE_Gamma).^2 );
for i=1:numel(names)
    errorbar(gamma_tot_v(i), i, gamma_tot_se(i), 'horizontal','k','LineStyle','none','CapSize',8);
end
xlabel('rate (s$^{-1}$)'); title('Total decoherence $\gamma_{\mathrm{tot}}=\gamma_{\mathrm{com}}+\gamma_{\mathrm{loc}}$');
if use_log_budget
    set(gca,'XScale','log'); xmin = max(min(gamma_tot_v/2), 1e-12);
    xmax = max(gamma_tot_v + gamma_tot_se); xlim([xmin, 10^(ceil(log10(xmax)))]);
end

% ------------- CHSH quick-check (Fig 5) -------------
fig_chsh = figure('Name',['Figure 5: CHSH feasibility (' scenario ')'],'Color','w');
bar(categorical(names), Smax_v); hold on; yline(2,'r--','LineWidth',1.4);
ylabel('$S_{\max}=2\sqrt{2}\,V_{\mathrm{corr}}$'); title('CHSH bound from fitted envelope');
text(0.5,0.92,sprintf('violate?  %s', join(string(violatesCHSH), ', ')), ...
    'Units','normalized','HorizontalAlignment','center');

% ------------- Scaling law (Fig 4) -----------------
fig_scal = scaling_check(G,hbar);
set(fig_scal,'Name',['Figure 4: Scaling check (' scenario ')']);

% ------------- Export figures & tables --------------
if export_figs
    if ~exist(outdir,'dir'), mkdir(outdir); end
    savepdf(fig_filter, fullfile(outdir, ['figA_filter_' scenario '.pdf']));
    savepdf(fig_vis,    fullfile(outdir, ['fig1_visibility_' scenario '.pdf']));
    savepdf(fig_budget, fullfile(outdir, ['fig2_budget_' scenario '.pdf']));
    savepdf(fig_total,  fullfile(outdir, ['fig3_total_' scenario '.pdf']));
    savepdf(fig_scal,   fullfile(outdir, ['fig4_scaling_' scenario '.pdf']));
    savepdf(fig_chsh,   fullfile(outdir, ['fig5_chsh_' scenario '.pdf']));
end

if export_tables
    T = table(names, C_hat_v, SE_C, Gamma_hat_v, SE_Gamma, chi_hat_v, SE_chi, ...
              gamma_com_v, gamma_loc_v, gamma_tot_v, Smax_v, t_star_v, violatesCHSH, ...
              'VariableNames',{'Platform','C','SE_C','Gamma','SE_Gamma','chi','SE_chi', ...
                               'gamma_com','gamma_loc','gamma_tot','Smax','t_star','CHSH_violates'});
    if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(T, fullfile(outdir, ['TableIII_' scenario '_' run_id '.csv']));
    % LaTeX table
    texfile = fullfile(outdir, ['TableIII_' scenario '_' run_id '.tex']);
    fid=fopen(texfile,'w');
    fprintf(fid,'\\begin{tabular}{lrrrrrrrrr}\n\\toprule\n');
    fprintf(fid,'Platform & $C$ & $\\Gamma$ (s$^{-1}$) & $\\chi$ (s$^{-1}$) & $\\gamma_{com}$ & $\\gamma_{loc}$ & $\\gamma_{tot}$ & $S_{\\max}$ & $t^*$ (s) & CHSH\\\\\\midrule\n');
    for i=1:height(T)
        fprintf(fid,'%s & %.2f$\\pm$%.2f & %.2e$\\pm$%.1e & %.2e$\\pm$%.1e & %.2e & %.2e & %.2e & %.2f & %.2g & %s \\\\\n', ...
            T.Platform{i}, T.C(i), T.SE_C(i), T.Gamma(i), T.SE_Gamma(i), T.chi(i), T.SE_chi(i), ...
            T.gamma_com(i), T.gamma_loc(i), T.gamma_tot(i), T.Smax(i), T.t_star(i), logical2yesno(T.CHSH_violates(i)));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n'); fclose(fid);
    fprintf('Saved Table III to:\n  %s\n  %s\n', ...
        fullfile(outdir, ['TableIII_' scenario '_' run_id '.csv']), texfile);
end
end
% ============================= HELPER FUNCS ===============================

function fig = filter_function_plot(S_eps, omega_g, t, outdir)
    % Plots S_e(f), |H(f)|^2, product, and prints the integral & gamma_com
    ff = logspace(-4, 3, 6000);
    H2 = (sin(pi*ff*t)./(pi*ff)).^2;                        % |H|^2
    if strcmpi(S_eps.type,'white')
        S = S_eps.S0 + 0*ff;
    else
        fcl = max(min(ff,S_eps.fhi),S_eps.flo);
        S = S_eps.k./max(fcl,1e-18);
    end
    integrand = S.*H2; Int = trapz(ff, integrand);
    lambda = omega_g/2; V = exp(-lambda^2*Int);
    gamma_com = max(-log(max(V,eps))/t, 0);

    fig = figure('Name','Figure A: Filter function & PSD','Color','w');
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    nexttile; loglog(ff,S,'k-'); grid on; xlabel('$f$ (Hz)'); ylabel('$S_\epsilon(f)$ (s)');
    title('Clock PSD');
    nexttile; loglog(ff,H2,'b-'); grid on; xlabel('$f$ (Hz)'); ylabel('$|H(f)|^2$');
    title(sprintf('Filter, t=%.2f s',t));
    nexttile; loglog(ff, integrand,'m-'); grid on; xlabel('$f$ (Hz)'); ylabel('$S\! \cdot |H|^2$');
    title(sprintf('Int=%.3e ; \\gamma_{com}=%.3e s^{-1}', Int, gamma_com));
    if nargin>3 && ~isempty(outdir) && exist(outdir,'dir')
        savepdf(fig, fullfile(outdir,'figA_filter_preview.pdf'));
    end
end

function gamma_com = gamma_from_psd(S_eps, omega_g, t)
    if strcmpi(S_eps.type,'white')
        gamma_com = (omega_g^2/8) * S_eps.S0;
    else
        ff = logspace(-4, 3, 6000);
        H2 = (sin(pi*ff*t)./(pi*ff)).^2;
        fcl = max(min(ff,S_eps.fhi),S_eps.flo);
        S = S_eps.k./max(fcl,1e-18);
        Int = trapz(ff, S.*H2);
        lambda = omega_g/2;
        V = exp(-lambda^2*Int);
        gamma_com = max(-log(max(V,eps))/t, 0);
    end
end

function savepdf(fig, path), print(fig, path, '-dpdf','-painters'); end

function s = logical2yesno(x), s = char("No" + (x~=0)*"Yes"); end

% ---- The rest (MLE, vis, scaling, etc.) are identical to demo7 --------
% (copied here verbatim for a single-file drop-in)

function [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat, okH] = fit_fringe_MLE_robust(theta, t, Nshots, nPlus, x0)
    options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5e4,'MaxIter',5e3);
    nll = @(x) fringe_nll(x, theta, t, Nshots, nPlus);
    xhat = fminsearch(nll, x0, options);
    [C_hat, phi_hat, Gamma_hat, chi_hat] = decode_params(xhat);
    H = num_hessian(nll, xhat) + 1e-10*eye(4);
    s = svd(H); condH = s(1)/max(s(end),eps); okH = isfinite(condH) && condH < 1e8;
    Cov = pinv(H);
    J = diag([C_hat*(1-C_hat), 1, Gamma_hat, chi_hat]);
    Cov_nat = J * Cov * J';
    SE = sqrt(max(diag(Cov_nat),0));
end

function Pse = parametric_SEs(B, theta, t, Nshots, C, phi0, Gamma, chi)
    est  = zeros(B,4);
    [Theta,T] = ndgrid(theta,t); %#ok<NASGU>
    for b=1:B
        p = 0.5*(1 + C*exp(-2*Gamma*T).*cos(2*chi*T + phi0 + theta(:)));
        p = min(max(p,1e-8),1-1e-8);
        nb = binornd(Nshots, p);
        x0 = [logit(C), phi0, log(Gamma), log(chi)];
        [C_b, phi_b, G_b, chi_b] = fit_fringe_MLE_only(theta, t, Nshots, nb, x0);
        est(b,:) = [C_b, phi_b, G_b, chi_b];
    end
    Pse = std(est,0,1).';
end

function [C_hat, phi_hat, Gamma_hat, chi_hat] = fit_fringe_MLE_only(theta, t, Nshots, nPlus, x0)
    options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5e4,'MaxIter',5e3);
    nll = @(x) fringe_nll(x, theta, t, Nshots, nPlus);
    xhat = fminsearch(nll, x0, options);
    [C_hat, phi_hat, Gamma_hat, chi_hat] = decode_params(xhat);
end

function val = fringe_nll(x, theta, t, Nshots, nPlus)
    [C, phi0, Gamma, chi] = decode_params(x);
    [~,T] = ndgrid(theta, t);
    p = 0.5*(1 + C*exp(-2*Gamma*T).*cos(2*chi*T + phi0 + theta(:)));
    p = min(max(p,1e-9),1-1e-9);
    val = -sum( nPlus(:).*log(p(:)) + (Nshots(:)-nPlus(:)).*log(1-p(:)) );
end

function H = num_hessian(fun, x)
    n = numel(x); H = zeros(n);
    h = 1e-5*(1+abs(x));
    for i=1:n
        ei = zeros(n,1); ei(i)=1;
        for j=i:n
            ej = zeros(n,1); ej(j)=1;
            fpp = fun(x + h(i)*ei + h(j)*ej);
            fpm = fun(x + h(i)*ei - h(j)*ej);
            fmp = fun(x - h(i)*ei + h(j)*ej);
            fmm = fun(x - h(i)*ei - h(j)*ej);
            H(i,j) = (fpp - fpm - fmp + fmm)/(4*h(i)*h(j));
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
    phi_u = unwrap(phi_t(:));
    t = t(:); t0 = mean(t); ts = std(t)+eps;
    X = [2*(t - t0)/ts, ones(numel(t),1)];
    W = diag(1./max(Vse(:),1e-3));
    b = (X'*W*X)\(X'*W*phi_u);
    chi  = max(b(1)/ts, 1e-8);
    phi0 = b(2) - 2*(-t0/ts)*b(1);
end

function G0 = robust_decay_init(V, t)
    if numel(V) < 2 || numel(unique(t)) < 2, G0 = 0; return; end
    V = min(max(V,1e-3),0.99);
    y = -log(V);
    t = t(:); t0 = mean(t);
    X = [2*(t - t0), ones(numel(t),1)];
    if rank(X) < 2, G0 = 0; return; end
    b = X\y(:);
    G0 = max(b(1), 0);
end

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
        [Vb, ~] = perTimeFn(nb, Nshots); %#ok<ASGLU>
        Vmat(b,:) = Vb;
    end
    Vlow  = quantile(Vmat, q, 1);
    Vhigh = quantile(Vmat, 1-q, 1);
end

function fig = scaling_check(G,hbar)
    m  = logspace(-32,-14,40);
    dx = logspace(-16,-7,40);
    r  = logspace(-6,-3,40);
    [M,DX,R] = ndgrid(m, dx, r);
    chi = G.*(M.^2).*(DX.^2)./(hbar.*(R.^3));
    y = log10(chi(:));
    lm = log10(M(:)); ldx = log10(DX(:)); lr = log10(R(:));
    lm_c  = lm  - mean(lm); ldx_c = ldx - mean(ldx); lr_c  = lr  - mean(lr);
    X = [ones(size(y)), lm_c, ldx_c, -lr_c];
    [Q,Rq] = qr(X,0); beta = Rq\(Q'*y);
    b0=beta(1); a=beta(2); b=beta(3); c=beta(4);
    resid = y - X*beta;
    s2 = sum(resid.^2)/(numel(y)-size(X,2));
    Cov = s2 * inv(Rq'*Rq);
    se  = sqrt(max(diag(Cov),0));
    CI  = 1.96 * se(2:4);

    fig = figure('Name','Figure 4: Scaling check','Color','w');
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(lm, y, 6, 'filled'); xx = sort(lm);
    plot(xx, b0 + a*(xx-mean(lm)), 'r-');
    xlabel('$\log_{10}(m/\mathrm{kg})$'); ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('a=%.2f \\pm %.2f (exp 2.0)', a, CI(1)));
    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(ldx, y, 6, 'filled'); xx = sort(ldx);
    plot(xx, b0 + b*(xx-mean(ldx)), 'r-');
    xlabel('$\log_{10}(\Delta x/\mathrm{m})$'); ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('b=%.2f \\pm %.2f (exp 2.0)', b, CI(2)));
    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(lr, y, 6, 'filled'); xx = sort(lr);
    plot(xx, b0 - c*(xx-mean(lr)), 'r-');
    xlabel('$\log_{10}(r/\mathrm{m})$'); ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('c=%.2f \\pm %.2f (exp 3.0)', c, CI(3)));
end
