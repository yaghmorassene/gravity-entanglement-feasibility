function entanglement_feasibility_demo7
rng(7,'twister');

% ========================= Presentation / Export =========================
use_log_budget = true;            % log x-axis for rate panels
export_figs    = true;            % save publication-ready PDFs
outdir         = 'entanglement_outputs';
compare_AB     = false;           % true: render A and B side-by-side

% Consistent, publication-ish defaults
set(groot,'defaultAxesFontName','Times', ...
          'defaultTextFontName','Times', ...
          'defaultAxesFontSize',11, ...
          'defaultTextInterpreter','latex', ...
          'defaultLegendInterpreter','latex', ...
          'defaultLineLineWidth',1.6);

% ============================= Physics setup =============================
G      = 6.67430e-11;
hbar   = 1.054e-34;
omega_g = 2*pi*1;        % 1 Hz gate frequency

% Platforms
P = repmat(struct('name',"", 'm',0, 'dx',0, 'r',0, ...
                  'gamma_loc',0, 'C0',0, 'phi0',0), 1, 3);
P(1) = struct('name',"QGEM",     'm',1e-14, 'dx',2.5e-7, 'r',2e-4,  'gamma_loc',3e-2,  'C0',0.85, 'phi0',0.0);
P(2) = struct('name',"MAQRO",    'm',1e-17, 'dx',1.0e-7, 'r',3.0e-4,'gamma_loc',8e-4, 'C0',0.60, 'phi0',0.0);
P(3) = struct('name',"Levitated",'m',1e-17, 'dx',8.0e-8, 'r',8.0e-4,'gamma_loc',2.5e-1,'C0',0.55, 'phi0',0.0);

% Acquisition
theta_K      = linspace(0,2*pi,12);
Nmean        = 150;
phase_jit_sd = 0.08;
contrast_jit_sd = 0.05;
B_boot       = 500;
SNR_min      = 1.5;

% ===================== Scenario selection / comparison ===================
if compare_AB
    run_scenario('A');
    run_scenario('B');
else
    % ======= Scenario toggle (pick 'A' or 'B') ==========================
    run_scenario('B');  % 'A' comparable  |  'B' subdominant (paper)
end

% ============================== Export ==================================
if export_figs
    if ~exist(outdir,'dir'), mkdir(outdir); end
    drawnow;
    figs = findobj('Type','figure');
    for k=1:numel(figs)
        nm = string(get(figs(k),'Name'));
        if strlength(nm)==0, nm = "figure_"+k; end
        print(figs(k), fullfile(outdir, sanitize_filename(nm)+".pdf"), '-dpdf','-painters');
    end
    fprintf('Saved PDFs to folder "%s".\n', outdir);
end

% ============================ Nested runner ==============================
    function run_scenario(scenario)
        % ===== PSD (white by default; pink supported) ===================
        S_eps.type = 'white';
        S_eps.S0   = 1e-32;      % good clock baseline (Scenario B)
        S_eps.k    = 1e-34;      % for pink
        S_eps.flo  = 1e-3; S_eps.fhi = 1e2;

        if scenario=='A'
            % Target comparable gamma_com
            gamma_com_target = 5e-1;                            % 0.5 s^-1
            S_eps.type = 'white';
            S_eps.S0   = 8*gamma_com_target/(omega_g^2);
        elseif scenario=='B'
            % keep baseline above
        else
            error('scenario must be ''A'' or ''B''');
        end
        fprintf('Scenario %s: using S0 = %.3e s -> expected gamma_com = %.3e s^-1\n',...
                scenario, S_eps.S0, (omega_g^2/8)*S_eps.S0 );

        suffix = " ("+scenario+")";

        % Collectors
        budget_names = strings(0);
        gamma_com_v  = [];  gamma_com_se = [];
        gamma_loc_v  = [];  gamma_loc_se = [];

        % ---------- Visibility by time ----------
        fig_vis = figure('Name','Visibility by time'+suffix,'Color','w');
        tl = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

        for u=1:numel(P)
            chi_true = G*(P(u).m^2)*(P(u).dx^2)/(hbar*P(u).r^3);

            % robust gamma_com
            switch lower(S_eps.type)
                case 'white'
                    gamma_com = (omega_g^2/8)*S_eps.S0;
                case 'pink'
                    H2 = @(f,t) (sin(pi*f*t)./(pi*f)).^2;
                    integrand = @(f,t) (S_eps.k./max(f,1e-18)).*H2(f,t);
                    Int = integral(@(f)integrand(f,1.0), S_eps.flo, S_eps.fhi,'RelTol',1e-10,'AbsTol',1e-20);
                    lambda = omega_g/2;
                    V = exp(-lambda^2 * Int);
                    gamma_com = max(-log(max(V,eps))/1.0, 0);
                otherwise
                    error('Unknown PSD type');
            end

            % time grid
            t_lim      = 0.6 * (pi/chi_true);
            base_grid  = linspace(0.05*t_lim, 1.2*t_lim, 6);
            dense_peak = linspace(0.4*t_lim, 0.9*t_lim, 6);
            tset       = unique(round([base_grid dense_peak],3));

            % simulate
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

            % per-time vis/phase
            [Vhat,Vse,theta0_hat] = per_time_complex_visibility(nPlus, Nshots, theta_K);
            snr  = Vhat./max(Vse,eps);
            keep = snr >= SNR_min;
            if sum(keep) < 3
                [~,ord] = sort(snr,'descend'); keep = false(size(snr));
                keep(ord(1:min(3,numel(ord)))) = true;
            end

            [chi0, phi00] = fit_phase_vs_time(theta0_hat(keep), tset(keep), Vse(keep));
            C0  = median(Vhat(keep));
            G0  = max(0, robust_decay_init(Vhat(keep), tset(keep)));

            x0 = [logit(C0), phi00, log(max(G0,1e-6)), log(max(chi0,1e-8))];
            [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat, okH] = ...
                fit_fringe_MLE_robust(theta_K, tset, Nshots, nPlus, x0);

            badSE = any(~isfinite(SE)) ...
                 || SE(1) > 1 ...
                 || SE(3) > 10*max(Gamma_hat,1e-12) ...
                 || SE(4) > 10*max(chi_hat,1e-12);
            if ~okH || badSE
                Pse = parametric_SEs(300, theta_K, tset, Nshots, C_hat, phi_hat, Gamma_hat, chi_hat);
                SE = Pse; Cov_nat = diag(Pse.^2);
            end
            SE = max(SE, [0.02; 1e-3; 0.20*max(Gamma_hat,1e-12); 0.20*max(chi_hat,1e-12)]);

            % curves + bands
            tt = linspace(min(tset), max(tset), 120);
            Vfit = C_hat*exp(-2*Gamma_hat*tt);
            band = visibility_band(tt, C_hat, Gamma_hat, Cov_nat([1 3],[1 3]));
            [Vlow,Vhigh] = bootstrap_band(B_boot, theta_K, tset, Nshots, nPlus, ...
                                          @(nk,Nk) per_time_complex_visibility(nk,Nk,theta_K), 0.16);

            nexttile(tl); hold on; box on;
            set(gca,'TickLabelInterpreter','latex');
            fill([tt, fliplr(tt)],[band.low, fliplr(band.high)], [0.8 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.35);
            fill([tset, fliplr(tset)],[Vlow, fliplr(Vhigh)], [0.6 0.8 1.0], 'EdgeColor','none','FaceAlpha',0.35);
            plot(tt, Vfit, 'k-');
            errorbar(tset(keep), Vhat(keep), Vse(keep), 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor',[0.9 0.4 0.1], 'Color',[0.9 0.4 0.1], ...
                'CapSize',0, 'LineStyle','none');
            scatter(tset(~keep), Vhat(~keep), 18, 'filled', 'MarkerFaceColor',[0.7 0.7 0.7]);
            xlabel('t (s)'); ylabel('Visibility'); ylim([0,1]); xlim([min(tset)*0.9, max(tset)*1.05]);
            title(P(u).name);

            % annotation
            txt = sprintf( ...
                ['$C=%.2f\\,\\pm\\,%.2f,\\; \\Gamma=%.2e\\,\\pm\\,%.2e\\;\\mathrm{s}^{-1}$\n' ...
                 '$\\chi=%.2e\\,\\pm\\,%.2e\\;\\mathrm{s}^{-1}$'], ...
                C_hat, SE(1), Gamma_hat, SE(3), chi_hat, SE(4));
            text(0.05, 0.08, txt, 'Units','normalized', ...
                'VerticalAlignment','bottom','FontSize',9);

            % budget collect
            budget_names(end+1,1) = string(P(u).name);
            gamma_com_v(end+1,1)  = gamma_com;
            gamma_com_se(end+1,1) = 0.05*gamma_com;
            gamma_loc_v(end+1,1)  = max(Gamma_hat - gamma_com, 0);
            gamma_loc_se(end+1,1) = SE(3);
        end

        % ---------- Dephasing budget ----------
        fig_budget = figure('Name','Dephasing budget with 68% error bars'+suffix,'Color','w');
        cats = categorical(budget_names);
        barh(cats, [gamma_com_v, gamma_loc_v], 'stacked'); hold on; box on;
        set(gca,'TickLabelInterpreter','latex');
        legend({'$\gamma_{\mathrm{com}}$','$\gamma_{\mathrm{loc}}$'},'Location','southeast');
        xlabel('rate (s$^{-1}$)'); title('Posterior dephasing budget (68\% error bars)');
        yt = 1:numel(budget_names);
        x_end = gamma_com_v + gamma_loc_v;
        for i=1:numel(budget_names)
            errorbar(x_end(i), yt(i), gamma_loc_se(i), 'horizontal', 'k', 'LineStyle','none','CapSize',8);
        end
        if use_log_budget
            set(gca,'XScale','log');
            xmin = max(min([gamma_com_v; x_end]/2), 1e-12);
            xmax = max(x_end + max(gamma_loc_se,1e-18));
            xlim([xmin, 10^(ceil(log10(xmax)))]);
        else
            xmax = max(x_end + gamma_loc_se, [], 'all'); if xmax<=0, xmax=1; end
            xlim([0, 1.2*xmax]);
        end

        % ---------- Total decoherence ----------
        gamma_tot    = gamma_com_v + gamma_loc_v;
        gamma_tot_se = sqrt( (gamma_com_se).^2 + (gamma_loc_se).^2 );
        fig_total = figure('Name','Total decoherence rate'+suffix,'Color','w');
        barh(cats, gamma_tot, 'FaceColor',[0.8 0.45 0.2]); hold on; box on;
        set(gca,'TickLabelInterpreter','latex');
        for i=1:numel(budget_names)
            errorbar(gamma_tot(i), i, gamma_tot_se(i), 'horizontal', 'k', 'LineStyle','none','CapSize',8);
        end
        xlabel('rate (s$^{-1}$)');
        title('Total decoherence $\gamma_{\mathrm{tot}}=\gamma_{\mathrm{com}}+\gamma_{\mathrm{loc}}$');
        if use_log_budget
            set(gca,'XScale','log');
            xmin = max(min(gamma_tot/2), 1e-12);
            xmax = max(gamma_tot + gamma_tot_se);
            xlim([xmin, 10^(ceil(log10(xmax)))]);
        end

        % ---------- Scaling law ----------
        fig_scal = scaling_check(G,hbar);
        set(fig_scal,'Name','Scaling check'+suffix);
    end
end

% ============================= Helpers ==================================
function [C_hat, phi_hat, Gamma_hat, chi_hat, SE, Cov_nat, okH] = fit_fringe_MLE_robust(theta, t, Nshots, nPlus, x0)
    options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5e4,'MaxIter',5e3);
    nll = @(x) fringe_nll(x, theta, t, Nshots, nPlus);
    xhat = fminsearch(nll, x0, options);
    [C_hat, phi_hat, Gamma_hat, chi_hat] = decode_params(xhat);

    H = num_hessian(nll, xhat) + 1e-10*eye(4);
    s = svd(H); condH = s(1)/max(s(end),eps);
    okH = isfinite(condH) && condH < 1e8;

    Cov = pinv(H);
    J = diag([C_hat*(1-C_hat), 1, Gamma_hat, chi_hat]); % d[natural]/d[x]
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
    [Q,Rq] = qr(X,0);
    beta = Rq\(Q'*y);
    b0=beta(1); a=beta(2); b=beta(3); c=beta(4);
    resid = y - X*beta;
    s2 = sum(resid.^2)/(numel(y)-size(X,2));
    Cov = s2 * inv(Rq'*Rq);
    se  = sqrt(max(diag(Cov),0));
    CI  = 1.96 * se(2:4);

    fig = figure('Name','Scaling check','Color','w');
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(lm, y, 6, 'filled'); xx = sort(lm);
    plot(xx, b0 + a*(xx-mean(lm)), 'r-');
    xlabel('$\log_{10}(m/\mathrm{kg})$');
    ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('a=%.2f \\pm %.2f (exp 2.0)', a, CI(1)));

    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(ldx, y, 6, 'filled'); xx = sort(ldx);
    plot(xx, b0 + b*(xx-mean(ldx)), 'r-');
    xlabel('$\log_{10}(\Delta x/\mathrm{m})$');
    ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('b=%.2f \\pm %.2f (exp 2.0)', b, CI(2)));

    nexttile; hold on; box on; set(gca,'TickLabelInterpreter','latex');
    scatter(lr, y, 6, 'filled'); xx = sort(lr);
    plot(xx, b0 - c*(xx-mean(lr)), 'r-');
    xlabel('$\log_{10}(r/\mathrm{m})$');
    ylabel('$\log_{10}(\chi/\mathrm{s}^{-1})$');
    title(sprintf('c=%.2f \\pm %.2f (exp 3.0)', c, CI(3)));
end

function s = sanitize_filename(t)
    s = regexprep(char(t),'[^A-Za-z0-9_\- ]','');
    s = strrep(strtrim(s),' ','_');
    if isempty(s), s='figure'; end
end
