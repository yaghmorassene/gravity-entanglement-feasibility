%% entanglement_figs_demo_fix.m
% Self-contained demo: visibility, dephasing budget, total rates
% with robust log-axis scaling so tiny gamma_com values are visible.

clear; clc;

% ---------- Physical constants
G     = 6.67430e-11;
hbar  = 1.054571817e-34;

% ---------- Platforms (QGEM, MAQRO, Levitated)
plat = {'QGEM','MAQRO','Levitated'};
m    = [1e-14,  1e-15,   1e-14];         % kg
dx   = [2.5e-7, 1.0e-7,  5.0e-8];        % m
r    = [2.0e-4, 5.0e-4,  1.0e-3];        % m
% Nominal local dephasing (illustrative)
gamma_loc = [1e-3, 1e-5, 1e-1];          % s^-1

% Clock-bounded white proper-time PSD at low f
S_eps0 = 1e-28;                          % Hz^-1 (conservative bound)

% ---------- Compute gravitational coupling & dephasing
omega_g    = (G./hbar) .* (m.^2) .* (dx ./ (r.^2));  % s^-1
gamma_com  = (omega_g.^2)/8 .* S_eps0;               % s^-1
gamma_tot  = gamma_loc + gamma_com;                  % s^-1

% ---------- Print numbers so you can cross-check
fprintf('\nDephasing rates (s^-1)\n');
for k = 1:numel(plat)
    fprintf('  %-9s | gamma_com = %.3e | gamma_loc = %.3e | Gamma_phi = %.3e\n', ...
        plat{k}, gamma_com(k), gamma_loc(k), gamma_tot(k));
end
fprintf('\n');

% ---------- FIGURE 1: Visibility vs time (one small panel per platform)
C0 = 0.85;  % baseline contrast
tGrid = { linspace(1.0, 2.0, 6), ...     % QGEM (s)
          linspace(10, 100, 7), ...      % MAQRO (s)
          linspace(1.0, 8.0, 6) };       % Levitated (s)

f1 = figure('Name','Figure 1: Visibility by time'); clf;
tiledlayout(f1,1,3,'TileSpacing','compact','Padding','compact');

for k = 1:3
    nexttile;
    t = tGrid{k};
    V = C0 * exp(-gamma_tot(k)*t);     % ideal exponential envelope
    % make simple "measurements"
    y = V + 0.02*(randn(size(V)));     % small shot noise
    y = min(max(y,0),1);

    plot(t, V, 'k-', 'LineWidth', 1.5); hold on;
    scatter(t, y, 25, 'filled');
    xlabel('t (s)'); ylabel('Visibility');
    title(plat{k});
    grid on;
    ylim([0,1]);
    txt = sprintf('C = %.2f,  \\Gamma = %.2g s^{-1}', C0, gamma_tot(k));
    text(0.02, 0.08, txt, 'Units','normalized','FontSize',9);
end

% ---------- FIGURE 2: Dephasing budget + Common-mode only
f2 = figure('Name','Figure 2: Dephasing budget'); clf;
tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');

% Left: stacked bars (gamma_com + gamma_loc)
nexttile; 
Ystack = [gamma_com(:), gamma_loc(:)];  % Nx2
bh = bar(Ystack,'stacked');
set(gca,'YScale','log');
xticklabels(plat);
ylabel('rate (s^{-1})');
title('Posterior dephasing budget (clock-bounded + local)');
legend({'\gamma_{com}','\gamma_{loc}'},'Location','northwest');
grid on; grid minor;
% Add tiny value labels
for k = 1:numel(plat)
    val = sum(Ystack(k,:));
    text(k, 1.15*max(Ystack(k,:)), sprintf('%.1e', val), ...
        'HorizontalAlignment','center','FontSize',8);
end

% Right: common-mode only, with robust log limits so bars are visible
nexttile;
bh2 = bar(gamma_com,'FaceColor',[0.2 0.4 0.8]);
set(gca,'YScale','log'); xticklabels(plat);
ylabel('\gamma_{com} (s^{-1})');
title('Common-mode only (clock-bounded)');
grid on; grid minor;
robustLogYLim(gca, gamma_com);

% ---------- FIGURE 3: Total dephasing rate
f3 = figure('Name','Figure 3: Total dephasing rate'); clf;
bar(gamma_tot,'FaceColor',[0.85 0.4 0.2]);
set(gca,'YScale','log'); xticklabels(plat);
ylabel('\Gamma_\phi (s^{-1})'); title('Total dephasing rate \Gamma_\phi');
grid on; grid minor;
robustLogYLim(gca, gamma_tot);

%% --------- Helper: robust y-limits for log axes
function robustLogYLim(ax, vals)
    vals = vals(:);
    vals = vals(isfinite(vals) & vals>0);
    if isempty(vals)
        % fallback if all zeros/invalid
        ylim(ax,[1e-45, 1e-35]);
        return;
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmin==vmax
        lo = vmin/10;
        hi = vmin*10;
    else
        lo = 10^(floor(log10(vmin))-1);     % one extra order of margin
        hi = 10^(ceil(log10(vmax))+1);
    end
    % keep some sanity bounds
    if ~isfinite(lo) || lo<=0, lo = vmin/10; end
    if ~isfinite(hi) || hi<=lo, hi = vmax*10; end
    ylim(ax,[lo, hi]);
end
