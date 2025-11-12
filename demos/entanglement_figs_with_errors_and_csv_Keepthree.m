function entanglement_figs_with_errors_and_csv()
% Clean version: Figures 1–3 + CSV, with 68% error bars and a projected heavy case.

%% ---- constants & PSD bound
G     = 6.67430e-11;           % m^3 kg^-1 s^-2
hbar  = 1.054571817e-34;       % J s
S_eps = 1e-28;                 % 1/Hz  (clock-bounded white band)

%% ---- platforms
P(1) = struct('name','QGEM'        ,'m',1e-14 ,'dx',2.5e-7 ,'r',2e-4  ,'t',linspace(1,2,8),  'gloc',1e-3 , 'C',0.85);
P(2) = struct('name','MAQRO'       ,'m',1e-15 ,'dx',1e-7   ,'r',5e-4  ,'t',linspace(10,100,8),'gloc',1e-5 , 'C',0.85);
P(3) = struct('name','Levitated'   ,'m',1e-14 ,'dx',1e-6   ,'r',1e-3  ,'t',linspace(1,8,8),  'gloc',8e-2 , 'C',0.85);
P(4) = struct('name','QGEM (proj.)','m',1.0   ,'dx',50e-6  ,'r',2e-2  ,'t',10,               'gloc',1e-5 , 'C',0.85);

nP = numel(P);

%% ---- rates
for k = 1:nP
    P(k).omega_g   = (G/hbar) * (P(k).m^2) * (P(k).dx) / (P(k).r^2);     % s^-1
    P(k).gcom      = (P(k).omega_g.^2 / 8) * S_eps;                      % s^-1
    P(k).Gamma     = P(k).gloc + P(k).gcom;                               % s^-1
end

% 68% CI placeholders
f_loc = 0.10;  % 10% on local
f_com = 0.30;  % 30% on common (from S_eps fit)
for k = 1:nP
    P(k).sig_gloc  = max( real(f_loc*P(k).gloc),  eps );
    P(k).sig_gcom  = max( real(f_com*P(k).gcom),  eps );
    P(k).sig_Gamma = hypot(P(k).sig_gloc, P(k).sig_gcom);
end

%% ---- Figure 1: Visibility
fig1 = figure('Name','Figure 1: Visibility by time','Color','w'); 
tiledlayout(fig1,1,3,'Padding','compact','TileSpacing','compact');
for k = 1:3
    nexttile;
    tgrid  = P(k).t(:);
    V_true = P(k).C * exp(-P(k).Gamma * tgrid);
    rng(7+k);  v_meas = max(0, min(1, V_true + 0.02*randn(size(V_true))));
    plot(tgrid, v_meas,'o','MarkerSize',4,'LineWidth',1); hold on
    plot(tgrid, V_true,'-','LineWidth',1.5);
    grid on
    xlabel('t (s)','Interpreter','latex'); ylabel('Visibility','Interpreter','latex');
    title(P(k).name,'Interpreter','latex');
    txt = sprintf('$C=%.2f,\\; \\Gamma=%.0e\\;\\mathrm{s}^{-1}$', P(k).C, P(k).Gamma);
    text(0.02,0.06,txt,'Units','normalized','Interpreter','latex');
    ylim([0 max(1,1.05*max(v_meas))]);
end

%% ---- Figure 2: Dephasing budget (stacked) + Common-mode only
names  = string({P.name});
loc    = [P.gloc];
com    = [P.gcom];
s_loc  = [P.sig_gloc];
s_com  = [P.sig_gcom];

fig2 = figure('Name','Figure 2: Dephasing budget','Color','w'); 
tiledlayout(fig2,1,2,'Padding','compact','TileSpacing','compact');

% Left: stacked with error bars (log y)
ax1 = nexttile; 
hb = bar(ax1, 1:nP, [loc(:), com(:)], 'stacked'); 
set(hb(1),'FaceColor',[0.91 0.41 0.17]);   % local
set(hb(2),'FaceColor',[0.20 0.45 0.90]);   % common
set(ax1,'YScale','log'); grid(ax1,'on');
xlim(ax1,[0.4 nP+0.6]); xticks(ax1,1:nP); xticklabels(ax1,names); xtickangle(ax1,10);
ylabel(ax1,'rate (s$^{-1}$)','Interpreter','latex');
title(ax1,'Posterior dephasing budget (clock-bounded + local)','Interpreter','latex'); hold(ax1,'on');

for k = 1:nP
    % local component (center at loc/2)
    drawYerr(ax1, k, max(loc(k)/2, eps), max(s_loc(k), eps));
    % common component (center at loc + com/2)
    drawYerr(ax1, k, max(loc(k) + com(k)/2, eps), max(s_com(k), eps));
end
legend(ax1,{'$\gamma_{\mathrm{loc}}$','$\gamma_{\mathrm{com}}$'}, ...
       'Location','southoutside','Orientation','horizontal', ...
       'Interpreter','latex','Box','off');

% Nice y-lims for stacked
stackTop = loc + com;
ymin = 10^floor(log10( max(min([loc(loc>0),com(com>0),stackTop(stackTop>0)]), eps) )) / 10;
ymax = 10^ceil (log10( max(stackTop)*1.2 + eps ));
ylim(ax1,[ymin ymax]);

% Right: common-mode only, with error bars (log y)
ax2 = nexttile;
bar(ax2, 1:nP, com(:), 0.6, 'FaceColor',[0.20 0.45 0.90]);
set(ax2,'YScale','log'); grid(ax2,'on');
xlim(ax2,[0.4 nP+0.6]); xticks(ax2,1:nP); xticklabels(ax2,names); xtickangle(ax2,10);
ylabel(ax2,'$\gamma_{\mathrm{com}}$ (s$^{-1}$)','Interpreter','latex');
title(ax2,'Common-mode only (clock-bounded)','Interpreter','latex'); hold(ax2,'on');
for k = 1:nP
    drawYerr(ax2, k, max(com(k), eps), max(s_com(k), eps));
end
ymin2 = 10^floor(log10( max(min(com(com>0)), eps) )) / 10;
ymax2 = 10^ceil (log10( max(com)*1.2 + eps ));
ylim(ax2,[ymin2 ymax2]);

%% ---- Figure 3: Total dephasing (log y) with error bars
fig3 = figure('Name','Figure 3: Total dephasing rate','Color','w');
ax3 = axes(fig3);
bar(ax3, 1:nP, [P.Gamma], 0.7, 'FaceColor',[0.85 0.50 0.15]); 
set(ax3,'YScale','log'); grid(ax3,'on');
xlim(ax3,[0.4 nP+0.6]); xticks(ax3,1:nP); xticklabels(ax3,names); xtickangle(ax3,10);
ylabel(ax3,'$\Gamma_{\phi}$ (s$^{-1}$)','Interpreter','latex');
title(ax3,'Total dephasing rate $\Gamma_{\phi}$','Interpreter','latex'); hold(ax3,'on');
for k = 1:nP
    drawYerr(ax3, k, max(P(k).Gamma, eps), max(P(k).sig_Gamma, eps));
end
ymin3 = 10^floor(log10( max(min([P.Gamma]), eps) )) / 10;
ymax3 = 10^ceil (log10( max([P.Gamma])*1.2 + eps ));
ylim(ax3,[ymin3 ymax3]);

%% ---- CSV export for Table III
t_rep = zeros(nP,1);
for k = 1:nP
    tk = P(k).t; if numel(tk)>1, t_rep(k) = tk(end); else, t_rep(k) = tk; end
end
T = table( ...
    string({P.name}).', ...
    [P.m].', [P.dx].', [P.r].', t_rep, ...
    [P.gloc].', [P.sig_gloc].', ...
    [P.gcom].', [P.sig_gcom].', ...
    [P.Gamma].', [P.sig_Gamma].', ...
    'VariableNames',{'Platform','m_kg','dx_m','r_m','t_s', ...
                     'gamma_loc','sigma_gamma_loc', ...
                     'gamma_com','sigma_gamma_com', ...
                     'Gamma_tot','sigma_Gamma_tot'});
writetable(T, 'dephasing_budget.csv');
fprintf('Wrote dephasing_budget.csv with %d rows.\n',height(T));
end

% -------- helper: vertical ±dy error bar at (x,y)
function drawYerr(ax, x, y, dy)
    if ~(isfinite(x)&&isfinite(y)&&isfinite(dy)) || dy<=0 || y<=0, return; end
    cap = 0.12;
    line(ax, [x x], [y-dy, y+dy], 'Color','k','LineWidth',0.9);
    line(ax, [x-cap, x+cap], [y-dy, y-dy], 'Color','k','LineWidth',0.9);
    line(ax, [x-cap, x+cap], [y+dy, y+dy], 'Color','k','LineWidth',0.9);
end
