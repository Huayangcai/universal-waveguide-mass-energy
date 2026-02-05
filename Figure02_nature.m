function Figure02_nature()
%% ============================================================
% Figure (2x3): Impact of feedback magnitude kappa and phase Theta on
% absorption partition, useful load power, internal dissipation, and efficiency
%
% S7-consistent definitions:
%   F = kappa * exp(i Theta)
%   Gamma_g = (Gamma_S + Gamma_L E)/(1 + Gamma_S Gamma_L E) with E = e^{-2K}
%   A = 1 - |Gamma_g|^2
%   eta_use = (1-|Gamma_L|^2)e^{-2R} / |1+F|^2   (so T_L = A * eta_use)
%   T_int = A - T_L
%   eta_use_total = T_L/A (A>0), else 0
%
% Passivity feasibility: |Gamma_L| <= 1  =>  kappa <= kappa_max = |Gamma_S| e^{-2R}
% ============================================================

clear; close all; clc;

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',2.0);

%% ==================== PARAMETERS ====================
delta_R = 0.1;
delta_G = 0.1;
Omega   = 1.0;

K0 = sqrt((delta_R + 1i*Omega) * (delta_G + 1i*Omega));
R0 = real(K0);

Eabs = exp(-2*R0); % |E| = e^{-2R}

Gamma_S_mag = 0.70;
phi_S       = 0.0;
GammaS      = Gamma_S_mag * exp(1i*phi_S);

kappa_max = Gamma_S_mag * Eabs;
kappa_vec = linspace(0.0, 1.10*kappa_max, 280);
Theta_vec = linspace(0, 2*pi, 241);
[KAPPA, THETA] = meshgrid(kappa_vec, Theta_vec);

GammaL_mag = (KAPPA .* exp(2*R0)) / Gamma_S_mag;
feasible   = (GammaL_mag <= 1);

F = KAPPA .* exp(1i*THETA);

%% ==================== CORE MAPS ====================
Gamma_g = (GammaS + F./GammaS) ./ (1 + F);

A = 1 - abs(Gamma_g).^2;
A = max(A, 0);

eta_use = ((1 - GammaL_mag.^2) * Eabs) ./ abs(1 + F).^2;
eta_use = max(min(real(eta_use),1),0);

T_L   = A .* eta_use;
T_int = A - T_L;

A(~feasible)       = NaN;
eta_use(~feasible) = NaN;
T_L(~feasible)     = NaN;
T_int(~feasible)   = NaN;

%% ==================== KEY POINTS (analytic, S7) ====================
Theta_res = pi;
Theta_ant = 0;

kappa_opt = (Gamma_S_mag^2) * exp(-4*R0);
kappa_cc  = (Gamma_S_mag^2);
GammaL_cc = Gamma_S_mag * exp(2*R0);
cc_feasible = (GammaL_cc <= 1);

%% ==================== FIGURE LAYOUT (2x3) ====================
width_cm  = 24.0;
height_cm = 18;
fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% ----- Panel (a): T_L -----
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
shade_infeasible(ax1, kappa_max, kappa_vec);
h1 = imagesc(ax1, kappa_vec, Theta_vec, T_L); set(ax1,'YDir','normal'); h1.AlphaData = ~isnan(T_L);
colormap(ax1, safe_turbo());
cb1 = colorbar(ax1); cb1.Label.String = '$\mathcal{T}_L$'; cb1.Label.Interpreter='latex'; cb1.TickLabelInterpreter='none';
caxis(ax1, nanclim(T_L, 2, 98));
format_colorbar_2dp(cb1);
xlabel(ax1,'$\kappa$'); ylabel(ax1,'$\Theta$ (rad)');
% title(ax1,'Useful load power fraction $\mathcal{T}_L(\kappa,\Theta)$');
add_overlays(ax1, kappa_max, kappa_opt, kappa_cc, cc_feasible);
% add_panel_label(ax1,'a');

% ----- Panel (b): A -----
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
shade_infeasible(ax2, kappa_max, kappa_vec);
h2 = imagesc(ax2, kappa_vec, Theta_vec, A); set(ax2,'YDir','normal'); h2.AlphaData = ~isnan(A);
colormap(ax2, parula(256));
cb2 = colorbar(ax2); cb2.Label.String = '$\mathcal{A}=1-|\bar{\Gamma}_g|^2$'; cb2.Label.Interpreter='latex'; cb2.TickLabelInterpreter='none';
caxis(ax2, [0 1]);
format_colorbar_2dp(cb2);
xlabel(ax2,'$\kappa$'); ylabel(ax2,'$\Theta$ (rad)');
% title(ax2,'One-port absorptivity $\mathcal{A}(\kappa,\Theta)$');
add_overlays(ax2, kappa_max, kappa_opt, kappa_cc, cc_feasible);
% add_panel_label(ax2,'b');

% ----- Panel (c): eta_use -----
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
shade_infeasible(ax3, kappa_max, kappa_vec);
h3 = imagesc(ax3, kappa_vec, Theta_vec, eta_use); set(ax3,'YDir','normal'); h3.AlphaData = ~isnan(eta_use);
colormap(ax3, safe_turbo());
cb3 = colorbar(ax3); cb3.Label.String = '$\eta_{\mathrm{use}}=\mathcal{T}_L/\mathcal{A}$'; cb3.Label.Interpreter='latex'; cb3.TickLabelInterpreter='none';
caxis(ax3, [0 1]);
format_colorbar_2dp(cb3);
xlabel(ax3,'$\kappa$'); ylabel(ax3,'$\Theta$ (rad)');
% title(ax3,'Useful-power efficiency $\eta_{\mathrm{use}}(\kappa,\Theta)$');
add_overlays(ax3, kappa_max, kappa_opt, kappa_cc, cc_feasible);
% add_panel_label(ax3,'c');

% ----- Panel (d): T_int -----
ax4 = nexttile(tlo,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
shade_infeasible(ax4, kappa_max, kappa_vec);
h4 = imagesc(ax4, kappa_vec, Theta_vec, T_int); set(ax4,'YDir','normal'); h4.AlphaData = ~isnan(T_int);
colormap(ax4, parula(256));
cb4 = colorbar(ax4); cb4.Label.String = '$\mathcal{T}_{\mathrm{int}}=\mathcal{A}-\mathcal{T}_L$'; cb4.Label.Interpreter='latex'; cb4.TickLabelInterpreter='none';
caxis(ax4, nanclim(T_int, 2, 98));
format_colorbar_2dp(cb4);
xlabel(ax4,'$\kappa$'); ylabel(ax4,'$\Theta$ (rad)');
% title(ax4,'Internal dissipation share $\mathcal{T}_{\mathrm{int}}(\kappa,\Theta)$');
add_overlays(ax4, kappa_max, kappa_opt, kappa_cc, cc_feasible);
% add_panel_label(ax4,'d');

% ----- Panel (e): Resonant cut Theta=pi -----
ax5 = nexttile(tlo,5); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');
[~, jpi] = min(abs(Theta_vec - pi));
TL_pi   = T_L(jpi,:);
A_pi    = A(jpi,:);
Tint_pi = T_int(jpi,:);
eta_pi  = eta_use(jpi,:);

plot(ax5, kappa_vec, A_pi,    '-',  'DisplayName','$\mathcal{A}$');
plot(ax5, kappa_vec, TL_pi,   '-',  'DisplayName','$\mathcal{T}_L$');
plot(ax5, kappa_vec, Tint_pi, '--', 'DisplayName','$\mathcal{T}_{\mathrm{int}}$');
plot(ax5, kappa_vec, eta_pi,  ':',  'DisplayName','$\eta_{\mathrm{use}}$');

xline(ax5, kappa_max, 'k-', 'LineWidth',1.0, 'DisplayName','$\kappa_{\max}$');
xline(ax5, kappa_opt, 'k:', 'LineWidth',1.2, 'DisplayName','$\kappa_{\mathrm{opt}}$');
if cc_feasible
    xline(ax5, kappa_cc, 'r:', 'LineWidth',1.2, 'DisplayName','$\kappa_{\mathrm{cc}}$');
else
    xline(ax5, kappa_cc, 'r:', 'LineWidth',1.2, 'DisplayName','$\kappa_{\mathrm{cc}}$ (inadmissible)');
end

xlabel(ax5,'$\kappa$'); ylabel(ax5,'fraction');
% title(ax5,'Resonant cut ($\Theta=\pi$): absorption partition and efficiency');
legend(ax5,'Location','eastoutside','Interpreter','latex');
legend boxoff
xlim(ax5,[min(kappa_vec) max(kappa_vec)]);
ylim(ax5,[0 1]);
% add_panel_label(ax5,'e');

% ----- Panel (f): Phase cut at kappa = kappa_opt -----
ax6 = nexttile(tlo,6); hold(ax6,'on'); box(ax6,'on'); grid(ax6,'on');
[~, iopt] = min(abs(kappa_vec - kappa_opt));
TL_k   = T_L(:,iopt);
A_k    = A(:,iopt);
eta_k  = eta_use(:,iopt);
Tint_k = T_int(:,iopt);

plot(ax6, Theta_vec, A_k,   '-',  'DisplayName','$\mathcal{A}$');
plot(ax6, Theta_vec, TL_k,  '-',  'DisplayName','$\mathcal{T}_L$');
plot(ax6, Theta_vec, Tint_k,'--', 'DisplayName','$\mathcal{T}_{\mathrm{int}}$');
plot(ax6, Theta_vec, eta_k, ':',  'DisplayName','$\eta_{\mathrm{use}}$');
xline(ax6, Theta_res, 'k--','LineWidth',1.1,'DisplayName','$\Theta=\pi$');
xline(ax6, Theta_ant, 'k:','LineWidth',1.1,'DisplayName','$\Theta=0$');

xlabel(ax6,'$\Theta$ (rad)'); ylabel(ax6,'fraction');
% title(ax6,'Phase cut at $\kappa=\kappa_{\mathrm{opt}}$');
legend(ax6,'Location','eastoutside','Interpreter','latex');
legend boxoff
xlim(ax6,[0 2*pi]);
xticks(ax6,[0 pi 2*pi]);
xticklabels(ax6,{'0','$\pi$','$2\pi$'});
ylim(ax6,[0 1]);

text(ax1, -0.25, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax2, -0.25, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax3, -0.25, 0.98, 'c', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax4, -0.25, 0.98, 'd', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');
text(ax5, -0.25, 0.98, 'e', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax6, -0.255, 0.98, 'f', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

set_theta_ticks(ax1);
set_theta_ticks(ax2);
set_theta_ticks(ax3);
set_theta_ticks(ax4);

text(ax2,0.2,3.5,'$\kappa_{opt}$','Interpreter','latex')
text(ax2,0.4,3.5,'$\kappa_{cc}$','Interpreter','latex')

% ----- Figure title -----
% title(tlo, sprintf(['Feedback control of absorbed-power partition: ' ...
%     '$\\Omega=%.2f$, $\\delta_R=%.2f$, $\\delta_G=%.2f$, $|\\Gamma_S|=%.2f$, $R=%.3f$; ' ...
%     '$\\kappa_{\\max}=|\\Gamma_S|e^{-2R}=%.3f$'], ...
%     Omega, delta_R, delta_G, Gamma_S_mag, R0, kappa_max), 'Interpreter','latex');

%% ---------------- Save ----------------
print(fig, '-dpng', '-r600', 'Figure_F_kappa_Theta_impact_2x3.png');
fprintf('Saved: Figure_F_kappa_Theta_impact_2x3.png\n');

%% ---------------- Console summary ----------------
fprintf('R = %.6f, e^{-2R}=%.6f, kappa_max=%.6f\n', R0, Eabs, kappa_max);
fprintf('kappa_opt=%.6f (Theta=pi), kappa_cc=%.6f (%s)\n', kappa_opt, kappa_cc, ternary(cc_feasible,'feasible','inadmissible'));

end

%% ==================== Helper functions ====================

function add_overlays(ax, kappa_max, kappa_opt, kappa_cc, cc_feasible)
hold(ax,'on');
xline(ax, kappa_max, 'k-', 'LineWidth',1.0);
yline(ax, pi,       'k--','LineWidth',1.0);
yline(ax, 0,        'k:','LineWidth',1.0);
plot(ax, kappa_opt, pi, 'ko', 'MarkerSize',7, 'MarkerFaceColor','w');
if cc_feasible
    plot(ax, kappa_cc, pi, 'rs', 'MarkerSize',7, 'MarkerFaceColor','r');
else
    plot(ax, kappa_cc, pi, 'rx', 'MarkerSize',9, 'LineWidth',2);
end
end

function shade_infeasible(ax, kappa_max, kappa_vec)
x0 = kappa_max;
x1 = max(kappa_vec);
patch(ax, [x0 x1 x1 x0], [0 0 2*pi 2*pi], [0.92 0.92 0.92], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
end

function clim = nanclim(X, pLow, pHigh)
v = X(:);
v = v(~isnan(v));
if isempty(v)
    clim = [0 1];
    return;
end
v = sort(v);
n = numel(v);
i1 = max(1, round((pLow/100)*n));
i2 = max(1, round((pHigh/100)*n));
i1 = min(i1,n); i2 = min(i2,n);
clim = [v(i1) v(i2)];
if clim(1) == clim(2)
    clim = [min(v) max(v)];
    if clim(1)==clim(2)
        clim = [clim(1) clim(1)+eps];
    end
end
end

function format_colorbar_2dp(cb)
drawnow; % ensure ticks are realized
ticks = cb.Ticks;
if isempty(ticks)
    return;
end
cb.TickLabels = compose('%.2f', ticks);
cb.TickLabelInterpreter = 'none';
end

function cmap = safe_turbo()
try
    cmap = turbo(256);
catch
    cmap = parula(256);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function set_theta_ticks(ax)
ylim(ax,[0 2*pi]);
yticks(ax,[0 pi 2*pi]);
yticklabels(ax,{'$0$','$\pi$','$2\pi$'});
end
