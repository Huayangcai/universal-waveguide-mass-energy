function Figure01_nature()
% -------------------------------------------------------------
% Figure 1: Generalized waveguide mass--energy relation and Cai--Smith representation with intrinsic asymmetry.
% (b) Cai-Smith disk: Gamma_g(Omega) Möbius map + xi radial deformation
% (a) U^2 = S^2 + |Gamma_g(zeta)|^2 with high-speed approximation
% -------------------------------------------------------------
clear; clc; close all;

%% ---------- Global style ----------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',2.0);

%% =========================================================
% Shared system parameters (used by both panels)
%% =========================================================
xi_vals = [0.1, 0.5, 0.9];
colors = [0.80 0.20 0.20;
          0.20 0.20 0.60;
          0.20 0.60 0.20];

% Caption-like system (panel b)
Gamma_S = 0.2;
Gamma_L = 0.6;
Phi = linspace(-pi/2, pi/2, 260); % phase sweep
R0 = 0.10; % attenuation scale for K
K = R0 + 1i*Phi;
E = exp(-2*K);
Gamma_bar = (Gamma_L.*E + Gamma_S) ./ (1 + Gamma_S*Gamma_L.*E);

% Option: derive Gamma_g0 for panel (a) from panel (b)
useDerivedGamma_g0 = true; % <-- set false to keep Gamma_g0 = 0.7
Gamma_g0_fixed = 0.7;
if useDerivedGamma_g0
    Gamma_g0 = median(abs(Gamma_bar),'omitnan');
else
    Gamma_g0 = Gamma_g0_fixed;
end

R = 0.2; % accumulated attenuation used in panel (a) caption text

%% =========================================================
% Create figure and tiled layout
%% =========================================================
fig = figure('Color','w','Units','centimeters','Position',[2 0 19.5 12.5]);
tlo = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

%% =========================================================
% PANEL (a): Mass-energy shells
%% =========================================================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

% |Gamma_g(zeta)| = |Gamma_g(Omega)| * exp{-2(2xi-1)Re(zeta)}
Gamma_abs = @(xi) abs(Gamma_g0) .* exp(-2*(2*xi-1)*R);

S = linspace(-2,2,600);

% ---- ADD: black dashed line U = |S| ----
plot(ax1, S, abs(S), 'k--', 'LineWidth', 1.8, ...
    'DisplayName', '$\mathcal{U}=|\mathcal{S}|$');

Smin = 0.06; % avoid the |S| in denominator for approximation
mask = abs(S) >= Smin;

for k = 1:numel(xi_vals)
    xi = xi_vals(k);
    m = Gamma_abs(xi); % mass term = |Gamma_g(zeta)|
    U_exact  = sqrt(S.^2 + m.^2);
    U_approx = nan(size(S));
    U_approx(mask) = abs(S(mask)) + (m.^2)./(2*abs(S(mask)));

    plot(ax1, S, U_exact, '-', 'Color', colors(k,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('$\\xi=%.1f$, $|\\Gamma_g|=%.2f$', xi, m));
    plot(ax1, S, U_approx, '--', 'Color', colors(k,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

xlabel(ax1,'$\mathcal{S}$','FontSize',14);
ylabel(ax1,'$\mathcal{U}$','FontSize',14);
% title(ax1,'(a) State-dependent mass-energy relation','FontSize',14,'FontWeight','bold', ...
%     'Position', [0, 2.5, 0]);
axis(ax1,[-2 2 0 2.5]);

legend1 = legend(ax1,'Location','southoutside','Box','off', ...
    'Orientation','horizontal','NumColumns',2);

%% =========================================================
% PANEL (b): Cai-Smith disk
%% =========================================================
ax2 = nexttile(tlo,2); hold(ax2,'on'); axis(ax2,'equal');

% Background scalar field: |Gamma_g| on unit disk
N = 450;
x = linspace(-1,1,N); y = linspace(-1,1,N);
[X,Y] = meshgrid(x,y);
G = X + 1i*Y;
magG = abs(G);
magG(magG>1) = NaN;
imagesc(ax2,x,y,magG,'AlphaData',~isnan(magG));
set(ax2,'YDir','normal');

% Colormap (fallback if turbo unavailable)
try
    colormap(ax2,turbo(256));
catch
    colormap(ax2,parula(256));
end

% Unit circle + constant-|Gamma| contours
th = linspace(0,2*pi,800);
plot(ax2,cos(th),sin(th),'k-','LineWidth',2.5,'HandleVisibility','off');
r_circles = [0.3 0.6 0.9];
for r = r_circles
    plot(ax2,r*cos(th),r*sin(th),'w--','LineWidth',1.2,'HandleVisibility','off');
end

% xi deformation: Gamma_g(zeta) = Gamma_g(Omega)*exp{-2(2xi-1)*R0}
plot_traj = @(xi,clr,name) ...
    plot(ax2, real( abs(Gamma_bar).*exp(-2*(2*xi-1)*R0).*exp(1i*angle(Gamma_bar)) ), ...
              imag( abs(Gamma_bar).*exp(-2*(2*xi-1)*R0).*exp(1i*angle(Gamma_bar)) ), ...
         '-', 'Color', clr, 'LineWidth', 2.5, 'DisplayName', name);

plot_traj(0.1,[0.8 0.2 0.2],'$\xi=0.1$');
plot_traj(0.5,[0.2 0.2 0.6],'$\xi=0.5$');
plot_traj(0.9,[0.2 0.6 0.2],'$\xi=0.9$');

% ---- Fix limits first (important for stable arrows/labels) ----
xlim(ax2,[-1.1 1.1]);
ylim(ax2,[-1.1 1.1]);
axis(ax2,'equal');

% ===== Coordinate arrows (keep) + axis off (no axes shown) =====
arrow_len = 1.0;
quiver(ax2, 0, 0, arrow_len, 0, 'k', 'LineWidth', 1.5, ...
    'MaxHeadSize', 0.3, 'AutoScale','off', 'HandleVisibility','off');
quiver(ax2, 0, 0, 0, arrow_len, 'k', 'LineWidth', 1.5, ...
    'MaxHeadSize', 0.3, 'AutoScale','off', 'HandleVisibility','off');

text(ax2, arrow_len+0.05, 0, '$\Re({\Gamma}_g)$', 'FontSize', 14, ...
    'HorizontalAlignment','left','VerticalAlignment','middle');
text(ax2, 0.05, arrow_len+0.00, '$\Im({\Gamma}_g)$', 'FontSize', 14, ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');

% Right panel: NO axes
axis(ax2,'off');
% title(ax2, '(b) Cai-Smith chart with asymmetry', 'FontSize', 14, 'FontWeight', 'bold', ...
%     'Position', [0, 1.23, 0]);

cb = colorbar(ax2,'southoutside');
cb.Label.String = 'Effective mass $|\Gamma_g|$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

%% ---------- Add (a)(b)labels at each tile top-left (robust)
% Nature-style panel labels (bold, stable, not affected by axis limits)
text(ax1, -0.2, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax2, 0.02, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

%% ---------- Export (stable) ----------
print(fig,'Figure01','-dpng','-r600');
fprintf('Saved: Figure01.png\n');
end

%% ===================== Helper: panel label using annotation =================
function h = add_panel_label(ax, str, dy_pt)
% Place panel label just ABOVE the top-left corner of the AXES FRAME (outside).
% Nature style: close to the corner, not centered.
%
% h = add_panel_label(ax,'a')
% h = add_panel_label(ax,'(a)', 3)   % dy in points (upward)

    if nargin < 3 || isempty(dy_pt), dy_pt = 3; end

    if ~isgraphics(ax,'axes')
        error('add_panel_label: ax must be a valid axes handle.');
    end
    fig = ancestor(ax,'figure');

    drawnow; % ensure TightInset is accurate (important with tiledlayout)

    % --- axes geometry in normalized units ---
    oldUnits = ax.Units;
    ax.Units = 'normalized';
    pos = ax.Position;
    ti  = ax.TightInset;      % [left bottom right top]
    ax.Units = oldUnits;

    % --- plot box (actual framed drawing region) ---
    px = pos(1) + ti(1);
    py = pos(2) + ti(2);
    pw = pos(3) - ti(1) - ti(3);
    ph = pos(4) - ti(2) - ti(4);

    % --- convert point offset -> normalized offset ---
    fig_old = fig.Units;
    fig.Units = 'points';
    figPos_pt = fig.Position;           % [x y w h] in points
    fig.Units = fig_old;

    dy = dy_pt / figPos_pt(4);

    % --- anchor: directly above the top-left corner of the plot box ---
    x = px-0.15;               % align to left edge
    y = py + ph + dy;     % just above top edge (outside)

    % textbox size (normalized); can be small because we fit to text
    wBox = 0.04; 
    hBox = 0.04;

    h = annotation(fig,'textbox',[x y wBox hBox], ...
        'String',str, ...
        'LineStyle','none', ...
        'FontName','Helvetica', ...
        'FontSize',14, ...
        'FontWeight','bold', ...
        'Interpreter','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FitBoxToText','on', ...
        'Margin',1);
end
