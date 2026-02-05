function FigureS4_nature()
% ======================================================================
% Figure S4:Waveguide relativity: velocity composition, rapidity linearisation, and Lorentz boosts in U-S space
% 2x2 visualization aligned with the waveguide mass–energy relation:
%   U^2 - S^2 = |Gamma_g|^2
% and the waveguide "velocity" parameter
%   beta_w = S/U,  |beta_w|<1,  eta = artanh(beta_w)
%
% Output:
%   Figure S4.png  (600 dpi)
% ======================================================================

clear; close all; clc;

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultTextFontSize',14);
set(groot,'DefaultLegendFontSize',14);
set(groot,'DefaultColorbarFontSize',14);

%% -------------------- Figure layout -------------------- 
width_cm = 22.5;
height_cm = 18;
fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

%% ======================================================================
% (a) Classical vs Relativistic velocity addition
% ======================================================================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

beta1 = linspace(-0.99, 0.99, 600);
beta2_list = [0.3, 0.5, 0.9];

C = lines(numel(beta2_list));

for i = 1:numel(beta2_list)
    beta2 = beta2_list(i);

    beta_classical = beta1 + beta2;

    denom = 1 + beta1.*beta2;
    denom = sign(denom).*max(abs(denom), 1e-12);  % numerical guard
    beta_rel = (beta1 + beta2) ./ denom;

    plot(ax1, beta1, beta_classical, '--', 'Color', C(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Classical ($\\beta_{w}^2=%.1f$)', beta2));
    plot(ax1, beta1, beta_rel, '-',  'Color', C(i,:), 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Relativistic ($\\beta_{w}^2=%.1f$)', beta2));
end

% Speed limit lines
plot(ax1, [-1,1], [ 1, 1], 'k:', 'LineWidth', 1.6, 'HandleVisibility','off');
plot(ax1, [-1,1], [-1,-1], 'k:', 'LineWidth', 1.6, 'HandleVisibility','off');

xlim(ax1, [-1, 1]);
ylim(ax1, [-2, 2]);
xlabel(ax1, '$\beta_{w}^{1}$ (first velocity-like parameter)');
ylabel(ax1, '$\beta_{w}^{\mathrm{total}}$');
% title(ax1, 'Velocity addition: classical vs Lorentz');

   text(ax1, -0.12, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

legend1=legend(ax1, 'Location','southeast', 'FontSize', 12, 'Box','off');
% set(legend1,...
%     'Position',[0.285690819405847 0.623888725207911 0.178096071879069 0.145616322755814],...
%     'FontSize',10);

%% ======================================================================
% (b) Rapidity–velocity relation
% ======================================================================
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

eta = linspace(-2.2, 2.2, 600);
beta_eta = tanh(eta);

plot(ax2, eta, beta_eta, '-', 'LineWidth', 2.6);

xlabel(ax2, '$\phi$');
ylabel(ax2, '$\beta_{w}=\tanh(\phi)$');
% title(ax2, 'Rapidity linearizes composition');

% Key formulas (kept compact)
text(ax2,0.5, 0.16, '$\phi_{w}^{\rm total}=\phi_1+\phi_2+\cdots$', 'Units','normalized', ...
    'FontSize', 12, 'FontWeight','bold', 'Color',[0.85 0 0]);
text(ax2, 0.5, 0.08, '$\phi_{w}^{\rm total}=\tanh(\phi_{\rm total})$', 'Units','normalized', ...
    'FontSize', 12, 'FontWeight','bold', 'Color',[0.85 0 0]);

% Mark example points
eta_pts = [0.3, 0.7, 1.2];
for k = 1:numel(eta_pts)
    e = eta_pts(k);
    b = tanh(e);
    plot(ax2, e, b, 'o', 'MarkerSize', 7.5, 'MarkerFaceColor',[0.85 0 0], ...
        'MarkerEdgeColor','k', 'LineWidth', 0.8);
    text(ax2, e-1.7, b+0.08, sprintf('$\\phi=%.1f,\\ \\beta_{w}=%.2f$', e, b), ...
        'FontSize', 12, 'Interpreter','latex');
end

   text(ax2, -0.12, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

ylim(ax2, [-1.05, 1.05]);

%% ======================================================================
% (c) Multi-section waveguide: derive eta_k from |Gamma_g,k| via U=1
%     Using U=1 => S = sqrt(1-|Gamma|^2), beta=S/U = sqrt(1-|Gamma|^2)
% ======================================================================
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');

Gamma_abs_k = [0.999, 0.995, 0.990, 0.985];  % example per-section "mass term"
Gamma_abs_k = min(max(Gamma_abs_k, 1e-12), 1-1e-12);

beta_k = sqrt(1 - Gamma_abs_k.^2);        % with U=1
eta_k  = atanh(beta_k);                    % rapidity per section

cum_eta  = cumsum(eta_k);
cum_beta = tanh(cum_eta);

% Bars: eta_k
b = bar(ax3, 1:numel(eta_k), eta_k, 0.55, 'FaceAlpha', 0.65, ...
    'EdgeColor',[0 0 0], 'LineWidth', 1.0, 'DisplayName','Section rapidity $\phi_k$');

% Line: cumulative eta
plot(ax3, 1:numel(eta_k), cum_eta, '-o', 'LineWidth', 2.2, 'MarkerSize', 7, ...
    'MarkerFaceColor',[0.85 0 0], 'DisplayName','Cumulative $\sum \phi_k$');

xlabel(ax3, 'Waveguide section index');
ylabel(ax3, '$\phi$');

% Labels: show |Gamma| and resulting beta
for i = 1:numel(eta_k)
    txt1 = sprintf('$|\\Gamma_{g,%d}|=%.3f$', i, Gamma_abs_k(i));
    txt2 = sprintf('$\\beta_{w}^{\\rm cum}=%.3f$', cum_beta(i));
    text(ax3, i, eta_k(i)+0.02, txt1, 'HorizontalAlignment','center', ...
        'FontSize', 12, 'Interpreter','latex');
    text(ax3, i, cum_eta(i)+0.05, txt2, 'HorizontalAlignment','center', ...
        'FontSize', 12, 'Color',[0.85 0 0], 'Interpreter','latex');
end

ylim(ax3, [0, max(cum_eta)*1.25]);
legend1=legend(ax3, 'Location','north', 'Box','off', 'FontSize', 12);
set(legend1,...
    'Position',[0.1508522566492 0.405988687710565 0.179346409217984 0.0592320273904239],...
    'FontSize',10);

   text(ax3, -0.12, 0.98, 'c', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

% title(ax3, 'Additivity of rapidity across sections');

%% ======================================================================
% (d) Lorentz boost in (U,S) space: invariant hyperbolae U^2 - S^2 = |Gamma|^2
% ======================================================================
ax4 = nexttile(tlo,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');

eta_boost = 0.6;
L = [cosh(eta_boost), -sinh(eta_boost); -sinh(eta_boost), cosh(eta_boost)];

% Plot light cone U=±S
Sline = linspace(-1.2, 1.2, 200);
plot(ax4, abs(Sline),  Sline, 'k:', 'LineWidth', 1.2, 'HandleVisibility','off'); % U=|S| (both cones)
plot(ax4, abs(Sline), -Sline, 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');

% Invariant hyperbolae for a few |Gamma|
Gamma_levels = [0.25, 0.5, 0.8];
S = linspace(-1.2, 1.2, 600);

for m = Gamma_levels
    U = sqrt(S.^2 + m^2);    % U >= 0 branch
    plot(ax4, U, S, 'g--', 'LineWidth', 2.0, 'HandleVisibility','off'); % original hyperbola

    % Sample points on this hyperbola and transform them
    idx = round(linspace(1, numel(S), 45));
    P  = [U(idx); S(idx)];
    Pp = L * P;

    scatter(ax4, P(1,:),  P(2,:),  22, 'b', 'filled', 'MarkerFaceAlpha', 0.55, ...
        'HandleVisibility','off');
    scatter(ax4, Pp(1,:), Pp(2,:), 22, 'r', 'filled', 'MarkerFaceAlpha', 0.55, ...
        'HandleVisibility','off');
end

% Add legend proxies
hp1 = scatter(ax4, NaN, NaN, 40, 'b', 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','Original states');
hp2 = scatter(ax4, NaN, NaN, 40, 'r', 'filled', 'MarkerFaceAlpha',0.55, 'DisplayName','Boosted states');
hp3 = plot(ax4, NaN, NaN, 'g--', 'LineWidth', 2.0, 'DisplayName','$\mathcal{U}^2-\mathcal{S}^2=|\Gamma_g|^2$');

xlabel(ax4, '$\mathcal{U}$');
ylabel(ax4, '$\mathcal{S}$');
% title(ax4, sprintf('Lorentz boost (hyperbolic rotation), $\\phi=%.2f$', eta_boost));

xlim(ax4, [0, 1.8]);
ylim(ax4, [-1.2, 1.2]);

legend(ax4, 'Location','east', 'Box','off', 'FontSize', 12);

    text(ax4, -0.12, 0.98, 'd', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

% Small note on invariance
% text(ax4, 0.52, 0.10, '$\mathcal{U}^{\prime 2}-\mathcal{S}^{\prime 2}=\mathcal{U}^2-\mathcal{S}^2$', ...
%     'Units','normalized', 'FontSize', 11, 'Interpreter','latex', 'Color',[0 0 0]);

%% ---------------- Save ----------------
% Set ALL font sizes right before output
% set(findall(gcf,'-property','FontSize'),'FontSize',14);
filename = 'Figure S4';
print(fig, '-dpng', '-r600', [filename, '.png']);
disp(['Figure saved as ', filename, '.png']);

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
