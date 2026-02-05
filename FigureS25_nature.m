% ============================================================
% Extended Data Fig. 25 (UPDATED, Scheme A / flux-based, ASYMMETRIC BOUNDARIES)
% ------------------------------------------------------------
% 2x1 vertical layout:
% (a)  \tilde{j} vs \tilde{\eta} for multiple asymmetry \xi
% (b)  log10(|\tilde{j}|) vs \tilde{\eta} for the same cases
%
% Core (single mode, Scheme A, flux-based):
%   Jox  = alpha * exp(  xi      * lambda * eta )
%   Jred = alpha * exp( -(1-xi)  * lambda * eta )
%   U = Jox + Jred
%   S = Jox - Jred
%   beta = S./U
%   rho = sqrt((1-|beta|)/(1+|beta|)) = exp(-|lambda*eta|/2) (independent of alpha, xi)
%   jtilde(display) = S ./ sqrt(1 + S.^2)
%
% Region boundaries are now computed PER-CURVE and PER-BRANCH:
%   Linear/Transition boundary:   |S| = S_thr1
%   Transition/Saturation boundary:|S| = S_thr2
% giving (eta1_neg, eta1_pos, eta2_neg, eta2_pos) for each \xi.
% These thresholds are overlaid as faint colored vertical lines.
% ============================================================

clear; close all; clc;

%% -------------------- publication style --------------------
set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 1);
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultTextInterpreter', 'latex');

%% -------------------- parameters --------------------
xi_values = [0.2, 0.5, 0.8];
colors    = lines(numel(xi_values));

eta_tilde = linspace(-5, 5, 2000);

% Scheme-A single-mode parameters (dimensionless demo)
alpha  = 1.0;
lambda = 1.0;
maxExp = 70;

% Region thresholds in |S|
S_thr1 = 1;   % |S| ~ 1 : linear -> transition
S_thr2 = 3;   % |S| ~ 3 : transition -> saturation

% y-limits
ylim_a = [-1.1, 1.1];
ylim_b = [-3, 0.1];

%% -------------------- compute curves and per-curve thresholds --------------------
nXi = numel(xi_values);
j_store = zeros(nXi, numel(eta_tilde));
S_store = zeros(nXi, numel(eta_tilde));

% thresholds: [eta1_neg, eta1_pos, eta2_neg, eta2_pos]
thr = nan(nXi, 4);

for i = 1:nXi
    xi = xi_values(i);

    [jtilde, S] = one_mode_schemeA_j(eta_tilde, alpha, xi, lambda, maxExp);
    j_store(i,:) = jtilde;
    S_store(i,:) = S;

    [eta1_neg, eta1_pos] = threshold_eta_from_S(eta_tilde, S, S_thr1);
    [eta2_neg, eta2_pos] = threshold_eta_from_S(eta_tilde, S, S_thr2);

    thr(i,:) = [eta1_neg, eta1_pos, eta2_neg, eta2_pos];
end

fprintf('=== Per-curve asymmetric boundaries (|S|=%.3g and %.3g) ===\n', S_thr1, S_thr2);
for i = 1:nXi
    fprintf('xi=%.2f: eta1=[%.3f, %.3f], eta2=[%.3f, %.3f]\n', ...
        xi_values(i), thr(i,1), thr(i,2), thr(i,3), thr(i,4));
end

%% -------------------- figure layout --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 16]);
tlo = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

%% ============================================================
% (a) \tilde{j} vs \tilde{\eta}
%% ============================================================
ax_a = nexttile(tlo, 1); hold on; grid on; box on;

% Light reference shading (optional): symmetric xi=0.5 boundaries
% (keeps your original "Linear/Transition/Saturation" visual blocks)
xi_ref = 0.5;
[~, Sref] = one_mode_schemeA_j(eta_tilde, alpha, xi_ref, lambda, maxExp);
[eta1r_neg, eta1r_pos] = threshold_eta_from_S(eta_tilde, Sref, S_thr1);
[eta2r_neg, eta2r_pos] = threshold_eta_from_S(eta_tilde, Sref, S_thr2);

patch([eta1r_neg, eta1r_pos, eta1r_pos, eta1r_neg], [ylim_a(1), ylim_a(1), ylim_a(2), ylim_a(2)], ...
      [0.9, 1, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta2r_neg, eta1r_neg, eta1r_neg, eta2r_neg], [ylim_a(1), ylim_a(1), ylim_a(2), ylim_a(2)], ...
      [1, 0.9, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta1r_pos, eta2r_pos, eta2r_pos, eta1r_pos], [ylim_a(1), ylim_a(1), ylim_a(2), ylim_a(2)], ...
      [1, 0.9, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta_tilde(1), eta2r_neg, eta2r_neg, eta_tilde(1)], [ylim_a(1), ylim_a(1), ylim_a(2), ylim_a(2)], ...
      [0.9, 0.9, 1], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta2r_pos, eta_tilde(end), eta_tilde(end), eta2r_pos], [ylim_a(1), ylim_a(1), ylim_a(2), ylim_a(2)], ...
      [0.9, 0.9, 1], 'EdgeColor','none', 'FaceAlpha',0.22);

text(0, 0.45, 'Linear', 'HorizontalAlignment','center', 'FontSize',12, 'BackgroundColor','w');
text(mean([eta2r_neg, eta1r_neg]), 0.45, 'Transition', 'HorizontalAlignment','center', 'FontSize',12, 'BackgroundColor','w');
text(mean([eta1r_pos, eta2r_pos]), 0.45, 'Transition', 'HorizontalAlignment','center', 'FontSize',12, 'BackgroundColor','w');
text(mean([eta_tilde(1), eta2r_neg]), 0.45, 'Saturation', 'HorizontalAlignment','center', 'FontSize',12, 'BackgroundColor','w');
text(mean([eta2r_pos, eta_tilde(end)]), 0.45, 'Saturation', 'HorizontalAlignment','center', 'FontSize',12, 'BackgroundColor','w');

% Main curves
h_curves_a = gobjects(nXi,1);
for i = 1:nXi
    h_curves_a(i) = plot(eta_tilde, j_store(i,:), '-', 'Color', colors(i,:), ...
        'DisplayName', sprintf('$\\xi = %.1f$', xi_values(i)));
end

% Centerline
xline(0, 'k-', 'LineWidth',1, 'HandleVisibility','off');

% Per-curve asymmetric boundary markers (faint, colored)
for i = 1:nXi
    c = 0.75*colors(i,:) + 0.25*[1 1 1]; % lighten
    draw_vline(ax_a, thr(i,1), ylim_a, ':', c, 1.3); % eta1_neg
    draw_vline(ax_a, thr(i,2), ylim_a, ':', c, 1.3); % eta1_pos
    draw_vline(ax_a, thr(i,3), ylim_a, '--', c, 1.0); % eta2_neg
    draw_vline(ax_a, thr(i,4), ylim_a, '--', c, 1.0); % eta2_pos
end

% Asymptotes
h_sat_a = plot(xlim, [1,1], 'r--', 'DisplayName', 'Saturation asymptote: $\tilde{j}=\pm 1$','LineWidth',2.5);
plot(xlim, [-1,-1], 'r--', 'HandleVisibility','off','LineWidth',2.5);

% Linear (first-order) asymptote: S ≈ alpha*lambda*eta, independent of xi
% To keep consistent bounding, show jtilde for that S.
eta_lin = linspace(eta1r_neg, eta1r_pos, 240);
S_lin   = (alpha*lambda) * eta_lin;
j_lin   = S_lin ./ sqrt(1 + S_lin.^2);
h_lin_a = plot(eta_lin, j_lin, 'k--', 'DisplayName', 'Linear asymptote (1st order): $\tilde{j}\approx \alpha\lambda\,\tilde{\eta}$','LineWidth',2.5);

ylabel('$\tilde{j}$');
xlim([eta_tilde(1), eta_tilde(end)]);
ylim(ylim_a);

text(ax_a, -0.10, 0.98, 'a', 'Units','normalized', ...
     'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
     'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

% Legend (keep clean: curves + asymptotes only)
% legend(ax_a, [h_curves_a; h_lin_a; h_sat_a], 'Location','southeast', 'Box','off');

% Small note explaining marker meaning (no legend clutter)
% text(ax_a, 0.02, 0.08, sprintf('Colored vertical lines: |S|=%.0f (:) and |S|=%.0f (--), per-\\xi and per-branch', S_thr1, S_thr2), ...
%     'Units','normalized','FontSize',11,'BackgroundColor','w');

%% ============================================================
% (b) log10(|\tilde{j}|) vs \tilde{\eta}
%% ============================================================
ax_b = nexttile(tlo, 2); hold on; grid on; box on;

eps0 = 10^(ylim_b(1));

% Same reference shading
patch([eta1r_neg, eta1r_pos, eta1r_pos, eta1r_neg], [ylim_b(1), ylim_b(1), ylim_b(2), ylim_b(2)], ...
      [0.9, 1, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta2r_neg, eta1r_neg, eta1r_neg, eta2r_neg], [ylim_b(1), ylim_b(1), ylim_b(2), ylim_b(2)], ...
      [1, 0.9, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta1r_pos, eta2r_pos, eta2r_pos, eta1r_pos], [ylim_b(1), ylim_b(1), ylim_b(2), ylim_b(2)], ...
      [1, 0.9, 0.9], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta_tilde(1), eta2r_neg, eta2r_neg, eta_tilde(1)], [ylim_b(1), ylim_b(1), ylim_b(2), ylim_b(2)], ...
      [0.9, 0.9, 1], 'EdgeColor','none', 'FaceAlpha',0.22);
patch([eta2r_pos, eta_tilde(end), eta_tilde(end), eta2r_pos], [ylim_b(1), ylim_b(1), ylim_b(2), ylim_b(2)], ...
      [0.9, 0.9, 1], 'EdgeColor','none', 'FaceAlpha',0.22);

% Curves
h_curves_b = gobjects(nXi,1);
for i = 1:nXi
    log_j = log10( max(abs(j_store(i,:)), eps0) );
    h_curves_b(i) = plot(eta_tilde, log_j, '-', 'Color', colors(i,:), ...
        'DisplayName', sprintf('$\\xi = %.1f$', xi_values(i)));
end

% Centerline
xline(0, 'k-', 'LineWidth',1, 'HandleVisibility','off');

% Per-curve asymmetric boundary markers
for i = 1:nXi
    c = 0.75*colors(i,:) + 0.25*[1 1 1];
    draw_vline(ax_b, thr(i,1), ylim_b, ':', c, 1.3);
    draw_vline(ax_b, thr(i,2), ylim_b, ':', c, 1.3);
    draw_vline(ax_b, thr(i,3), ylim_b, '--', c, 1.0);
    draw_vline(ax_b, thr(i,4), ylim_b, '--', c, 1.0);
end

% Saturation asymptote: log10(1)=0
h_sat_b = plot(xlim, [0,0], 'r--', 'DisplayName', 'Saturation asymptote','LineWidth',2.5);

% Linear asymptote in log form: log10(|alpha*lambda*eta|)
eta_lin2 = linspace(eta1r_neg, eta1r_pos, 240);
log_lin  = log10( max(abs(alpha*lambda*eta_lin2), eps0) );
h_lin_b  = plot(eta_lin2, log_lin, 'k--', 'DisplayName', 'Linear asymptote','LineWidth',2.5);

xlabel('$\tilde{\eta}=\lambda f\,\eta_{\rm eff}$');
ylabel('$\log_{10}(|\tilde{j}|)$');
xlim([eta_tilde(1), eta_tilde(end)]);
ylim(ylim_b);

text(ax_b, -0.10, 0.98, 'b', 'Units','normalized', ...
     'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
     'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');

legend(ax_b, [h_curves_b; h_lin_b; h_sat_b], 'Location','southeast', 'Box','off');

% Save
outPng = 'FigureS25.png';
print(fig, '-dpng', '-r600', outPng);
disp(['Saved: ' outPng]);

%% ============================================================
% Local functions
%% ============================================================
function [jtilde, S] = one_mode_schemeA_j(eta, alpha, xi, lambda, maxExp)
    a =  xi      * lambda;
    b = -(1 - xi)* lambda;

    ea = exp( max(min(a.*eta,  maxExp), -maxExp) );
    eb = exp( max(min(b.*eta,  maxExp), -maxExp) );

    Jox  = alpha .* ea;
    Jred = alpha .* eb;

    S = Jox - Jred;
    jtilde = S ./ sqrt(1 + S.^2);
end

function [eta_neg, eta_pos] = threshold_eta_from_S(eta, S, Sthr)
    % Find eta<0 s.t. S(eta)=-Sthr and eta>0 s.t. S(eta)=+Sthr
    % using monotone interpolation on the provided grid.
    eta = eta(:);
    S   = S(:);

    % POSITIVE branch (eta>=0)
    mpos = (eta >= 0);
    etap = eta(mpos);
    Sp   = S(mpos);

    if Sthr <= Sp(1)
        eta_pos = etap(1);
    elseif Sthr >= Sp(end)
        eta_pos = etap(end);
    else
        eta_pos = interp1(Sp, etap, Sthr, 'linear');
    end

    % NEGATIVE branch (eta<=0): target is -Sthr
    mneg = (eta <= 0);
    etan = eta(mneg);
    Sn   = S(mneg);
    target = -Sthr;

    if target <= Sn(1)          % more negative than min(Sn)
        eta_neg = etan(1);
    elseif target >= Sn(end)    % closer to 0 than max(Sn)=0
        eta_neg = etan(end);
    else
        eta_neg = interp1(Sn, etan, target, 'linear');
    end
end

function draw_vline(ax, x0, ylims, ls, col, lw)
    plot(ax, [x0 x0], [ylims(1) ylims(2)], ls, 'Color', col, 'LineWidth', lw, 'HandleVisibility','off');
end
