function FigureS2_nature()
%% ============================================================
% Figure S2:Asymmetry-controlled snapshots and extreme envelopes in a lossy waveguide
%% ============================================================
clear; close all; clc;

%% ------------------ Parameters -------------------
delta_R = 1e-1;
delta_G = 3e-1;
Omega   = 10.0;

xi_list = [0.30, 0.50, 0.80];
Lhat    = 1;

Gamma_L = 0.1 * exp(1i*pi/2);
Gamma_S = 0.0 * exp(-1i*pi/6);

theta_list = [0, pi/2, pi, 3*pi/2, 2*pi];

%% ------------------ Dispersion (xi-independent) ---------------
K0  = sqrt((delta_R + 1i*Omega) * (delta_G + 1i*Omega));
K   = Lhat * K0;
E2K = exp(-2*K);

Gammabar_g = (Gamma_L * E2K + Gamma_S) / (1 + Gamma_S * Gamma_L * E2K);
rho   = abs(Gammabar_g);
phi_g = angle(Gammabar_g);

%% ------------------ Spatial grid -----------------------------
Ns = 800;
s  = linspace(0, 1, Ns);
zeta = K .* s;
x    = real(zeta);

%% ------------------ Normalization ----------------------------
Aabs = 1/sqrt(2);
A    = Aabs;
B    = Gammabar_g * A;

Nt   = numel(theta_list);
cols = lines(Nt);
lw_snap = 1.8;
lw_env  = 2.8;

%% ------------------ Figure: 3x2 ------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 18]);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

panelChars = {'a','b','c','d','e','f'};
panelIdx = 0;

% store legend handles/text from FIRST ROW ONLY
legH = [];
legT = {};
axV1 = [];   % the axis on which we create the legend (first-row voltage)

for rr = 1:numel(xi_list)
    xi = xi_list(rr);

    % forward/backward waves
    v_plus  = A .* exp(-2*xi     .* zeta);
    v_minus = B .* exp( 2*(1-xi) .* zeta);

    v_ph = v_plus + v_minus;
    i_ph = v_plus - v_minus;

    % envelopes
    Emax = Aabs * ( exp(-2*xi*x) + rho * exp(2*(1-xi)*x) );
    Emin = Aabs * abs( exp(-2*xi*x) - rho * exp(2*(1-xi)*x) );

    % snapshots
    v_inst = zeros(Nt, Ns);
    i_inst = zeros(Nt, Ns);
    for k = 1:Nt
        th = theta_list(k);
        v_inst(k,:) = real( v_ph .* exp(1i*th) );
        i_inst(k,:) = real( i_ph .* exp(1i*th) );
    end

    %% ========== Col 1: Voltage ==========
    axV = nexttile((rr-1)*2 + 1); hold(axV,'on'); box(axV,'on'); grid(axV,'on');
    axV.GridAlpha=0.2; axV.MinorGridAlpha=0.1; axV.XMinorGrid='on'; axV.YMinorGrid='on';

    hSnapV = gobjects(Nt,1);
    for k = 1:Nt
        hSnapV(k) = plot(axV, s, v_inst(k,:), 'LineWidth', lw_snap, 'Color', cols(k,:));
    end

    fill(axV, [s fliplr(s)], [Emax fliplr(Emin)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
    fill(axV, [s fliplr(s)], [-Emin fliplr(-Emax)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');

    hEmax = plot(axV, s, +Emax, 'k--', 'LineWidth', lw_env);
    plot(axV, s, -Emax, 'k--', 'LineWidth', lw_env, 'HandleVisibility','off');
    hEmin = plot(axV, s, +Emin, 'k:',  'LineWidth', 2.0);
    plot(axV, s, -Emin, 'k:',  'LineWidth', 2.0, 'HandleVisibility','off');

    ylabel(axV,'$v(s,t)$','Interpreter','latex');
    xlim(axV,[0 1]);

    panelIdx = panelIdx + 1;
    text(axV, -0.12, 0.94, panelChars{panelIdx}, 'Units','normalized', ...
        'FontWeight','bold', 'Interpreter','none');
    text(axV, 0.98, 0.94, sprintf('$\\xi=%.2f$', xi), 'Units','normalized', ...
        'HorizontalAlignment','right','Interpreter','latex');

    if rr < numel(xi_list)
        axV.XTickLabel = [];
    else
        xlabel(axV,'$s=x/L$','Interpreter','latex');
    end

    % ---- collect legend ONLY from first row ----
    if rr == 1
        axV1 = axV;
        legPhase = {'$\theta=0$','$\theta=\pi/2$','$\theta=\pi$','$\theta=3\pi/2$','$\theta=2\pi$'};
        legH = [hSnapV(:); hEmax; hEmin];
        legT = [legPhase(:); {'$E_{\max}$'; '$E_{\min}$'}];
    else
        set(hSnapV,'HandleVisibility','off');
        set(hEmax,'HandleVisibility','off');
        set(hEmin,'HandleVisibility','off');
    end

    %% ========== Col 2: Current ==========
    axI = nexttile((rr-1)*2 + 2); hold(axI,'on'); box(axI,'on'); grid(axI,'on');
    axI.GridAlpha=0.2; axI.MinorGridAlpha=0.1; axI.XMinorGrid='on'; axI.YMinorGrid='on';

    for k = 1:Nt
        h = plot(axI, s, i_inst(k,:), 'LineWidth', lw_snap, 'Color', cols(k,:));
        set(h,'HandleVisibility','off'); % never in legend
    end

    fill(axI, [s fliplr(s)], [Emax fliplr(Emin)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
    fill(axI, [s fliplr(s)], [-Emin fliplr(-Emax)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');

    plot(axI, s, +Emax, 'k--', 'LineWidth', lw_env, 'HandleVisibility','off');
    plot(axI, s, -Emax, 'k--', 'LineWidth', lw_env, 'HandleVisibility','off');
    plot(axI, s, +Emin, 'k:',  'LineWidth', 2.0,   'HandleVisibility','off');
    plot(axI, s, -Emin, 'k:',  'LineWidth', 2.0,   'HandleVisibility','off');

    ylabel(axI,'$i(s,t)$','Interpreter','latex');
    xlim(axI,[0 1]);

    panelIdx = panelIdx + 1;
    text(axI, -0.12, 0.94, panelChars{panelIdx}, 'Units','normalized', ...
        'FontWeight','bold', 'Interpreter','none');

    if rr < numel(xi_list)
        axI.XTickLabel = [];
    else
        xlabel(axI,'$s=x/L$','Interpreter','latex');
    end
end
%% -------- One global legend (top-center, single row) ----------
% Create legend on first-row voltage axis, then manually move it to top center.
lg = legend(axV1, legH, legT, 'Interpreter','latex', 'Orientation','horizontal');
lg.Box = 'off';
lg.FontSize = 11;
lg.NumColumns = numel(legT);     % force single row (if too long, MATLAB may still wrap)

% Move legend to top-center of the FIGURE (normalized figure coordinates)
lg.Units = 'normalized';
lg.Position = [0.10, 0.965, 0.80, 0.03];  % [x y w h], adjust if needed

% Give some extra top margin so legend doesn't overlap tiles (older MATLAB)
fig.Position(4) = fig.Position(4) + 1.5;  % add a bit of height (cm)

%% ------------------ Final font + save ------------------------
set(findall(gcf,'-property','FontSize'),'FontSize',14);
filename = 'Figure S2';
print(fig,'-dpng','-r600',[filename,'.png']);
fprintf('Saved: %s.png\n', filename);

%% ------------------ Summary -------------------
fprintf('\n=== Summary (common across rows) ===\n');
fprintf('K0 = %.4f + i%.4f\n', real(K0), imag(K0));
fprintf('K  = %.4f + i%.4f  (Lhat=%g)\n', real(K), imag(K), Lhat);
fprintf('|Gammabar_g| = %.4f, angle = %.2f deg\n', rho, rad2deg(phi_g));
end

