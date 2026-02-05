function FigureS6_nature()
% ===============================================================
% Figure S5 (UPDATED): Resonance-aware useful load power landscape
% and its maximum at fixed attenuation (1x2, Nature-style)
%
% Left:  T_L(R,|Gamma_L|) on the odd-pi (Theta=pi) resonance class
%        + gray dashed absorptivity contours A=0.1,0.5,0.9
%        + overlay |Gamma_L|_opt and (feasible) |Gamma_L|_cc
% Right: analytic maximum T_L,max(R) and comparison to "critical-coupling"
%        branch (where feasible)
%
% Core definitions (S7-consistent):
%   E = exp(-2K) = exp(-2R) exp(-i2Phi)
%   F = Gamma_S Gamma_L E
%   Gamma_g = (Gamma_L E + Gamma_S)/(1 + F)
%   A = 1 - |Gamma_g|^2
%   T_L = A * (1-|Gamma_L|^2) exp(-2R) / |1+F|^2
%   eta_use = T_L/A = (1-|Gamma_L|^2) exp(-2R) / |1+F|^2
%   T_int = A - T_L
%
% Analytic (Theta=pi):
%   |Gamma_L|_opt = |Gamma_S| exp(-2R)
%   T_L,max(R) = (1-|Gamma_S|^2) exp(-2R) / (1 - |Gamma_S|^2 exp(-4R))
%   |Gamma_L|_cc  = |Gamma_S| exp(+2R) (may be infeasible if >1)
% ===============================================================

clear; clc; close all;

%% -------------------- Global style --------------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultTextFontSize',14);
set(groot,'DefaultLineLineWidth',2.2);

%% -------------------- User parameters --------------------
gS    = 0.70;      % |Gamma_S|
phi_S = 0.0;       % arg(Gamma_S) (take 0 for clean resonance slice)
Phi0  = 0.0;       % fix Phi for the slice (Theta enforced via phi_L)
Rmin  = 0.0;
Rmax  = 1.2;       % show beyond feasibility of critical coupling for gS=0.7

Gamma_S = gS * exp(1i*phi_S);

% Enforce odd-pi class: Theta = arg(Gamma_S Gamma_L) - 2Phi0 = pi
phi_L_res = wrapToPi_local(pi + 2*Phi0 - phi_S);

% Critical-coupling feasibility boundary: |Gamma_L|_cc = gS e^{2R} <= 1
Rcc = 0.5*log(1/max(gS,1e-12));

%% -------------------- Figure layout --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 23.5 14]);
tlo = tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');

%% ===============================================================
% (a) Useful load power map T_L over (R, |Gamma_L|) on Theta=pi slice
% ===============================================================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on');

NR  = 520;
NGL = 420;
R   = linspace(Rmin, Rmax, NR);

gL_max = 0.9995;
gL     = linspace(0, gL_max, NGL);
[RR, GL] = meshgrid(R, gL);

% E = exp(-2K) with Phi fixed
E = exp(-2*RR) .* exp(-1i*2*Phi0);

% Complex load reflection on the resonance class
Gamma_L = GL .* exp(1i*phi_L_res);

% Feedback factor and composed reflection
F = Gamma_S .* Gamma_L .* E;
Gamma_g = (Gamma_L .* E + Gamma_S) ./ (1 + F);

% Absorptivity
A = 1 - abs(Gamma_g).^2;
A = max(A, 0);

% Useful load power fraction
eta_use = (1 - abs(Gamma_L).^2) .* exp(-2*RR) ./ abs(1 + F).^2;
eta_use = max(min(real(eta_use),1),0);
T_L = A .* eta_use;

% Mask any numerical non-physical values
maskPhys = isfinite(T_L) & (T_L >= 0);
Tplot    = nan(size(T_L));
Tplot(maskPhys) = T_L(maskPhys);

% Choose caxis matched to reachable set, rounded to 2 decimals
TLmin = min(Tplot(:),[],'omitnan');
TLmax = max(Tplot(:),[],'omitnan');
clim  = [floor(100*TLmin)/100, ceil(100*TLmax)/100];
if clim(1) == clim(2), clim = [clim(1) clim(1)+0.01]; end

fprintf('TL color scale matched to reachable set: [%.5f, %.5f]\n', TLmin, TLmax);

% Heatmap
hImg = imagesc(ax1, R, gL, Tplot);
set(ax1,'YDir','normal');
alphaImg = zeros(size(Tplot)); alphaImg(maskPhys) = 1;
set(hImg,'AlphaData',alphaImg);

try colormap(ax1, turbo(256)); catch colormap(ax1, parula(256)); end
caxis(ax1, clim);

cb = colorbar(ax1,'eastoutside');
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'none';
cb.Label.String = '$\mathcal{T}_L(R,|\Gamma_L|)$';
format_colorbar_2dp(cb);

xlabel(ax1,'$R=\Re[K(\Omega)]$');
ylabel(ax1,'$|\Gamma_L(\Omega)|$');

% --- Gray dashed contours of absorptivity A at 0.1,0.5,0.9
% levelsA = [0.1 0.5 0.9];
% [C,hc] = contour(ax1, R, gL, A, levelsA, 'LineColor',[0.55 0.55 0.55], ...
%     'LineStyle','--','LineWidth',1.4);
% hc.HandleVisibility = 'off';
% 
% % Place labels near lower-left automatically (no manual clicking)
% place_contour_labels(ax1, R, gL, A, levelsA);

% --- Overlay: |Gamma_L|_opt and |Gamma_L|_cc (feasible segment), and |Gamma_L|=1
gL_opt = gS * exp(-2*R);      % always feasible
gL_cc  = gS * exp( 2*R);      % may exceed 1
idxFeas = (gL_cc <= 1);

h_opt = plot(ax1, R, gL_opt, 'k-', 'LineWidth', 2.8);
h_cc  = plot(ax1, R(idxFeas), gL_cc(idxFeas), 'w--', 'LineWidth', 2.8);
% h_sat = plot(ax1, R, ones(size(R)), 'k:', 'LineWidth', 2.0);

% Mark feasibility boundary for critical coupling
if Rcc >= Rmin && Rcc <= Rmax
    xline(ax1, Rcc, 'w:', 'LineWidth', 2.0, 'HandleVisibility','off');
    text(ax1, Rcc, 0.97, '$|\Gamma_L|_{\rm cc}=1$', ...
        'Color','w','Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'BackgroundColor',[0 0 0 0.25],'EdgeColor','none');
end

xlim(ax1,[Rmin Rmax]); ylim(ax1,[0 1]);
grid(ax1,'on'); grid(ax1,'minor');

legend(ax1, [h_opt h_cc], ...
    { ...
    '$|\Gamma_L|_{\rm opt}=|\Gamma_S|e^{-2R}$', ...
    '$|\Gamma_L|_{\rm cc}=|\Gamma_S|e^{2R}$ (feasible)', ...
    }, ...
    'Location','east','Box','off','Interpreter','latex');

% Panel label
text(ax1, -0.18, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

%% ===============================================================
% (b) Maximum useful load power vs R (analytic) + comparison
% ===============================================================
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on');

Rline = linspace(max(Rmin,1e-6), Rmax, 900);

% Analytic maximum on Theta=pi slice
TL_max = (1 - gS^2) .* exp(-2*Rline) ./ (1 - gS^2 .* exp(-4*Rline));

% Value at critical coupling branch (where feasible): |Gamma_L|_cc=gS e^{2R}
gL_cc_line = gS .* exp(2*Rline);
TL_cc = nan(size(Rline));
idx = (gL_cc_line <= 1);
gLtmp = gL_cc_line(idx);
kappa = gS .* gLtmp .* exp(-2*Rline(idx));  % kappa = |Gamma_S||Gamma_L|e^{-2R}
TL_cc(idx) = (1-gS^2).*(1-gLtmp.^2).*exp(-2*Rline(idx)) ./ (1 - kappa).^2;

h1 = plot(ax2, Rline, TL_max, 'k-', 'LineWidth', 3.0);
h1.DisplayName = '$\mathcal{T}_{L,\max}(R)$ (analytic, $\Theta=\pi$)';

h2 = plot(ax2, Rline, TL_cc, '--', 'Color',[0.80 0.20 0.20], 'LineWidth', 2.4);
h2.DisplayName = '$\mathcal{T}_L$ at critical coupling (if feasible)';

% Mark Rcc
if Rcc >= Rmin && Rcc <= Rmax
    xline(ax2, Rcc, 'k:', 'LineWidth', 1.8, 'HandleVisibility','off');
    text(ax2, Rcc+0.1, 0.05, '$R_{\rm cc}=\frac{1}{2}\ln(1/|\Gamma_S|)$', ...
        'Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'BackgroundColor',[1 1 1 0.70],'EdgeColor','none');
end

% Mark global maximum over the plotted range
[TLpeak, ip] = max(TL_max);
Rpeak = Rline(ip);
% plot(ax2, Rpeak, TLpeak, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7);
% text(ax2, Rpeak, TLpeak, sprintf('  max = %.2f', TLpeak), ...
%     'Interpreter','none','FontSize',11,'VerticalAlignment','middle');

grid(ax2,'on'); grid(ax2,'minor');
xlim(ax2,[Rmin Rmax]);
ylim(ax2,[0 1.05]);

xlabel(ax2,'$R=\Re[K(\Omega)]$');
ylabel(ax2,'$\mathcal{T}_L$');

legend(ax2,'Location','northeast','Box','off','Interpreter','latex','FontSize',11);

% Panel label
text(ax2, -0.18, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

%% -------------------- Export --------------------
filename = 'Figure S6';
print(fig,'-dpng','-r600',[filename,'.png']);
fprintf('Saved: %s.png\n', filename);

end

%% ===================== helpers =====================
function x = wrapToPi_local(x)
x = mod(x + pi, 2*pi) - pi;
end

function format_colorbar_2dp(cb)
drawnow;
ticks = cb.Ticks;
if isempty(ticks), return; end
cb.TickLabels = compose('%.2f', ticks);
cb.TickLabelInterpreter = 'none';
end

function place_contour_labels(ax, R, gL, A, levelsA)
% Find near-lower-left points where A ~ level and place labels there
[RR, GL] = meshgrid(R, gL);

% focus region (left-lower): you can tighten/loosen these
mask = (RR <= (min(R)+0.45*(max(R)-min(R)))) & (GL <= 0.55) & isfinite(A);

for k = 1:numel(levelsA)
    lvl = levelsA(k);
    D = abs(A - lvl);
    D(~mask) = inf;
    [~, idx] = min(D(:));
    if isinf(D(idx)), continue; end
    rp = RR(idx); gp = GL(idx);

    % small marker and label
    plot(ax, rp, gp, 'o', 'MarkerSize',3.5, 'MarkerFaceColor',[0.55 0.55 0.55], ...
        'MarkerEdgeColor',[0.45 0.45 0.45], 'HandleVisibility','off');
    text(ax, rp, gp, sprintf('  %.1f', lvl), ...
        'Color',[0.45 0.45 0.45], 'FontSize',11, 'Interpreter','none', ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'BackgroundColor',[1 1 1 0.55], 'EdgeColor','none');
end
end
