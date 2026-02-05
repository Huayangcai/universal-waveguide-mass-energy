function FigureS9_nature()
% ===============================================================
% Figure S8 (v2, FULL, RUNNABLE) — aligned with Figure24_DecayDelay_Main_nature_v3/v4
% Key updates vs v1:
%   (1) ThetaSigma unified with Figure24: ThetaSigma = pi/12
%   (2) Archetype markers are NOT hard-coded; they are extracted from the
%       actual F(Omega)=GammaS*GammaL*exp(-2K(Omega)) trajectories at the
%       corridor-selected resonance Omega_res (same logic as Figure24).
%   (3) Keeps manual dual-y overlay (NO yyaxis) for right panel.
%   (4) Quarter-wave label forced to TOP of all labels.
%   (5) "Fabry--Pérot" rendered with Unicode é (Interpreter='none').
% Export: FigureS8_v2.png (600 dpi)
% ===============================================================

clear; clc; close all;

%% ---------- Global style (match Figure24) ----------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultLineLineWidth',2.0);

%% ---------- Parameters ----------
N = 700;
max_tau_tilde = 10;
eps_mag = 1e-12;

% ---- UNIFIED corridor width (matches Figure24) ----
ThetaSigma = pi/12;
alpha_base = 0.22;
alpha_gain = 0.78;

delta_min = 1e-4;
delta_max = 0.999;

kappa_bounds = [0.1, 0.3, 0.7, 0.95, 0.9999];
region_colors = {[0.25 0.80 0.25],[0.85 0.85 0.15],[1.00 0.65 0.20], ...
                 [1.00 0.35 0.10],[0.95 0.20 0.20]};

omega_trt_over4 = pi;    % represents (ω t_rt/4) in SI; constant for demo
Qscale = 1e4;            % plot Q/Qscale
asymShift = 1.05;        % visual shift for asymptotic guide only

contour_levels = [1,2,5];

% ---------- Layout (avoid clipping; give room to right y-axis labels) ----------
layout.leftFrac = 0.44;
layout.gap = 0.085;
layout.margins = [0.09 0.10 0.12 0.07];   % [L R B T]
layout.cbarH = 0.045;
layout.cbarGap = 0.032;

%% =========================================================
% CONFIGS (same as Figure24)
%% =========================================================
cfg = {
    'Quarter-wave',    0.80,              1.00;
    'Fabry--Pérot',    0.80,              0.80;
    'Helmholtz-type',  0.20*exp(1i*pi/2), 0.60;
    'Tunneling-like',  0.85*exp(1i*pi),   0.85;
};
nCfg = size(cfg,1);

% Marker/color scheme for the archetypes (only aesthetics)
mkList  = {'o','s','d','^'};
colList = {
    [0.90 0.15 0.15];   % Quarter-wave
    [0.15 0.75 0.25];   % Fabry--Pérot
    [0.10 0.45 0.95];   % Helmholtz-type
    [0.75 0.20 0.80]    % Tunneling-like
};

% Label offsets (data units on disk)
offList = {
    [ 0.20, -0.22];     % Quarter-wave
    [ 0.28,  0.10];     % Fabry--Pérot
    [ 0.30,  0.14];     % Helmholtz-type
    [ 0.30, -0.10]      % Tunneling-like
};

%% =========================================================
% Compute archetype points from Figure24-like F(Omega) at Omega_res
%% =========================================================
% Use the same dispersion + grid as Figure24 (for strict correspondence)
delta_R = 0.20;
delta_G = 0.10;
Omega_min = 0.05;
Omega_max = 10.0;
N_Omega   = 2200;
Omega = linspace(Omega_min, Omega_max, N_Omega);

[R_Om, Phi_Om] = Ksplit(Omega, delta_R, delta_G);
K = R_Om + 1i*Phi_Om;
E = exp(-2*K);

wrapPi = @(x) mod(x + pi, 2*pi) - pi;

minSep_Om  = 0.40;
maxMarkers = 1;

F_res = zeros(nCfg,1);
Omega_res = zeros(nCfg,1);

for i = 1:nCfg
    GammaS = cfg{i,2};
    GammaL = cfg{i,3};

    F = GammaS.*GammaL.*E;     % Figure24 convention
    D = 1 + F;

    Lam = 1 ./ max(abs(D), eps_mag);

    theta = angle(F);
    det   = abs(wrapPi(theta - pi));

    % --- pick resonance in corridor (use SAME corridor width) ---
    [Ores, ires] = pick_resonances_corridor(Omega, det, Lam, maxMarkers, minSep_Om, 25, ThetaSigma);
    if isempty(ires)
        [~,idx] = max(Lam);
        ires = idx;
        Ores = Omega(idx);
    end

    idxR = ires(1);
    Omega_res(i) = Ores(1);
    F_res(i) = F(idxR);
end

% Reorder to match the intended aesthetic order in colList/mkList:
% We want: Quarter-wave, Fabry--Pérot, Helmholtz-type, Tunneling-like
% cfg already follows that order, so keep as is.

%% =========================================================
% Left panel: F-plane grid (same as v1, but now markers from F_res)
%% =========================================================
[ReF, ImF] = meshgrid(linspace(-1,1,N), linspace(-1,1,N));
kappa = hypot(ReF, ImF);
Theta = atan2(ImF, ReF);
maskDisk = (kappa <= 1);

kappa_safe = kappa;
kappa_safe(~maskDisk) = NaN;
kappa_safe = min(max(kappa_safe, eps_mag), 1-eps_mag);

DeltaTheta = atan2(sin(Theta - pi), cos(Theta - pi));

tauA_tilde = 1 ./ max(-log(kappa_safe), eps_mag);
tauA_tilde(~maskDisk) = NaN;

tau_plot = tauA_tilde;
tau_plot(tau_plot > max_tau_tilde) = max_tau_tilde;

gate = exp(-(DeltaTheta./ThetaSigma).^2);
alphaMap = alpha_base + alpha_gain*gate;
alphaMap(~maskDisk) = 0;
alphaMap = max(0, min(1, alphaMap));

%% =========================================================
% Right panel: δ-domain curves (same as v1)
%% =========================================================
delta = logspace(log10(delta_min), log10(delta_max), 1800);
kappa1 = 1 - delta;
kappa1 = min(max(kappa1, eps_mag), 1-eps_mag);

tauA_1D   = 1 ./ max(-log(kappa1), eps_mag);
tauA_asym = asymShift*(1./delta);

Q_exact = omega_trt_over4 * tauA_1D;
Q_plot  = Q_exact / Qscale;
Q_asym  = asymShift*(omega_trt_over4./delta)/Qscale;

example_k   = [0.1,0.3,0.5,0.7,0.9,0.95,0.99];
example_d   = max(min(1-example_k, delta_max), delta_min);
example_tau = 1 ./ max(-log(example_k), eps_mag);
example_Q   = (omega_trt_over4 * example_tau) / Qscale;

%% =========================================================
% Create figure & axes (manual overlay layout)
%% =========================================================
fig = figure('Color','w','Units','centimeters','Position',[2 0 19.5 12.5]);

ax1       = axes('Parent',fig);  % left
ax2_left  = axes('Parent',fig);  % right: tau
ax2_right = axes('Parent',fig);  % right: Q overlay

cbar1 = [];

apply_layout();
fig.SizeChangedFcn = @(~,~) apply_layout();

%% =========================================================
% PANEL (a): Feedback disk heatmap + archetype markers (from F_res)
%% =========================================================
axes(ax1); cla(ax1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');

hImg = imagesc(ax1, [-1,1], [-1,1], tau_plot);
set(ax1,'YDir','normal');
set(hImg,'AlphaData',alphaMap);

try colormap(ax1, turbo(256)); catch colormap(ax1, parula(256)); end
caxis(ax1,[0 max_tau_tilde]);

th = linspace(0,2*pi,1200);
plot(ax1, cos(th), sin(th), 'k-', 'LineWidth',2.2, 'HandleVisibility','off');

% κ circles
radii = [0.2 0.5 0.8];
for r = radii
    plot(ax1, r*cos(th), r*sin(th), 'w--', 'LineWidth',0.8, 'HandleVisibility','off');
end

% τ̃ contours
tau_cont = tauA_tilde; tau_cont(tau_cont > max_tau_tilde) = NaN;
[~,hC] = contour(ax1, ReF, ImF, tau_cont, contour_levels, ...
    'LineColor',[1 1 1]*0.95,'LineWidth',0.9);
hC.HandleVisibility = 'off';

% resonance direction Θ≈π
plot(ax1,[0 -1.02],[0 0],'w-','LineWidth',1.0,'HandleVisibility','off');

% archetype markers + labels (NOW: use F_res)
ht_all = gobjects(nCfg,1);
hp_all = gobjects(nCfg,1);

for i = 1:nCfg
    pos = [real(F_res(i)), imag(F_res(i))];
    mk  = mkList{i};
    col = colList{i};
    name = cfg{i,1};
    off = offList{i};

    plot(ax1, pos(1), pos(2), mk, 'MarkerSize',10, ...
        'MarkerFaceColor',col,'MarkerEdgeColor','k','LineWidth',1.1);

    labelStr = sprintf('%s', name);
    [ht,hp] = boxed_label(ax1, pos(1)+off(1), pos(2)+off(2), labelStr, ...
        'FontSize',10,'FaceAlpha',0.35,'Interpreter','none');

    ht_all(i)=ht; hp_all(i)=hp;
end

% force Quarter-wave label to top
idxQ = find(strcmp(cfg(:,1),'Quarter-wave'));
if ~isempty(idxQ)
    uistack(hp_all(idxQ),'top');
    uistack(ht_all(idxQ),'top');
end

xlabel(ax1,'$\mathrm{Re}(F)$','FontSize',14);
ylabel(ax1,'$\mathrm{Im}(F)$','FontSize',14);
xlim(ax1,[-1.05 1.05]); ylim(ax1,[-1.05 1.05]); axis(ax1,'square');

% panel label a
text(ax1, -0.18, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

% colorbar (manual position)
cbar1 = colorbar(ax1,'Location','southoutside');
cbar1.Label.Interpreter = 'latex';
cbar1.Label.FontSize = 14;
cbar1.Label.String = '$\tilde{\tau}_A=1/(-\ln\kappa)$';
cbar1.Ticks = [0 5 10];
% ---------- move Helmholtz-type label to right-half real axis ----------
% Find index of Helmholtz-type
idxH = find(strcmp(cfg(:,1),'Helmholtz-type'));
if ~isempty(idxH)
    % Desired label location (right half disk)
    xLab = 0.55;
    yLab = 0.02;

    % Delete old label objects for Helmholtz-type (text + patch)
    if isgraphics(ht_all(idxH)), delete(ht_all(idxH)); end
    if isgraphics(hp_all(idxH)), delete(hp_all(idxH)); end

    % Recreate label at new location
    [htH,hpH] = boxed_label(ax1, xLab, yLab, 'Helmholtz-type', ...
        'FontSize',10,'FaceAlpha',0.35,'Interpreter','none');
    ht_all(idxH) = htH; hp_all(idxH) = hpH;

    % Add a subtle leader line from actual point to label (optional but helpful)
    posH = [real(F_res(idxH)), imag(F_res(idxH))];
    plot(ax1, [posH(1) xLab], [posH(2) yLab], '-', ...
        'Color',[1 1 1]*0.9, 'LineWidth',1.0, 'HandleVisibility','off');

    % Keep Quarter-wave label on top, and also keep Helmholtz label visible
    idxQ = find(strcmp(cfg(:,1),'Quarter-wave'));
    if ~isempty(idxQ)
        if isgraphics(hp_all(idxQ)), uistack(hp_all(idxQ),'top'); end
        if isgraphics(ht_all(idxQ)), uistack(ht_all(idxQ),'top'); end
    end
    if isgraphics(hp_all(idxH)), uistack(hp_all(idxH),'top'); end
    if isgraphics(ht_all(idxH)), uistack(ht_all(idxH),'top'); end
end

apply_layout();

%% =========================================================
% PANEL (b): δ-domain scaling with overlay dual-y axes
%% =========================================================
% bottom: τ̃_A
axes(ax2_left); cla(ax2_left); hold(ax2_left,'on'); box(ax2_left,'on'); grid(ax2_left,'on');
ax2_left.XScale = 'log'; ax2_left.YScale = 'log';
ax2_left.YColor = [0 0.35 0.95];

xlabel(ax2_left,'$\delta=1-\kappa$','FontSize',14);
ylabel(ax2_left,'$\tilde{\tau}_A=1/(-\ln\kappa)$','FontSize',14);

xlim(ax2_left,[delta_min delta_max]);
Tmin = max(min(tauA_1D(tauA_1D>0))*0.8, 1e-3);
Tmax = max(tauA_1D)*1.2;
ylim(ax2_left,[Tmin Tmax]);
yl = ylim(ax2_left);

% shaded κ bands
kappa_edges = [0 kappa_bounds];
for i = 1:numel(kappa_bounds)
    k1 = kappa_edges(i); k2 = kappa_edges(i+1);
    d1 = 1-k2; d2 = 1-k1;
    d1 = max(min(d1,delta_max),delta_min);
    d2 = max(min(d2,delta_max),delta_min);
    if d2 <= d1, continue; end
    patch(ax2_left,[d1 d2 d2 d1],[yl(1) yl(1) yl(2) yl(2)], ...
        region_colors{i},'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');
end

h1 = plot(ax2_left, delta, tauA_1D, '-', 'LineWidth',2.6, 'Color',[0 0.35 0.95]);
plot(ax2_left, delta, tauA_asym, ':', 'LineWidth',1.8, 'Color',[0.15 0.15 0.15], 'HandleVisibility','off');
text(ax2_left, 3e-2, asymShift*(1/3e-2)*1.08, '$\tilde{\tau}_A\sim\delta^{-1}$', ...
    'FontSize',14,'Color',[0.15 0.15 0.15],'BackgroundColor',[1 1 1 0.75]);

plot(ax2_left, example_d, example_tau, 'o', 'MarkerSize',7.5, ...
    'MarkerFaceColor',[0 0.35 0.95], 'MarkerEdgeColor','k','LineWidth',1.0);

% top: Q overlay
axes(ax2_right); cla(ax2_right); hold(ax2_right,'on');
ax2_right.Color = 'none';
ax2_right.Box = 'off';
ax2_right.XScale = 'log'; ax2_right.YScale = 'log';
ax2_right.YAxisLocation = 'right';
ax2_right.YColor = [0.90 0.15 0.15];
ax2_right.XTick = []; ax2_right.XLabel = [];
ax2_right.XGrid = 'off'; ax2_right.YGrid = 'on';

apply_layout();

xlim(ax2_right,[delta_min delta_max]);
Qmin = max(min(Q_plot(Q_plot>0))*0.9, 1e-4);
Qmax = max(Q_plot)*1.2;
ylim(ax2_right,[Qmin Qmax]);

ylabel(ax2_right, sprintf('$Q/10^{%d}$', round(log10(Qscale))), 'FontSize',14);

h2 = plot(ax2_right, delta, Q_plot, '--', 'LineWidth',2.2, 'Color',[0.90 0.15 0.15]);
plot(ax2_right, delta, Q_asym, ':', 'LineWidth',1.8, 'Color',[0.15 0.15 0.15], 'HandleVisibility','off');

text(ax2_right, 0.5e-2, max(Qmin*1.25,1e-4)+0.001, '$Q\sim(\omega t_{\rm rt}/4)\,\delta^{-1}$', ...
    'FontSize',14,'Color',[0.15 0.15 0.15],'BackgroundColor',[1 1 1 0.75]);

plot(ax2_right, example_d, example_Q, 's', 'MarkerSize',7.0, ...
    'MarkerFaceColor',[0.90 0.15 0.15], 'MarkerEdgeColor','k','LineWidth',1.0);

% legend (proxy for Q curve)
h2_proxy = plot(ax2_left, NaN, NaN, '--', 'LineWidth',2.2, 'Color',[0.90 0.15 0.15]);
legend(ax2_left, [h1 h2_proxy], {'Ring-down $\tilde{\tau}_A$', sprintf('$Q/10^{%d}$', round(log10(Qscale)))}, ...
    'Location','northeast','FontSize',14,'Box','off');

% panel label b
text(ax2_left, -0.15, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

apply_layout();

%% ---------- Export ----------
print(fig,'FigureS9','-dpng','-r600');
fprintf('Saved: FigureS9.png\n');

%% ===================== nested layout =====================
function apply_layout()
    L = layout.margins(1); R = layout.margins(2);
    B = layout.margins(3); T = layout.margins(4);
    gap = layout.gap;
    cbH = layout.cbarH;
    cbGap = layout.cbarGap;

    W = 1 - L - R;
    H = 1 - B - T;

    w1 = W*layout.leftFrac;
    w2 = W - w1 - gap;

    axY = B + cbH + cbGap;
    axH = H - (cbH + cbGap);

    set(ax1,'Units','normalized','Position',[L, axY, w1, axH]);

    x2 = L + w1 + gap;
    set(ax2_left ,'Units','normalized','Position',[x2, axY, w2, axH]);
    set(ax2_right,'Units','normalized','Position',[x2, axY, w2, axH]);

    if ~isempty(cbar1) && isgraphics(cbar1)
        cbar1.Units = 'normalized';
        cbar1.Position = [L, B, w1, cbH];
    end
    drawnow limitrate;
end

end

%% ===================== helper: Ksplit (same as Figure24) =====================
function [R, Phi] = Ksplit(Omega, delta_R, delta_G)
P = delta_R*delta_G - Omega.^2;
Q = Omega.*(delta_R + delta_G);
Zabs = sqrt(P.^2 + Q.^2);
R   = sqrt((Zabs + P)/2);
Phi = sqrt((Zabs - P)/2);
R = real(R); Phi = real(Phi);
end

%% ===================== helper: pick resonances with unified corridor =====================
function [Ores, ires] = pick_resonances_corridor(Omega, detAbs, Lam, maxNum, minSep, win, corridor)
Ores = []; ires = [];
cand = find(islocalmax(Lam));
if isempty(cand), return; end

cand = cand(detAbs(cand) < corridor);
if isempty(cand), return; end

[~,ord] = sort(Lam(cand),'descend');
cand = cand(ord);

for k = 1:numel(cand)
    j0 = cand(k);
    j1 = max(1, j0-win);
    j2 = min(numel(Omega), j0+win);
    [~,jj] = max(Lam(j1:j2));
    j = j1 + jj - 1;

    w = Omega(j);
    if isempty(Ores) || all(abs(w - Ores) >= minSep)
        Ores(end+1,1) = w; %#ok<AGROW>
        ires(end+1,1) = j; %#ok<AGROW>
    end
    if numel(Ores) >= maxNum, break; end
end
end

%% ===================== helper: boxed label =====================
function [ht, hp] = boxed_label(ax, x, y, str, varargin)
p = inputParser;
p.addParameter('FontSize',10);
p.addParameter('FaceAlpha',0.35);
p.addParameter('Interpreter','none');
p.parse(varargin{:});
fs = p.Results.FontSize;
fa = p.Results.FaceAlpha;
interp = p.Results.Interpreter;

ht = text(ax, x, y, str, 'Units','data', 'Interpreter',interp, ...
    'FontSize',fs, 'FontWeight','bold', ...
    'Color','w', 'HorizontalAlignment','center', 'VerticalAlignment','middle');

drawnow;
ext = ht.Extent;

dx = 0.08 * ext(3);
dy = 0.25 * ext(4);

x1 = ext(1)-dx; x2 = ext(1)+ext(3)+dx;
y1 = ext(2)-dy; y2 = ext(2)+ext(4)+dy;

hp = patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], [0 0 0], ...
    'EdgeColor','none', 'FaceAlpha',fa);

uistack(hp,'bottom');
uistack(ht,'top');
end
