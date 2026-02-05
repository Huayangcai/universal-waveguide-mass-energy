function Figure03_nature()
% ===============================================================
% Figure (Main) for Section 2.4 (v4):
% - tau_g computed via d/dOmega ln D = D'/D (avoids complex-log branch jumps)
% - Resonance corridor threshold unified: uses ThetaSigma everywhere
% - Removes prctile dependency (custom percentile via sort)
% Export: Figure24_DecayDelay_Main_v4.png (600 dpi)
% ===============================================================

clear; close all; clc;

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultTextFontSize',14);
set(groot,'DefaultLineLineWidth',2.4);

%% ---------------- Resonator settings (demo) ----------------
% [Name, GammaS, GammaL]
cfg = {
    'Quarter-wave',    0.80,              1.00;
    'Fabry--P\''erot', 0.80,              0.80;
    'Helmholtz',       0.20*exp(1i*pi/2), 0.60;
    'Tunnelling',      0.85*exp(1i*pi),   0.85;
};
nCfg = size(cfg,1);

colors = [
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange/red
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560   % purple
];
styles = {'-','--','-.',':'};  % main styles
lws    = [3.2, 2.8, 2.4, 2.2]; % main line widths

%% ---------------- Dispersion parameters ----------------
delta_R = 0.20;
delta_G = 0.10;

%% ---------------- Frequency grid ----------------
Omega_min = 0.05;
Omega_max = 10.0;
N_Omega   = 2200;
Omega = linspace(Omega_min, Omega_max, N_Omega);

[R_Om, Phi_Om] = Ksplit(Omega, delta_R, delta_G);
K = R_Om + 1i*Phi_Om;
E = exp(-2*K);

wrapPi = @(x) mod(x + pi, 2*pi) - pi;

%% ---------------- Resonance corridor controls ----------------
ThetaSigma = pi/12;   % corridor around Theta ~ pi
minSep_Om  = 0.40;
maxMarkers = 3;
eps_mag    = 1e-12;

%% ---------------- Precompute per-config ----------------
F_all     = cell(nCfg,1);
D_all     = cell(nCfg,1);
Lam_all   = cell(nCfg,1);
kap_all   = cell(nCfg,1);
det_all   = cell(nCfg,1);
tauG_all  = cell(nCfg,1);
Ores_all  = cell(nCfg,1);
ires_all  = cell(nCfg,1);

Lam_max = 0;
tauG_abs_max = 0;

for i = 1:nCfg
    GammaS = cfg{i,2};
    GammaL = cfg{i,3};

    F = GammaS.*GammaL.*E;
    D = 1 + F;

    kap   = abs(F);
    theta = angle(F);
    det   = abs(wrapPi(theta - pi));

    Lam = 1./abs(D);

    % -------- tau_g = -Im d/dOmega ln D = -Im(D'/D) -----------
    dD = gradient(D, Omega);

    % Regularize division to avoid |D| ~ 0 blow-up
    Dreg = D;
    m0 = abs(Dreg) < eps_mag;
    if any(m0)
        Dreg(m0) = eps_mag .* exp(1i*angle(Dreg(m0) + 1i*eps_mag));
    end
    dlogD = dD ./ Dreg;
    tauG  = -imag(dlogD);

    % Mild cap to avoid single-point numerical spikes (no toolbox)
    cap = percentile_abs(tauG, 99.7);
    tauG = max(min(tauG, cap), -cap);

    [Ores, ires] = pick_resonances(Omega, det, Lam, maxMarkers, minSep_Om, 25, ThetaSigma);

    F_all{i}    = F;
    D_all{i}    = D;
    Lam_all{i}  = Lam;
    kap_all{i}  = kap;
    det_all{i}  = det;
    tauG_all{i} = tauG;
    Ores_all{i} = Ores;
    ires_all{i} = ires;

    Lam_max = max(Lam_max, max(Lam));
    tauG_abs_max = max(tauG_abs_max, max(abs(tauG)));
end

%% =========================================================
% Panel (a) background on F-disk: tau_tilde = 1/(-ln kappa)
%% =========================================================
N = 700;
max_tau_tilde = 10;

[ReF, ImF] = meshgrid(linspace(-1,1,N), linspace(-1,1,N));
kappa_grid = hypot(ReF, ImF);
Theta_grid = atan2(ImF, ReF);
maskDisk = (kappa_grid <= 1);

kappa_safe = kappa_grid;
kappa_safe(~maskDisk) = NaN;
kappa_safe = min(max(kappa_safe, eps_mag), 1-eps_mag);

DeltaTheta = wrapPi(Theta_grid - pi);

tau_tilde = 1 ./ max(-log(kappa_safe), eps_mag);
tau_plot = tau_tilde;
tau_plot(~maskDisk) = NaN;
tau_plot(tau_plot > max_tau_tilde) = max_tau_tilde;

alpha_base = 0.18;
alpha_gain = 0.82;
gate = exp(-(DeltaTheta./ThetaSigma).^2);
alphaMap = alpha_base + alpha_gain*gate;
alphaMap(~maskDisk) = 0;
alphaMap = max(0, min(1, alphaMap));

%% ---------------- Figure layout ----------------
width_cm  = 24.0;
height_cm = 18.0;
fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

%% ===================== (a) F-disk =====================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); axis(ax1,'equal'); grid(ax1,'on');
hImg = imagesc(ax1, [-1,1], [-1,1], tau_plot);
set(ax1,'YDir','normal');
set(hImg,'AlphaData',alphaMap);

try colormap(ax1, turbo(256)); catch, colormap(ax1, parula(256)); end
caxis(ax1,[0 max_tau_tilde]);

th = linspace(0,2*pi,1200);
plot(ax1, cos(th), sin(th), 'k-', 'LineWidth',2.0, 'HandleVisibility','off');
plot(ax1,[0 -1.02],[0 0],'w-','LineWidth',1.0,'HandleVisibility','off'); % Theta ~ pi

for i = 1:nCfg
    F = F_all{i};
    plot(ax1, real(F), imag(F), styles{i}, 'Color',colors(i,:), ...
        'LineWidth',lws(i), 'DisplayName', cfg{i,1});

    idx = ires_all{i};
    if ~isempty(idx)
        plot(ax1, real(F(idx)), imag(F(idx)), 'o', ...
            'MarkerSize',7.0,'MarkerFaceColor',colors(i,:), ...
            'MarkerEdgeColor','w','LineWidth',0.9,'HandleVisibility','off');
    end
end

xlabel(ax1,'$\mathrm{Re}(F)$'); ylabel(ax1,'$\mathrm{Im}(F)$');
xlim(ax1,[-1.05 1.05]); ylim(ax1,[-1.05 1.05]);

cb1 = colorbar(ax1,'eastoutside');
cb1.Label.Interpreter = 'latex';
cb1.Label.String = '$\tilde{\tau_\mathrm{A}}=1/(-\ln\kappa)$';
cb1.Ticks = [0 5 10];

text(ax1, -0.18, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');
hold(ax1,'off');

%% ===================== (b) Lambda + kappa (yyaxis) =====================
ax2 = nexttile(tlo,2); box(ax2,'on'); grid(ax2,'on'); grid(ax2,'minor');

% ---- left axis: Lambda
yyaxis(ax2,'left'); hold(ax2,'on');
ax2.YColor = [0 0 0];
for i = 1:nCfg
    Lam = Lam_all{i};
    plot(ax2, Omega, Lam, styles{i}, 'Color',colors(i,:), 'LineWidth',lws(i));
    Ores = Ores_all{i};
    for k = 1:numel(Ores)
        xline(ax2, Ores(k), ':', 'LineWidth',0.9, ...
            'Color', lighten_color(colors(i,:),0.60), 'HandleVisibility','off');
    end
end
ylabel(ax2,'$\Lambda(\Omega)=1/|1+F(\Omega)|$');
ylim(ax2,[0 Lam_max*1.08]);

% ---- right axis: kappa
yyaxis(ax2,'right'); hold(ax2,'on');
ax2.YColor = [0.25 0.25 0.25];
for i = 1:nCfg
    kap = kap_all{i};
    plot(ax2, Omega, kap, '-', 'Color', lighten_color(colors(i,:),0.45), ...
        'LineWidth',1.3, 'HandleVisibility','off');
end
ylabel(ax2,'$\kappa(\Omega)=|F(\Omega)|$');
ylim(ax2,[0 1]);
yticks(ax2,0:0.2:1);

% ---- common x
xlabel(ax2,'$\Omega$');
xlim(ax2,[Omega_min Omega_max]);

yyaxis(ax2,'left');
text(ax2, -0.15, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');
hold(ax2,'off');

%% ===================== (c) Ring-down: echoes + envelope =====================
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on'); grid(ax3,'minor');

N_echo = 28;
t_cont = linspace(0,N_echo,1200);

for i = 1:nCfg
    kap  = kap_all{i};
    ires = ires_all{i};
    if isempty(ires)
        [~, idx] = max(Lam_all{i});
    else
        idx = ires(1);
    end

    k0 = max(min(kap(idx),1-eps_mag),eps_mag);
    tauT = 1 / max(-log(k0), eps_mag);

    n = (0:N_echo);
    a_n = (k0).^n;
    env = exp(-t_cont/tauT);

    stem(ax3, n, a_n, 'Marker','none', 'LineWidth',1.0, ...
        'Color', lighten_color(colors(i,:),0.12), 'HandleVisibility','off');

    plot(ax3, t_cont, env, styles{i}, 'Color',colors(i,:), 'LineWidth',lws(i), ...
        'DisplayName', sprintf('%s: $\\tilde{\\tau}_A=%.2f$', cfg{i,1}, tauT));
end

set(ax3,'YScale','log');
ylim(ax3,[1e-4 1.2]);
xlim(ax3,[0 N_echo]);
xlabel(ax3,'Normalised time $\tilde{t}=t/T_{\rm rt}$');
ylabel(ax3,'Ring-down envelope (normalised)');
legend1=legend(ax3,'Location','northeast','Box','off','FontSize',13);
set(legend1,...
    'Position',[0.163867935779504 0.307741495469557 0.258036370315944 0.186846415201823]);
text(ax3, -0.15, 0.98, 'c', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');
hold(ax3,'off');

%% ===================== (d) |tau_g| vs delta, with x-range restriction =====================
ax4 = nexttile(tlo,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on'); grid(ax4,'minor');

delta_corr_all = zeros(0,1);

for i = 1:nCfg
    kap = kap_all{i}(:);
    del = max(1-kap, 1e-6);

    tauG = abs(tauG_all{i}(:));
    det  = det_all{i}(:);

    plot(ax4, del, tauG, styles{i}, 'Color', lighten_color(colors(i,:),0.30), ...
        'LineWidth',1.4, 'HandleVisibility','off');

    idxC = (det < ThetaSigma);
    plot(ax4, del(idxC), tauG(idxC), styles{i}, 'Color', colors(i,:), ...
        'LineWidth', lws(i), 'DisplayName', cfg{i,1});

    if any(idxC)
        delta_corr_all = [delta_corr_all; del(idxC)]; %#ok<AGROW>
    end
end

use_log_delta = false;
if isempty(delta_corr_all)
    all_del = zeros(0,1);
    for i = 1:nCfg
        all_del = [all_del; max(1-kap_all{i}(:), 1e-6)]; %#ok<AGROW>
    end
    dmin = max(min(all_del)*0.90, 1e-3);
    dmax = min(max(all_del)*1.05, 1);
else
    dmin = max(min(delta_corr_all)*0.90, 1e-3);
    dmax = min(max(delta_corr_all)*1.05, 1);
end

if use_log_delta
    set(ax4,'XScale','log');
    xlim(ax4,[max(dmin,1e-3) dmax]);
else
    set(ax4,'XScale','linear');
    xlim(ax4,[dmin dmax]);
end

ylim(ax4,[0 tauG_abs_max*1.08]);
xlabel(ax4,'$\delta(\Omega)=1-\kappa(\Omega)$');
ylabel(ax4,'$|\tau_g(\Omega)|$ (arb.)');
text(ax4, -0.15, 0.98, 'd', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');
hold(ax4,'off');

%% ---------------- Export ----------------
set(findall(gcf,'-property','FontSize'),'FontSize',14);
% print(fig, '-dpng', '-r600', 'Figure24_DecayDelay_Main_v4.png');
% fprintf('Saved: Figure24_DecayDelay_Main_v4.png\n');
    filename = 'Figure03_nature';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);

end

%% ============================================================
function [R, Phi] = Ksplit(Omega, delta_R, delta_G)
P = delta_R*delta_G - Omega.^2;
Q = Omega.*(delta_R + delta_G);
Zabs = sqrt(P.^2 + Q.^2);
R   = sqrt((Zabs + P)/2);
Phi = sqrt((Zabs - P)/2);
R = real(R); Phi = real(Phi);
end

%% ============================================================
function [Ores, ires] = pick_resonances(Omega, detAbs, Lam, maxNum, minSep, win, corridor)
% Picks up to maxNum local maxima of Lam within the phase corridor detAbs<corridor,
% enforcing minimum separation minSep in Omega.
Ores = []; ires = [];

cand = find(islocalmax(Lam));
if isempty(cand), return; end

gate = detAbs(cand) < corridor;
cand = cand(gate);
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

%% ============================================================
function c2 = lighten_color(c, f)
c2 = (1-f)*c + f*[1 1 1];
c2 = max(0, min(1, c2));
end

%% ============================================================
function cap = percentile_abs(x, p)
% Toolbox-free percentile of |x| using sorting (p in [0,100]).
a = abs(x(:));
a = a(isfinite(a));
if isempty(a)
    cap = 0;
    return;
end
a = sort(a);
n = numel(a);
p = max(0, min(100, p));
k = max(1, min(n, round((p/100)*n)));
cap = a(k);
end
