function FigureS11_nature()
% ============================================================
% 8-panel dashboard (a-h), with consistent color/linestyle mapping:
%   - NG (nongeneric): BLUE
%   - G  (generic):    ORANGE
%   - fixed: SOLID
%   - SVD:   DASHED
%
% Panels:
%   a: channel "mass" rho = sqrt(P_out): fixed vs SVD (NG & G)
%   b: Fig.3B-like one-channel output P_out: fixed vs SVD (NG & G)
%   c: log-log near-zero scaling of TWO-CHANNEL total output: fixed vs SVD (NG & G)
%   d: TWO-CHANNEL total output vs detuning: fixed vs SVD (NG & G)
%   e: Nongeneric two-port outputs + total: fixed vs SVD
%   f: Generic   two-port outputs + total: fixed vs SVD
%   g: Cai–Smith disk (nongeneric): Gamma_fixed vs Gamma_svd
%   h: Cai–Smith disk (generic):   Gamma_fixed vs Gamma_svd
%
% Definitions:
%   - fixed coherent setting: a0 = v_min from SVD of S(0), held constant.
%   - SVD-probe: a*(delta)=v_min(delta) at each delta.
%   - one-channel P_out: ||S a||^2  (fixed)  vs sigma_min^2 (SVD).
%   - Cai–Smith proxy complex amplitude: Gamma(delta)=a^H S a (Rayleigh quotient).
% ============================================================

clc; close all;

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultLineLineWidth',1.8);

omega0 = 0;
Nd = 5001;

% style (consistent mapping)
cNG = [0 0.4470 0.7410];   % MATLAB default blue
cG  = [0.8500 0.3250 0.0980]; % MATLAB default orange
lw  = 2.0;

%% -------------------- ONE-CHANNEL (Fig.3B-like) parameters --------------------
g1 = 400; g2 = 1200;

% nongeneric (symmetric)
gcN = (g1+g2)/2;
kappaN = abs(g1-g2)/4;
cfgBN = make_cfg(omega0,g1,g2,gcN,gcN,kappaN);

% generic (asymmetric)
gc1G = 1200;
gc2G = (g1+g2) - gc1G;
kappaG = abs(g1-gc1G)/2;
cfgBG = make_cfg(omega0,g1,g2,gc1G,gc2G,kappaG);

dMaxB = 650;
deltaB = linspace(-dMaxB,dMaxB,Nd).';

% one-channel outputs + Cai–Smith proxy Gamma for fixed/SVD
ONE_N = eval_onechannel_fixed_svd(deltaB,cfgBN);
ONE_G = eval_onechannel_fixed_svd(deltaB,cfgBG);

% exponents (display on subplot b)
expo1N_fixed = fit_exponent(deltaB, ONE_N.P_fixed, [1,80]);
expo1N_svd   = fit_exponent(deltaB, ONE_N.P_svd,   [1,80]);
expo1G_fixed = fit_exponent(deltaB, ONE_G.P_fixed, [1,80]);
expo1G_svd   = fit_exponent(deltaB, ONE_G.P_svd,   [1,80]);

%% -------------------- TWO-CHANNEL (Fig.3D-like) parameters --------------------
g1D = 800; g2D = 1600;

% nongeneric
gcD_N = (g1D+g2D)/2;
kappaD_N = abs(g1D-g2D)/4;
cfgDN = make_cfg(omega0,g1D,g2D,gcD_N,gcD_N,kappaD_N);

% generic
gc1D_G = 1800;
gc2D_G = (g1D+g2D) - gc1D_G;
kappaD_G = abs(g1D-gc1D_G)/2;
cfgDG = make_cfg(omega0,g1D,g2D,gc1D_G,gc2D_G,kappaD_G);

dMaxD = 4000;
deltaD = linspace(-dMaxD,dMaxD,Nd).';

mismatch = struct('amp',0.00,'phase_deg',0.0);

% FIXED: a0=vmin(S(0))
RES_DN_fixed = eval_outputs_fixed(deltaD,cfgDN,mismatch,'svd_at0');
RES_DG_fixed = eval_outputs_fixed(deltaD,cfgDG,mismatch,'svd_at0');

% SVD-probe: a*(delta)=vmin(S(delta))
RES_DN_svd   = eval_case_outputs(deltaD,cfgDN,mismatch,'svd');
RES_DG_svd   = eval_case_outputs(deltaD,cfgDG,mismatch,'svd');

% exponents (two-channel totals)
expoDN_fixed = fit_exponent(deltaD,RES_DN_fixed.Pout_tot,[5,150]);
expoDG_fixed = fit_exponent(deltaD,RES_DG_fixed.Pout_tot,[5,150]);
expoDN_svd   = fit_exponent(deltaD,RES_DN_svd.Pout_tot,  [5,150]);
expoDG_svd   = fit_exponent(deltaD,RES_DG_svd.Pout_tot,  [5,150]);

%% ============================================================
% 8-panel figure (4x2)
% ============================================================
fig = figure('Color','w','Units','centimeters','Position',[2 2 23 21]);
tlo = tiledlayout(fig, 4, 2, 'Padding','compact', 'TileSpacing','compact');
nTiles = 8;
letters = char('a' + (0:(nTiles-1)));  % a..l

% (a) channel mass rho = sqrt(P_out): fixed vs SVD (NG & G)
ax1 = nexttile(tlo, 1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(deltaB, sqrt(max(ONE_N.P_fixed,0)), '--','Color',cNG, 'LineWidth',lw); hold on;
plot(deltaB, sqrt(max(ONE_N.P_svd,  0)), '-','Color',cNG, 'LineWidth',lw);
plot(deltaB, sqrt(max(ONE_G.P_fixed,0)), '--','Color',cG,  'LineWidth',lw);
plot(deltaB, sqrt(max(ONE_G.P_svd,  0)), '-','Color',cG,  'LineWidth',lw);
grid on;
xlabel('$\Delta \omega$'); ylabel('$\rho(\Delta \omega)$');
% title('Waveguide channel mass $\rho$: fixed (solid) vs SVD (dashed)');
legend('NG fixed','NG SVD','G fixed','G SVD','Location','north');
panel_letter_local(ax1, letters(1));

% (b) Fig.3B-like one-channel output P_out: fixed vs SVD (NG & G)
ax2 = nexttile(tlo, 2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(deltaB, ONE_N.P_fixed, '--','Color',cNG,'LineWidth',lw); hold on;
plot(deltaB, ONE_N.P_svd,   '-','Color',cNG,'LineWidth',lw);
plot(deltaB, ONE_G.P_fixed, '--','Color',cG, 'LineWidth',lw);
plot(deltaB, ONE_G.P_svd,   '-','Color',cG, 'LineWidth',lw);
grid on;
xlabel('$\Delta \omega$'); ylabel('$P_{\mathrm{out}}/P_{\mathrm{in}}$');
% title('Fig.3B-like one-channel output: fixed vs SVD');
legend( ...
    sprintf('NG fixed (%.2f)',expo1N_fixed), ...
    sprintf('NG SVD (%.2f)',expo1N_svd), ...
    sprintf('G fixed (%.2f)',expo1G_fixed), ...
    sprintf('G SVD (%.2f)',expo1G_svd), ...
    'Location','north');
panel_letter_local(ax2, letters(2));
% (c) Log-log near-zero scaling: TWO-CHANNEL total outputs (fixed vs SVD; NG & G)
ax3 = nexttile(tlo, 3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
mask = abs(deltaD)>=5 & abs(deltaD)<=200;

loglog(abs(deltaD(mask)), max(RES_DN_fixed.Pout_tot(mask),1e-30), '--','Color',cNG,'LineWidth',1.4); hold on;
loglog(abs(deltaD(mask)), max(RES_DN_svd.Pout_tot(mask),  1e-30), '-','Color',cNG,'LineWidth',1.4);
loglog(abs(deltaD(mask)), max(RES_DG_fixed.Pout_tot(mask),1e-30), '--','Color',cG, 'LineWidth',1.4);
loglog(abs(deltaD(mask)), max(RES_DG_svd.Pout_tot(mask),  1e-30), '-','Color',cG, 'LineWidth',1.4);

% slope guides
xline = logspace(log10(5),log10(200),160);
x0 = 40;
y0_2 = median(max(RES_DN_fixed.Pout_tot(abs(deltaD)>30 & abs(deltaD)<50),1e-30));
y0_4 = median(max(RES_DG_svd.Pout_tot(  abs(deltaD)>30 & abs(deltaD)<50),1e-30));
loglog(xline, y0_2*(xline/x0).^2, ':', 'Color',[0.2 0.2 0.2], 'LineWidth',1.2);
loglog(xline, y0_4*(xline/x0).^4, '-.', 'Color',[0.2 0.2 0.2], 'LineWidth',2.0);

grid on;
xlabel('$\Delta \omega$'); ylabel('$P_{\mathrm{out}}/P_{\mathrm{in}}$');
% title('Near-zero scaling (two-channel totals): fixed vs SVD');
legend( ...
    sprintf('NG fixed (%.2f)',expoDN_fixed), ...
    sprintf('NG SVD (%.2f)',expoDN_svd), ...
    sprintf('G fixed (%.2f)',expoDG_fixed), ...
    sprintf('G SVD (%.2f)',expoDG_svd), ...
    'slope=2','slope=4','Location','northwest');
panel_letter_local(ax3, letters(3));
% (d) TWO-CHANNEL total output vs detuning: fixed vs SVD (NG & G)
ax4 = nexttile(tlo, 4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
plot(deltaD, RES_DN_fixed.Pout_tot, '--','Color',cNG,'LineWidth',lw); hold on;
plot(deltaD, RES_DN_svd.Pout_tot,   '-','Color',cNG,'LineWidth',lw);
plot(deltaD, RES_DG_fixed.Pout_tot, '--','Color',cG, 'LineWidth',lw);
plot(deltaD, RES_DG_svd.Pout_tot,   '-','Color',cG, 'LineWidth',lw);
grid on;
xlabel('$\Delta \omega$'); ylabel('$P_{\mathrm{out}}/P_{\mathrm{in}}$');
% title('Fig.3D-like total output: fixed vs SVD (NG \& G)');
legend('NG fixed','NG SVD','G fixed','G SVD','Location','north');
panel_letter_local(ax4, letters(4));
% (e) Nongeneric: port outputs + total, fixed vs SVD (line style rule)
ax5 = nexttile(tlo, 5); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');
% port colors fixed across modes for readability
cP1 = [0.2 0.6 0.9];
cP2 = [0.9 0.4 0.2];
plot(deltaD, RES_DN_fixed.Pout_comp1, '--','Color',cP1,'LineWidth',1.8); hold on;
plot(deltaD, RES_DN_svd.Pout_comp1,   '-','Color',cP1,'LineWidth',1.6);
plot(deltaD, RES_DN_fixed.Pout_comp2,'-', 'Color',cP2,'LineWidth',1.8);
plot(deltaD, RES_DN_svd.Pout_comp2,   '--','Color',cP2,'LineWidth',1.6);
plot(deltaD, RES_DN_fixed.Pout_tot,   'k--','LineWidth',2.2);
plot(deltaD, RES_DN_svd.Pout_tot,     'k-','LineWidth',2.0);
grid on;
xlabel('$\Delta \omega$'); ylabel('$P_{\mathrm{out}}/P_{\mathrm{in}}$');
% title('Nongeneric: two-port outputs + total (solid=fixed, dashed=SVD)');
legend1=legend('$|\mathbf{b_1}|^2$ fixed','$|\mathbf{b_1}|^2$ SVD','$|\mathbf{b_2}|^2$ fixed','$|\mathbf{b_2}|^2$ SVD','total fixed','total SVD', ...
       'Location','north');
set(legend1,...
    'Position',[0.205682517727068 0.356411532129316 0.116748148678271 0.132381753096753]);
panel_letter_local(ax5, letters(5));
% (f) Generic: port outputs + total, fixed vs SVD (same line style rule)
ax6 = nexttile(tlo, 6); hold(ax6,'on'); box(ax6,'on'); grid(ax6,'on');
plot(deltaD, RES_DG_fixed.Pout_comp1, '--','Color',cP1,'LineWidth',1.8); hold on;
plot(deltaD, RES_DG_svd.Pout_comp1,   '-','Color',cP1,'LineWidth',1.6);
plot(deltaD, RES_DG_fixed.Pout_comp2, '-','Color',cP2,'LineWidth',1.8);
plot(deltaD, RES_DG_svd.Pout_comp2,   '--','Color',cP2,'LineWidth',1.6);
plot(deltaD, RES_DG_fixed.Pout_tot,   'k--','LineWidth',2.2);
plot(deltaD, RES_DG_svd.Pout_tot,     'k-','LineWidth',2.0);
grid on;
xlabel('$\Delta \omega$'); ylabel('$P_{\mathrm{out}}/P_{\mathrm{in}}$');
% title('Generic: two-port outputs + total (solid=fixed, dashed=SVD)');
legend1=legend('$|\mathbf{b_1}|^2$ fixed','$|\mathbf{b_1}|^2$ SVD','$|\mathbf{b_2}|^2$ fixed','$|\mathbf{b_2}|^2$ SVD','total fixed','total SVD', ...
       'Location','north');
set(legend1,...
    'Position',[0.691195768099409 0.35976831054894 0.116748148678271 0.132381753096753]);
panel_letter_local(ax6, letters(6));
% (g) Cai–Smith disk: nongeneric (fixed vs SVD trajectories)
ax7 = nexttile(tlo, 7); hold(ax7,'on'); box(ax7,'on'); grid(ax7,'on');
plot(real(ONE_N.Gamma_fixed), imag(ONE_N.Gamma_fixed), ':','Color',cNG,'LineWidth',1.8); hold on;
plot(real(ONE_N.Gamma_svd),   imag(ONE_N.Gamma_svd),   '-','Color',cNG,'LineWidth',1.8);
viscircles([0 0], 1, 'LineStyle','-', 'Color',[0.85 0 0], 'LineWidth',1.6);
axis equal; grid on;
xlabel('$\Re(\Gamma)$'); ylabel('$\Im(\Gamma)$');
% title('Cai--Smith disk (NG): $\Gamma=a^\dagger S a$ (solid=fixed, dashed=SVD)');
legend('NG fixed','NG SVD','unit circle','Location','best');
panel_letter_local(ax7, letters(7));
% (h) Cai–Smith disk: generic (fixed vs SVD trajectories)
ax8 = nexttile(tlo, 8); hold(ax8,'on'); box(ax8,'on'); grid(ax8,'on');
plot(real(ONE_G.Gamma_fixed), imag(ONE_G.Gamma_fixed), ':','Color',cG,'LineWidth',1.8); hold on;
plot(real(ONE_G.Gamma_svd),   imag(ONE_G.Gamma_svd),   '-','Color',cG,'LineWidth',1.8);
viscircles([0 0], 1, 'LineStyle','-', 'Color',[0.85 0 0], 'LineWidth',1.6);
axis equal; grid on;
xlabel('$\Re(\Gamma)$'); ylabel('$\Im(\Gamma)$');
% title('Cai--Smith disk (G): $\Gamma=a^\dagger S a$ (solid=fixed, dashed=SVD)');
legend('G fixed','G SVD','unit circle','Location','best');
panel_letter_local(ax8, letters(8));

print(fig, '-dpng', '-r600', 'Figure S11');
disp(['Saved: Figure S11']);
end

%% ===================== Helper functions =========================
function cfg = make_cfg(omega0,g1,g2,gc1,gc2,kappa)
cfg = struct();
cfg.omega0 = omega0;
cfg.g1 = g1; cfg.g2 = g2;
cfg.gc1 = gc1; cfg.gc2 = gc2;
cfg.kappa = kappa;
end

function S = S_TCMT_2x2(omega,cfg)
D1 = (omega-cfg.omega0) + 1i*(cfg.g1+cfg.gc1)/2;
D2 = (omega-cfg.omega0) + 1i*(cfg.g2+cfg.gc2)/2;
den = (D1.*D2 - cfg.kappa^2);
M = [ cfg.gc1.*D2, sqrt(cfg.gc1*cfg.gc2).*cfg.kappa; ...
      sqrt(cfg.gc1*cfg.gc2).*cfg.kappa, cfg.gc2.*D1 ];
S = eye(2) - 1i*(M./den);
end

function OUT = eval_onechannel_fixed_svd(delta,cfg)
% Returns:
%   P_fixed(delta)=||S(delta)a0||^2, where a0=vmin(S(0))
%   P_svd(delta)=sigma_min(delta)^2
%   Gamma_fixed(delta)=a0^H S(delta) a0
%   Gamma_svd(delta)=a*(delta)^H S(delta) a*(delta), a*(delta)=vmin(S(delta))

Nd = numel(delta);

S0 = S_TCMT_2x2(cfg.omega0,cfg);
[~,Sig0,V0] = svd(S0);
s0 = diag(Sig0);
[~,kmin0] = min(s0);
a0 = V0(:,kmin0); a0 = a0/norm(a0);

OUT.P_fixed = zeros(Nd,1);
OUT.P_svd   = zeros(Nd,1);
OUT.Gamma_fixed = complex(zeros(Nd,1));
OUT.Gamma_svd   = complex(zeros(Nd,1));

for i=1:Nd
    S = S_TCMT_2x2(cfg.omega0 + delta(i), cfg);

    % fixed
    b = S*a0;
    OUT.P_fixed(i) = real(b'*b);
    OUT.Gamma_fixed(i) = a0'*(S*a0);

    % svd (reoptimize)
    [~,Sig,V] = svd(S);
    svals = diag(Sig);
    [smin,kmin] = min(svals);
    a = V(:,kmin); a = a/norm(a);

    OUT.P_svd(i) = smin^2;
    OUT.Gamma_svd(i) = a'*(S*a);
end
end

function RES = eval_case_outputs(delta,cfg,mismatch,mode)
Nd = numel(delta);
RES.Pout_comp1 = zeros(Nd,1);
RES.Pout_comp2 = zeros(Nd,1);
RES.Pout_tot   = zeros(Nd,1);

for i=1:Nd
    w = cfg.omega0 + delta(i);
    S = S_TCMT_2x2(w,cfg);

    if strcmpi(mode,'svd')
        [~,Sig,V] = svd(S);
        svals = diag(Sig);
        [~,kmin] = min(svals);
        a = V(:,kmin); a = a/norm(a);
    else
        error('This script uses only SVD-probe here.');
    end

    if mismatch.amp~=0, a(2) = (1+mismatch.amp)*a(2); end
    if mismatch.phase_deg~=0, a(2) = a(2)*exp(1i*mismatch.phase_deg*pi/180); end
    a = a/norm(a);

    b = S*a;
    RES.Pout_comp1(i) = abs(b(1))^2;
    RES.Pout_comp2(i) = abs(b(2))^2;
    RES.Pout_tot(i)   = RES.Pout_comp1(i) + RES.Pout_comp2(i);
end
end

function RES = eval_outputs_fixed(delta,cfg,mismatch,howPick)
% Fixed coherent setting: choose a0 at delta=0 and keep fixed across delta.

S0 = S_TCMT_2x2(cfg.omega0,cfg);

if strcmpi(howPick,'svd_at0')
    [~,Sig,V] = svd(S0);
    svals = diag(Sig);
    [~,kmin] = min(svals);
    a0 = V(:,kmin); a0 = a0/norm(a0);
else
    error('Unknown howPick.');
end

% apply mismatch once (fixed OA/EOM settings)
if mismatch.amp~=0, a0(2) = (1+mismatch.amp)*a0(2); end
if mismatch.phase_deg~=0, a0(2) = a0(2)*exp(1i*mismatch.phase_deg*pi/180); end
a0 = a0/norm(a0);

Nd = numel(delta);
RES.Pout_comp1 = zeros(Nd,1);
RES.Pout_comp2 = zeros(Nd,1);
RES.Pout_tot   = zeros(Nd,1);

for i=1:Nd
    S = S_TCMT_2x2(cfg.omega0 + delta(i), cfg);
    b = S*a0;
    RES.Pout_comp1(i) = abs(b(1))^2;
    RES.Pout_comp2(i) = abs(b(2))^2;
    RES.Pout_tot(i)   = RES.Pout_comp1(i) + RES.Pout_comp2(i);
end
end

function expo = fit_exponent(delta,Pout,winAbs)
mask = abs(delta) >= winAbs(1) & abs(delta) <= winAbs(2);
x = log(abs(delta(mask)));
y = log(max(Pout(mask),1e-30));
p = polyfit(x,y,1);
expo = p(1);
end

function panel_letter_local(ax, letter)
text(ax, -0.12, 0.98, letter, 'Units','normalized', 'FontSize',14, 'FontWeight','bold', ...
    'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');
end
