function FigureS7_nature()
% ===============================================================
% v4: Contrast-fixed TL colorbar + sparse dashed A-contours
%  - TL colorbar matches *reachable* TL range (robust quantile scaling)
%  - A contours ONLY at [0.1 0.5 0.9], dashed gray
% ===============================================================

clear; clc; close all;

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',13);

%% ---------------- Parameters (EDIT HERE) ----------------
GammaS_mag = 0.65;                 % |Gamma_S|
GammaS_ph  = 0.20*pi;              % arg(Gamma_S)
GammaS     = GammaS_mag*exp(1i*GammaS_ph);

R   = 0.35;                         % loss budget
Phi = 0.90*pi;                      % one-way phase
K   = R + 1i*Phi;
E   = exp(-2*K);                    % round-trip factor

% Sampling Gamma_L disk
nRad = 200; nAng = 420;
rL   = linspace(0,1,nRad);
aL   = linspace(-pi,pi,nAng);
[RR,TT] = meshgrid(rL,aL);
GammaL  = RR.*exp(1i*TT);

%% ---------------- Core maps ----------------
F    = GammaS.*GammaL.*E;
D    = 1 + F;
Gbar = (GammaL.*E + GammaS)./D;     % \bar{Gamma}_g

A_one = 1 - abs(Gbar).^2;           % A = 1 - |Gbar|^2

% Useful load-delivered power fraction
TL = (1-abs(GammaS)^2).*(1-abs(GammaL).^2).*abs(E)./abs(D).^2;

%% ---------------- Admissibility filter ----------------
eps_disk = 1e-10;
mask = isfinite(real(Gbar)) & isfinite(imag(Gbar)) & isfinite(TL) ...
     & (abs(Gbar) <= 1+eps_disk);

x_sc  = real(Gbar(mask));
y_sc  = imag(Gbar(mask));
TL_sc = TL(mask);

[TLmax_raw, iMax] = max(TL_sc);
zTLmax = x_sc(iMax) + 1i*y_sc(iMax);

[~, iMin] = min(abs(x_sc + 1i*y_sc));
zMin = x_sc(iMin) + 1i*y_sc(iMin);

%% ---------------- Grid on Cai–Smith disk + reachable mask ----------------
Ng = 520;
xv = linspace(-1,1,Ng);
yv = linspace(-1,1,Ng);
[X,Y] = meshgrid(xv,yv);
inside_disk = (X.^2 + Y.^2) <= 1;

Agrid = 1 - (X.^2 + Y.^2);
Agrid(~inside_disk) = NaN;

Ftl = scatteredInterpolant(x_sc,y_sc,TL_sc,'natural','none');
TLgrid = Ftl(X,Y);

shp = alphaShape(x_sc,y_sc);
shp.Alpha = 1.8*shp.Alpha;         % modest smoothing
in_reach = inShape(shp,X,Y);

TLgrid(~inside_disk | ~in_reach) = NaN;

%% ---------------- Robust TL color scaling (avoid "all blue") ----------------
% Use quantiles over reachable grid values (not the full disk)
TLvals = TLgrid(isfinite(TLgrid));
if isempty(TLvals)
    error('TLgrid has no finite values: reachable set too small.');
end

TL_lo = prctile(TLvals, 5);        % lower bound (clip near-zeros)
TL_hi = prctile(TLvals, 99.5);     % upper bound (clip rare spikes)
if TL_hi <= TL_lo
    TL_lo = min(TLvals);
    TL_hi = max(TLvals);
end

% Levels for filled contours
nLv = 8;
Lv = linspace(TL_lo, TL_hi, nLv);

%% ---------------- Resonance-tuned path Theta = pi ----------------
thetaL_res = pi - angle(GammaS) - angle(E);  % enforce Theta=pi
x_path = linspace(0,1,500);
GammaL_path = x_path .* exp(1i*thetaL_res);

D_path = 1 + GammaS.*GammaL_path.*E;
Gbar_path = (GammaL_path.*E + GammaS)./D_path;

x_opt = abs(GammaS)*exp(-2*R);
x_cc  = abs(GammaS)*exp(+2*R);

if x_opt<=1
    GammaL_opt = x_opt*exp(1i*thetaL_res);
    D_opt = 1 + GammaS*GammaL_opt*E;
    z_opt = (GammaL_opt*E + GammaS)/D_opt;
else
    z_opt = NaN;
end

%% ---------------- Figure layout ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 0 19.5 12.5]);
tlo = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

%% ===================== (a) Cai–Smith disk =====================
ax1 = nexttile(tlo,1); hold(ax1,'on'); axis(ax1,'equal'); box(ax1,'on');
xlabel(ax1,'$\Re(\bar{\Gamma}_g)$'); ylabel(ax1,'$\Im(\bar{\Gamma}_g)$');
xlim(ax1,[-1.05 1.05]); ylim(ax1,[-1.05 1.05]);

% Unit circle
th = linspace(-pi,pi,800);
plot(ax1,cos(th),sin(th),'k-','LineWidth',1.2);

% ---- TL filled contours + colorbar (2-decimal caxis) ----
contourf(ax1, X, Y, TLgrid, Lv, 'LineStyle','none');
colormap(ax1, turbo);

% Robust TL range (already computed as TL_lo, TL_hi)
% -> round caxis limits to 2 decimal places
cLo = round(TL_lo, 2);
cHi = round(TL_hi, 2);

% Safety: avoid equal limits after rounding
if cHi <= cLo
    cLo = floor(TL_lo*100)/100;
    cHi = ceil (TL_hi*100)/100;
    if cHi <= cLo
        cHi = cLo + 0.01;
    end
end

caxis(ax1, [cLo cHi]);

cb = colorbar(ax1);
cb.Label.String = '$T_L$';
cb.Label.Interpreter = 'latex';

% Set ticks and force 2-decimal tick labels
cb.Ticks = linspace(cLo, cHi, 6);
cb.TickLabels = arrayfun(@(v) sprintf('%.2f', v), cb.Ticks, 'UniformOutput', false);


% Resonance path
plot(ax1,real(Gbar_path),imag(Gbar_path),'k-','LineWidth',1.8);

% ---- Gray A contours ONLY at 0.1/0.5/0.9, dashed ----
Alevels = [0.1 0.5 0.9];
contour(ax1, X, Y, Agrid, Alevels, ...
    'LineWidth',1.0, ...
    'LineColor',[0.60 0.60 0.60], ...
    'LineStyle','--');

% ---- Place A labels near bottom-left (programmatic, non-interactive) ----
theta_lab = -3*pi/4;  % bottom-left direction (225 deg)
for a = Alevels
    r = sqrt(max(0, 1 - a));       % since A = 1 - r^2
    xlab = r*cos(theta_lab);
    ylab = r*sin(theta_lab);

    % small offset so text sits slightly inside the contour line
    xlab = xlab + 0.02;
    ylab = ylab + 0.02;

    text(ax1, xlab, ylab, sprintf('$\\mathcal{A}=%.1f$', a), ...
        'Interpreter','latex', ...
        'Color',[0.60 0.60 0.60], ...
        'FontSize',12, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'BackgroundColor','w', ...   % keep readable over TL fill
        'Margin',1);
end


% Markers
plot(ax1,real(zTLmax),imag(zTLmax),'ko','MarkerSize',6,'LineWidth',1.6);
text(ax1,real(zTLmax)+0.03,imag(zTLmax),'$T_L$ max','FontSize',12);

plot(ax1,real(zMin),imag(zMin),'ks','MarkerSize',6,'LineWidth',1.6);
text(ax1,real(zMin)+0.03,imag(zMin),'min $|\bar{\Gamma}_g|$','FontSize',12);

if ~isnan(z_opt)
    plot(ax1,real(z_opt),imag(z_opt),'kd','MarkerSize',7,'LineWidth',1.6);
    text(ax1,real(z_opt)-0.4,imag(z_opt)-0.15,'$|\Gamma_L|_{\rm opt}$ on $\Theta=\pi$','FontSize',12);
end

% title(ax1,'Cai--Smith disk: dashed $\mathcal A$ contours + filled $T_L$ + resonance path');

%% ===================== (b) 1D TL vs |GammaL| =====================
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'$|\Gamma_L|$'); ylabel(ax2,'$T_L$ (phase tuned: $\Theta=\pi$)');

x = linspace(0,1,700);
kappa = abs(GammaS).*x.*exp(-2*R);
TL_1D = (1-abs(GammaS)^2).*(1-x.^2).*exp(-2*R)./(1-kappa).^2;
plot(ax2,x,TL_1D,'k-','LineWidth',2.0);

yl = ylim(ax2);
plot(ax2,[x_opt x_opt],yl,'k--','LineWidth',1.2);
text(ax2,x_opt+0.015,yl(2)*0.95,'$|\Gamma_L|_{\rm opt}=|\Gamma_S|e^{-2R}$','FontSize',10);

if x_cc<=1
    plot(ax2,[x_cc x_cc],yl,'k:','LineWidth',1.8);
    text(ax2,x_cc+0.015,yl(2)*0.80,'$|\Gamma_L|_{\rm cc}=|\Gamma_S|e^{2R}$','FontSize',10);
else
    text(ax2,0.02,yl(2)*0.80, ...
        sprintf('$|\\Gamma_L|_{\\rm cc}=|\\Gamma_S|e^{2R}=%.2f$ (inadmissible)',x_cc), ...
        'FontSize',10);
end

% title(ax2,'Resonant tuning: $T_L$-optimum differs from critical coupling');

text(ax1, -0.2, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax2, -0.2, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

%% ---------------- Console summary ----------------
fprintf('TL color scale matched to reachable set: [%.4g, %.4g]\n', TL_lo, TL_hi);
fprintf('|GammaL|_opt = %.4f, |GammaL|_cc = %.4f, TLmax_raw = %.4f\n', x_opt, x_cc, TLmax_raw);
%% ---------- Export (stable) ----------
print(fig,'Figure S07','-dpng','-r600');
fprintf('Saved: FigureS07.png\n');

end
