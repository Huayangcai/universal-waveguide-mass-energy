function FigureS5_nature()
% ===============================================================
% Figure S5: Resonance-aware maximum absorbed power principle at fixed attenuation
% (1x2 figure, Nature-style)
%
% Theory (Supplementary Methods S7, resonance direction under D = 1 + F):
%   K(Ω) = R(Ω) + iΦ(Ω),  E = exp(-2K) = exp(-2R) exp(-i2Φ)
%   Γ_g = (Γ_L E + Γ_S) / (1 + Γ_S Γ_L E)          (Möbius boundary composition)
%   Eq. (S7-17):  P_abs = exp(-4 ξ R) - |Γ_g|^2 exp(4(1-ξ)R)
%
% Resonant (constructive) interference class (D = 1 + F):
%   Θ = arg(Γ_S Γ_L) - 2Φ  ≈ (2m+1)π  -> here we enforce Θ = π at each R.
%
% Critical coupling (perfect cancellation of the numerator => Γ_g = 0):
%   Γ_S = - Γ_L E  (at Ω = Ω_res)
%   => amplitude: |Γ_S| = |Γ_L| exp(-2R)  <=>  |Γ_L|_cc = |Γ_S| exp(+2R)
%   => phase:     arg Γ_S = arg Γ_L - 2Φ + π  (mod 2π)
%
% Practical note:
%   Since passive loads satisfy |Γ_L| ≤ 1, critical coupling is achievable
%   only for R ≤ R_cc = 0.5 ln(1/|Γ_S|). Beyond this, the best attainable
%   point saturates at |Γ_L| = 1 (with phase still on the resonant class).
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
gS    = 0.1;     % |Gamma_S|
xi    = 0.50;     % asymmetry coefficient ξ ∈ [0,1]
phi_S = 0;        % arg(Gamma_S)
Phi0  = 0;        % set Φ=0 for the resonant slice (Θ enforced by choosing phi_L)
Rmin  = 0; 
Rmax  = 2;        % range of one-way attenuation R = Re[K(Ω)]

% Resonant (odd-π) direction under D=1+F: Θ = π  => choose phi_L accordingly:
% Θ = (phi_S + phi_L) - 2Phi0 = π  (mod 2π)
phi_L_res = wrapToPi(pi + 2*Phi0 - phi_S);  % enforce Θ=π for the plotted slice

% Critical-coupling feasibility boundary (|Gamma_L|_cc <= 1)
Rcc = 0.5*log(1/max(gS,1e-12));

%% -------------------- Figure + tiledlayout --------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 23.5 14]);
tlo = tiledlayout(fig, 1, 2, 'Padding','compact','TileSpacing','compact');

%% ===============================================================
% (a) Resonant absorbed-power map over (R, |Gamma_L|)
% Using Eq. (S7-17) with Γ_g computed from boundary composition.
% We mask regions where Eq. (S7-17) yields non-physical negative values
% (i.e., where it would violate passivity for the chosen slice).
% ===============================================================
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on');

NR  = 520;
NGL = 420;
R   = linspace(Rmin, Rmax, NR);

gL_max = 0.9995;
gL     = linspace(0, gL_max, NGL);
[RR, GL] = meshgrid(R, gL);

% Propagation factor at the chosen resonant slice
E = exp(-2*RR) .* exp(-1i*2*Phi0);   % E = exp(-2K) with Φ=Phi0

% Complex boundary reflections on the resonant (odd-π) class
Gamma_S = gS * exp(1i*phi_S);
Gamma_L = GL .* exp(1i*phi_L_res);

% Boundary-composed generalized reflection Γ_g
Gamma_g = (Gamma_L .* E + Gamma_S) ./ (1 + Gamma_S .* Gamma_L .* E);

% Eq. (S7-17): absorbed-power metric
P_abs = exp(-4*xi*RR) - abs(Gamma_g).^2 .* exp(4*(1-xi)*RR);

% Mask negative values instead of clamping (clearer & theory-consistent)
maskPhys = (P_abs > 0) & isfinite(P_abs);
P_plot   = nan(size(P_abs));
P_plot(maskPhys) = P_abs(maskPhys);

% Plot as log10 for dynamic range
Z = nan(size(P_plot));
Z(maskPhys) = log10(P_plot(maskPhys));

hImg = imagesc(ax1, R, gL, Z);
set(ax1,'YDir','normal');

% Alpha mask: show only physically valid region
alphaImg = zeros(size(Z));
alphaImg(maskPhys) = 1;
set(hImg,'AlphaData',alphaImg);

% Colormap
try colormap(ax1, turbo(256)); catch colormap(ax1, parula(256)); end

% Robust caxis
Zvec = Z(maskPhys);
if ~isempty(Zvec)
    cl = prctile(Zvec, [2 98]);
    if cl(1) == cl(2)
        cl = [cl(1)-1, cl(2)+1];
    end
    caxis(ax1, cl);
end

cb = colorbar(ax1,'eastoutside');
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.Label.String = '$\log_{10} P_{\rm abs}(R,|\Gamma_L|)$';

xlabel(ax1,'$R=\Re[K(\Omega)]$');
ylabel(ax1,'$|\Gamma_L(\Omega)|$');

% --- Critical coupling line: |Gamma_L|_cc = |Gamma_S| exp(+2R)
gL_cc = gS * exp(2*R);

% Plot only the feasible segment where gL_cc <= 1
idxFeas = (gL_cc <= 1);
h_cc = plot(ax1, R(idxFeas), gL_cc(idxFeas), 'w--', 'LineWidth', 3.0);

% Mark the feasibility boundary (where |Gamma_L|_cc hits 1)
if Rcc >= Rmin && Rcc <= Rmax
    xline(ax1, Rcc, 'w:', 'LineWidth', 2.0, 'HandleVisibility','off');
    text(ax1, Rcc, 0.97, '$|\Gamma_L|_{\rm cc}=1$', ...
        'Color','w','Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'BackgroundColor',[0 0 0 0.25],'EdgeColor','none');
end

% Limiting references
h_lim0 = plot(ax1, R, gS*ones(size(R)), '--', 'Color',[0.25 0.65 1.0], 'LineWidth', 2.2);
h_sat  = plot(ax1, R, ones(size(R)), '--', 'Color',[1.00 0.45 0.15], 'LineWidth', 2.2);

xlim(ax1,[Rmin Rmax]); ylim(ax1,[0 1]);
grid(ax1,'on'); grid(ax1,'minor');

legend(ax1, [h_cc h_lim0 h_sat], ...
    { ...
    '$|\Gamma_L|_{\rm cc}=|\Gamma_S|e^{2R}$', ...
    '$R\to0:\ |\Gamma_L|_{\rm cc}\to|\Gamma_S|$', ...
    '$|\Gamma_L|=1$' ...
    }, ...
    'Location','east','Box','off','Interpreter','latex');

% Phase/interference-class note (this panel is on Θ=π slice)
% text(ax1, 0.52, 0.12, sprintf('$\\Theta=\\pi$ (odd-$\\pi$ constructive class), $\\phi_S=%.0f^\\circ$, $\\phi_L=%.0f^\\circ$', ...
%     phi_S*180/pi, phi_L_res*180/pi), ...
%     'Units','normalized','FontSize',11, ...
%     'BackgroundColor',[1 1 1 0.70], 'EdgeColor','k', ...
%     'Interpreter','latex','HorizontalAlignment','center');

% Panel label
text(ax1, -0.18, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

%% ===============================================================
% (b) Maximum absorbed power vs R (with passive constraint |Gamma_L|<=1)
%   - Ideal (if critical coupling feasible): P_max = exp(-4 ξ R) (Γ_g=0)
%   - Constrained maximum: choose |Γ_L| = min(|Γ_S|e^{2R}, 1) on Θ=π slice
%   - Also show "amplitude-only but wrong phase" to illustrate phase requirement
% ===============================================================
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on');

Rline = linspace(max(Rmin,1e-6), Rmax, 900);

% Ideal maximum (critical coupling, only meaningful where feasible)
P_ideal = exp(-4*xi*Rline);
h_ideal = plot(ax2, Rline, P_ideal, 'k-', 'LineWidth', 3.0);
h_ideal.DisplayName = '$P_{\rm abs,max}^{\rm (cc)}(R)=e^{-4\xi R}$ (if feasible)';

% Constrained maximum with |Gamma_L|<=1 on the resonant slice Θ=π
gL_star = min(gS*exp(2*Rline), 1);
Gamma_L_star = gL_star .* exp(1i*phi_L_res);
Eline = exp(-2*Rline) .* exp(-1i*2*Phi0);
Gamma_g_star = (Gamma_L_star .* Eline + Gamma_S) ./ (1 + Gamma_S .* Gamma_L_star .* Eline);
P_star = exp(-4*xi*Rline) - abs(Gamma_g_star).^2 .* exp(4*(1-xi)*Rline);
P_star = max(P_star, 0);  % constrain to passive range for plotting
h_star = plot(ax2, Rline, P_star, '-', 'Color',[0.15 0.55 0.20], 'LineWidth', 2.6);
h_star.DisplayName = 'Constrained optimum ($|\Gamma_L|\le1$, $\Theta=\pi$)';

% Wrong-phase demonstration (same amplitude choice, but set Θ=0 by choosing phi_L = -phi_S + 2Phi0)
phi_L_wrong = wrapToPi(0 + 2*Phi0 - phi_S); % enforce Θ=0 (destructive class under D=1+F)
Gamma_L_wrong = gL_star .* exp(1i*phi_L_wrong);
Gamma_g_wrong = (Gamma_L_wrong .* Eline + Gamma_S) ./ (1 + Gamma_S .* Gamma_L_wrong .* Eline);
P_wrong = exp(-4*xi*Rline) - abs(Gamma_g_wrong).^2 .* exp(4*(1-xi)*Rline);
P_wrong = max(P_wrong, 0);
h_wrong = plot(ax2, Rline, P_wrong, '--', 'Color',[0.80 0.20 0.20], 'LineWidth', 2.2);
h_wrong.DisplayName = 'Amplitude-matched but wrong phase ($\Theta=0$)';

% Feasibility boundary line
if Rcc >= Rmin && Rcc <= Rmax
    xline(ax2, Rcc, 'k:', 'LineWidth', 1.8, 'HandleVisibility','off');
    text(ax2, Rcc, 0.26, '$R_{\rm cc}=\frac{1}{2}\ln(1/|\Gamma_S|)$', ...
        'Interpreter','latex','FontSize',11, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'BackgroundColor',[1 1 1 0.70],'EdgeColor','none');
end

% Small-R expansion guide for the ideal curve
P_taylor = 1 - 4*xi*Rline;
h_taylor = plot(ax2, Rline, P_taylor, ':', 'Color',[0.25 0.65 1.0], 'LineWidth', 2.0);
h_taylor.DisplayName = '$R\to0:\ e^{-4\xi R}\approx 1-4\xi R$';

grid(ax2,'on'); grid(ax2,'minor');
xlim(ax2,[Rmin Rmax]);
ylim(ax2,[0 1.05]);

xlabel(ax2,'$R=\Re[K(\Omega)]$');
ylabel(ax2,'$P_{\rm abs}$');

legend(ax2,'Location','northwest','Box','off','Interpreter','latex','FontSize',11);

% text(ax2, 0.05, 0.83, sprintf('$|\\Gamma_S|=%.2f,\\ \\xi=%.2f$', gS, xi), ...
%     'Units','normalized','FontSize',13, ...
%     'BackgroundColor',[1 1 1 0.70],'EdgeColor','none','Interpreter','latex');

% text(ax2, 0.05, 0.72, sprintf('Resonant class: $\\Theta=\\pi$ (set by $\\phi_L=%.0f^\\circ$)', phi_L_res*180/pi), ...
%     'Units','normalized','FontSize',11, ...
%     'BackgroundColor',[1 1 1 0.70],'EdgeColor','k','Interpreter','latex');

% Panel label
text(ax2, -0.18, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

%% -------------------- Export --------------------
% set(findall(gcf,'-property','FontSize'),'FontSize',14);
filename = 'Figure S5';
print(fig,'-dpng','-r600',[filename,'.png']);
fprintf('Saved: %s.png\n', filename);

end
