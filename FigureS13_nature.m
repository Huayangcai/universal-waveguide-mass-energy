clear; close all; clc;

% ============================================================
% CPA solve (2D) + fixed@Omega0 vs per-frequency optimal input
% ------------------------------------------------------------
% K(Ω) = sqrt((iΩ+δR)(iΩ+δG)), enforce Re(K)>=0
% 1) Solve CPA by 2D root: N(Ω0)=0 for (|Γ_L|, Ω0), phase(Γ_L)=π
% 2) Local scan around Ω0
% 3) Compare TWO coherent-input strategies on every subplot:
%    (i) fixed@Ω0: use a_min(Ω0) for all Ω
%   (ii) per-frequency: use a_min(Ω) for each Ω
%
% Panels (2x2):
% (a) σmin/σmax (clipped dB) + effective ||Sa|| for fixed@Ω0
% (b) optimal coherent input ratio+phase: fixed@Ω0 vs per-frequency (NO -t/r1)
% (c) output powers P1,P2,Ptot: fixed@Ω0 vs per-frequency
% (d) scaling (right of Ω0): fixed@Ω0 vs per-frequency with separate fitted slopes
% ============================================================

%% -------------------- style --------------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultAxesFontSize',12);

fprintf('=== CPA solve (2D) + fixed@Omega0 vs per-frequency optimal input ===\n');

%% -------------------- parameters --------------------
delta_R = 1e-5;
delta_G = 1e-5;

con     = sqrt(-2 + sqrt(5));     % ≈0.485868...
Gamma_S = con * exp(1i*0);
phi_L   = pi;                      % enforce Γ_L phase pi

Omega_lo = 1e-6;
Omega_hi = 3.5;                    % allow wide; solution may land near pi or pi/2 depending on params
Omega_target = pi/2;               % initial guess helper

% local scan
dWin = 0.25;
Np   = 200000;                     % vectorized

% scaling fit window (right side)
dMin = 1e-4;
dMax = 0.025;

%% -------------------- derived --------------------
tau_S = sqrt(max(1 - abs(Gamma_S)^2, 0));

%% ============================================================
% Solve CPA: N(Ω0)=0 for (|Γ_L|, Ω0)
%% ============================================================
magL0 = con;
Omega0_guess = Omega_target;

p0 = log(magL0/(1-magL0)); % inverse sigmoid
q0 = log((Omega0_guess - Omega_lo)/(Omega_hi - Omega0_guess)); % inverse affine-sigmoid
x0 = [p0; q0];

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-14, ...
    'StepTolerance',1e-14, ...
    'OptimalityTolerance',1e-14, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',6000);

fprintf('\n--- fsolve starts ---\n');
[xsol, ~, exitflag] = fsolve(@(x) CPA_residual_2eq(x, Gamma_S, tau_S, delta_R, delta_G, Omega_lo, Omega_hi, phi_L), x0, opts);

[magL_sol, Omega0] = decode_x(xsol, Omega_lo, Omega_hi);
Gamma_L = magL_sol * exp(1i*phi_L);
tau_L   = sqrt(max(1 - magL_sol^2, 0));

N0 = N_of_Omega(Gamma_S, Gamma_L, tau_S, tau_L, delta_R, delta_G, Omega0);

fprintf('\n=== CPA solution ===\n');
fprintf('exitflag = %d\n', exitflag);
fprintf('|Gamma_L| = %.15f, phase=%.3f rad\n', magL_sol, phi_L);
fprintf('Omega0    = %.15f\n', Omega0);
fprintf('|N(Omega0)| = %.3e\n', abs(N0));

%% ============================================================
% Local scan around Ω0
%% ============================================================
Om_lo = max(Omega_lo, Omega0 - dWin);
Om_hi = min(Omega_hi, Omega0 + dWin);
Omega = linspace(Om_lo, Om_hi, Np).';   % (Np×1)
[~, idx0] = min(abs(Omega - Omega0));

% propagation
K   = K_sqrt_branch(Omega, delta_R, delta_G);
E   = exp(-2*K);
emK = exp(-K);

% denominator and S-params
D = 1 + Gamma_S*Gamma_L.*E;
D(abs(D) < 1e-14) = 1e-14 + 1i*1e-14;

r1 = (Gamma_S + Gamma_L.*E)./D;
r2 = (Gamma_L + Gamma_S.*E)./D;
t  = (emK .* tau_S .* tau_L)./D;

% N(Ω)
N = Gamma_S*Gamma_L + (Gamma_S^2 + Gamma_L^2 - (tau_S*tau_L)^2).*E + Gamma_S*Gamma_L.*E.^2;

%% ============================================================
% Singular values (2×2 closed-form via A=S'*S)
%% ============================================================
a11 = abs(r1).^2 + abs(t).^2;
a22 = abs(r2).^2 + abs(t).^2;
a12 = conj(r1).*t + conj(t).*r2;

trA  = a11 + a22;
detA = a11.*a22 - abs(a12).^2;
disc = max(trA.^2 - 4*detA, 0);
root = sqrt(disc);

lam_max = 0.5*(trA + root);
lam_min = 0.5*(trA - root);

sigma_max = sqrt(max(lam_max, 0));
sigma_min = sqrt(max(lam_min, 0));

%% ============================================================
% Per-frequency optimal input v_min(Ω) from eigenvector of A=S'*S (lam_min)
% v = [a12; lam_min - a11]
%% ============================================================
v1_pf = a12;
v2_pf = lam_min - a11;

tiny = (abs(v1_pf) + abs(v2_pf)) < 1e-30;
v1_pf(tiny) = 1; v2_pf(tiny) = 0;

nv = sqrt(abs(v1_pf).^2 + abs(v2_pf).^2);
a1_pf = v1_pf ./ nv;
a2_pf = v2_pf ./ nv;

% fixed@Ω0 input = per-frequency optimum evaluated at Ω0, then held constant
a0 = [a1_pf(idx0); a2_pf(idx0)];
a0 = a0 / norm(a0);
a1_fx = a0(1);
a2_fx = a0(2);

%% ============================================================
% (b) coherent input ratio+phase for both strategies
%% ============================================================
ratio_pf = a1_pf ./ a2_pf;
mag_pf   = abs(ratio_pf);
dphi_pf  = wrap_pi(unwrap(angle(ratio_pf)));

ratio_fx = (a1_fx / a2_fx) * ones(Np,1);
mag_fx   = abs(ratio_fx);
dphi_fx  = wrap_pi(angle(ratio_fx)) .* ones(Np,1);  % constant

%% ============================================================
% (c) Output powers using both strategies
% b = S a
%% ============================================================
% per-frequency
b1_pf = r1.*a1_pf + t.*a2_pf;
b2_pf = t.*a1_pf + r2.*a2_pf;
P1_pf = abs(b1_pf).^2;
P2_pf = abs(b2_pf).^2;
Ptot_pf = P1_pf + P2_pf;  % should equal sigma_min^2 (numerically)

% fixed@Ω0
b1_fx = r1*a1_fx + t*a2_fx;
b2_fx = t*a1_fx + r2*a2_fx;
P1_fx = abs(b1_fx).^2;
P2_fx = abs(b2_fx).^2;
Ptot_fx = P1_fx + P2_fx;

% Effective output norms (for subplot a comparison)
sigEff_pf = sqrt(max(Ptot_pf, 0));   % should ~ sigma_min
sigEff_fx = sqrt(max(Ptot_fx, 0));   % ||S a0||




%% ============================================================
% Plot 2x2 big figure (each panel compares fixed vs per-frequency)
%% ============================================================
dbFloor = -120;
sigFloor = 10^(dbFloor/20);

sigmin_dB = 20*log10(max(sigma_min, sigFloor));
sigmax_dB = 20*log10(max(sigma_max, sigFloor));
sigEff_fx_dB = 20*log10(max(sigEff_fx, sigFloor));  % effective output for fixed input

ratioFloor = 1e-12;
mag_pf_dB = 20*log10(max(mag_pf, ratioFloor));
mag_fx_dB = 20*log10(max(mag_fx, ratioFloor));

Pfloor = 1e-14;

% --- Scaling: compare both with separate fits (right of Ω0) ---
dOmega = Omega - Omega0;  % relative frequency

% Adjust dOmega to make sure it's within a valid range
dMin_eff = max(dMin, 5*(Omega(2) - Omega(1)));  % Ensure effective window
mask_fit = (dOmega > dMin_eff) & (dOmega < dMax);  % Select valid data points
if nnz(mask_fit) < 80
    warning('Fit points too few; expanding dMax to 0.10');
    dMax = 0.10;
    mask_fit = (dOmega > dMin_eff) & (dOmega < dMax);
end

dOmega = Omega - Omega0;
dMin_eff = max(dMin, 5*(Omega(2)-Omega(1)));
mask_fit = (dOmega > dMin_eff) & (dOmega < dMax);
if nnz(mask_fit) < 80
    warning('Fit points too few; expanding dMax to 0.10');
    dMax = 0.10;
    mask_fit = (dOmega > dMin_eff) & (dOmega < dMax);
end

% Fit per-frequency
px = log10(dOmega(mask_fit));
py_pf = log10(Ptot_pf(mask_fit) + 1e-300);
pf_pf = polyfit(px, py_pf, 1);
slope_pf = pf_pf(1);
int_pf   = pf_pf(2);

% Fit fixed
py_fx = log10(Ptot_fx(mask_fit) + 1e-300);
pf_fx = polyfit(px, py_fx, 1);
slope_fx = pf_fx(1);
int_fx   = pf_fx(2);

xline_ref = logspace(log10(min(dOmega(mask_fit))), log10(max(dOmega(mask_fit))), 200).';
yfit_pf = 10.^(int_pf + slope_pf*log10(xline_ref));
yfit_fx = 10.^(int_fx + slope_fx*log10(xline_ref));

% Update: Independent reference line positioning for Scaling
% Adjusted position based on the geometric mean of dOmega and the power range
x_ref = sqrt(min(dOmega(mask_fit)) * max(dOmega(mask_fit)));
y_ref = sqrt(min([Ptot_pf(mask_fit); Ptot_fx(mask_fit)]) * ...
             max([Ptot_pf(mask_fit); Ptot_fx(mask_fit)]));
y_ref = y_ref * 0.3; % Adjust scale

% Reference lines for ΔΩ^2 and ΔΩ^4
yref4 = y_ref * (xline_ref./x_ref).^4;  % ΔΩ⁴ reference line
yref2 = y_ref * (xline_ref./x_ref).^2;  % ΔΩ² reference line

fprintf('\n=== Scaling fit (right of Ω0) ===\n');
fprintf('Fit window: dOmega in (%.3e, %.3e), points=%d\n', dMin_eff, dMax, nnz(mask_fit));
fprintf('slope (per-frequency) = %.4f\n', slope_pf);
fprintf('slope (fixed@Ω0)      = %.4f\n', slope_fx);
% Plotting
width_cm = 22; height_cm = 18;
fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% (a) Singular values
ax1 = nexttile; hold on; box on; grid on;
yyaxis left
plot(Omega, sigmin_dB, 'b-', 'LineWidth', 2.2, 'DisplayName','$\sigma_{\min}$');
plot(Omega, sigEff_fx_dB, 'k--', 'LineWidth', 2.0, 'DisplayName','$\|\mathbf{S} \mathbf{a}_\mathbf{fixed}\|$');
ylabel('$20\log_{10}(\cdot)$ (dB)');
yyaxis right
plot(Omega, sigmax_dB, 'r--', 'LineWidth', 2.2, 'DisplayName','$\sigma_{\max}$');
ylabel('$20\log_{10}\sigma_{\max}$ (dB)');
xline(Omega0,'k--','LineWidth',1.6,'HandleVisibility','off');
xlabel('$\Omega$');
xlim([min(Omega), max(Omega)]);
% title('(a) Singular values and effective output norm');
legend('Location','southwest'); legend boxoff
ylim(ax1, [dbFloor, 5]);

% (b) Optimal coherent input: fixed@Ω0 vs per-frequency (NO -t/r1)
ax2 = nexttile; hold on; box on; grid on;
yyaxis left
plot(Omega, mag_pf_dB, 'k-', 'LineWidth', 2.0, 'DisplayName','$20\log_{10}|\mathbf{a_1}/\mathbf{a_2}|$ (per-freq)');
plot(Omega, mag_fx_dB, 'k--','LineWidth', 2.0, 'DisplayName','$20\log_{10}|\mathbf{a_1}/\mathbf{a_2}|$ (fixed@$\Omega_0$)');
ylabel('$20\log_{10}||\mathbf{a_1}/\mathbf{a_2}||$ (dB)');
yyaxis right
plot(Omega, dphi_pf, 'm-', 'LineWidth', 2.0, 'DisplayName','$\Delta\phi$ (per-freq)');
plot(Omega, dphi_fx, 'm--','LineWidth', 2.0, 'DisplayName','$\Delta\phi$ (fixed@$\Omega_0$)');
ylabel('$\Delta\psi$ (rad)');
xline(Omega0,'k--','LineWidth',1.6,'HandleVisibility','off');
xlabel('$\Omega$');
xlim([min(Omega), max(Omega)]);
% title('(b) Optimal coherent input: per-frequency vs fixed');
legend('Location','southwest'); legend boxoff

% (c) Output powers: per-frequency vs fixed@Ω0
ax3 = nexttile; hold on; box on; grid on;
semilogy(Omega, max(Ptot_pf, Pfloor), 'k-',  'LineWidth', 2.4, 'DisplayName','$P_{\rm tot}$ (per-freq)');
semilogy(Omega, max(Ptot_fx, Pfloor), 'k--', 'LineWidth', 2.0, 'DisplayName','$P_{\rm tot}$ (fixed@$\Omega_0$)');

semilogy(Omega, max(P1_pf, Pfloor), 'b-',  'LineWidth', 1.8, 'DisplayName','$P_1$ (per-freq)');
semilogy(Omega, max(P1_fx, Pfloor), 'b--', 'LineWidth', 1.6, 'DisplayName','$P_1$ (fixed)');

semilogy(Omega, max(P2_pf, Pfloor), 'r:',  'LineWidth', 1.8, 'DisplayName','$P_2$ (per-freq)');
semilogy(Omega, max(P2_fx, Pfloor), 'r--', 'LineWidth', 1.6, 'DisplayName','$P_2$ (fixed)');

xline(Omega0,'k--','LineWidth',1.6,'HandleVisibility','off');
xlabel('$\Omega$');
ylabel('$P_{\rm out}$');
xlim([min(Omega), max(Omega)]);
% title('(c) Output powers: per-frequency vs fixed@$\Omega_0$');
legend('Location','best'); legend boxoff
ylim([Pfloor, max([Ptot_pf;Ptot_fx])*1.2]);

% (d) Scaling: compare both with separate fits
ax4 = nexttile; hold on; box on; grid on;
loglog(dOmega(mask_fit), Ptot_pf(mask_fit), 'bo', 'MarkerSize', 4.5, 'DisplayName','data (per-freq)');
loglog(dOmega(mask_fit), Ptot_fx(mask_fit), 'ko', 'MarkerSize', 3.8, 'DisplayName','data (fixed)');

loglog(xline_ref, yfit_pf, 'r-', 'LineWidth', 2.2, 'DisplayName', sprintf('fit per-freq slope=%.2f', slope_pf));
loglog(xline_ref, yfit_fx, 'g-', 'LineWidth', 2.2, 'DisplayName', sprintf('fit fixed slope=%.2f', slope_fx));

% loglog(xline_ref, yref4, 'r--', 'LineWidth', 1.8, 'DisplayName','ref $\propto \Delta\Omega^4$');
% loglog(xline_ref, yref2, 'g--', 'LineWidth', 1.8, 'DisplayName','ref $\propto \Delta\Omega^2$');

xlabel('$\Delta\Omega=\Omega-\Omega_0$');
ylabel('$P_{\rm tot}$');
% title('(d) Scaling (right of $\Omega_0$): per-freq vs fixed');
legend('Location','best'); legend boxoff

text(ax1, -0.14, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax2, -0.14, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax3, -0.14, 0.98, 'c', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');

text(ax4, -0.14, 0.98, 'd', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');
% sgtitle(sprintf('Compare coherent inputs: $\\Omega_0=%.6f$, $|N(\\Omega_0)|=%.2e$', Omega0, abs(N0)), ...
%     'FontSize',13);

% print(fig, '-dpng', '-r600', 'CPAEP_fixed_vs_perfreq.png');
% fprintf('\nSaved: CPAEP_fixed_vs_perfreq.png\n');
% 保存图形为PNG文件
saveas(fig, 'FigureS13.png', 'png');
fprintf('\nSaved: FigureS13.png\n');


fprintf('\n=== DONE ===\n');




%% ============================================================
%                    Local functions
% ============================================================

function F = CPA_residual_2eq(x, Gamma_S, tau_S, delta_R, delta_G, Om_lo, Om_hi, phi_L)
    [magL, Omega] = decode_x(x, Om_lo, Om_hi);
    Gamma_L = magL * exp(1i*phi_L);
    tau_L   = sqrt(max(1 - magL^2, 0));
    N = N_of_Omega(Gamma_S, Gamma_L, tau_S, tau_L, delta_R, delta_G, Omega);
    F = [real(N); imag(N)];
end

function [magL, Omega] = decode_x(x, Om_lo, Om_hi)
    p = x(1);
    q = x(2);
    s1 = sigmoid(p);
    s2 = sigmoid(q);
    epsm = 1e-10;
    magL  = epsm + (1-2*epsm)*s1;
    Omega = Om_lo + (Om_hi-Om_lo)*s2;
end

function y = sigmoid(z)
    y = 1 ./ (1 + exp(-z));
end

function K = K_sqrt_branch(Omega, delta_R, delta_G)
    Z = (1i*Omega + delta_R) .* (1i*Omega + delta_G);
    K = sqrt(Z);
    flip = real(K) < 0;
    K(flip) = -K(flip);
end

function N = N_of_Omega(Gamma_S, Gamma_L, tau_S, tau_L, delta_R, delta_G, Omega)
    K = K_sqrt_branch(Omega, delta_R, delta_G);
    E = exp(-2*K);
    N = Gamma_S*Gamma_L + (Gamma_S^2 + Gamma_L^2 - (tau_S*tau_L)^2).*E + Gamma_S*Gamma_L.*E.^2;
end

function y = wrap_pi(x)
    y = mod(x + pi, 2*pi) - pi;
end
