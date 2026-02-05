function FigureS12_nature()
% ==============================================================
% Solve_z_and_verify_TCMT_vs_Waveguide_v2
%
% Fixes vs v1:
%   (1) Do NOT wrap v into [-pi,pi] during optimization (prevents jumps).
%   (2) Track u(omega), v(omega) directly; compute K_eff continuously:
%         z = exp(-u + i v),  K_eff = u/2 - i*(unwrap(v)/2)
%       (avoid principal-log branch jumps).
%   (3) Use sqrtz = exp(-u/2 + i v/2) consistently everywhere (avoid sqrt branch).
%   (4) Plot log10|z| (or semilogy |z|) so it's visible when |z|<<1.
%   (5) Optional: verify CPA input vector consistency at omega=0.
%
% ==============================================================

clc; close all;
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultAxesFontSize',12);
set(groot,'DefaultLineLineWidth',1.8);

%% ------------------ TCMT parameters (example) ------------------
par = struct();
par.w01 = 0;
par.w02 = 0;

par.g1  = 800;
par.g2  = 1600;

par.gc1 = 1800;
par.gc2 = (par.g1 + par.g2) - par.gc1;     % real-zero balance

par.kappa = abs(par.g1 - par.gc1)/2;       % CPA-EP coupling condition (w01=w02)

fprintf('TCMT params:\n');
fprintf('  g1=%g, g2=%g, gc1=%g, gc2=%g, kappa=%g\n', par.g1, par.g2, par.gc1, par.gc2, par.kappa);
fprintf('  balance check: gc1+gc2=%.3f vs g1+g2=%.3f\n', par.gc1+par.gc2, par.g1+par.g2);

%% ------------------ frequency grid ------------------
wMax = 4000;
Nw   = 1001;
w    = linspace(-wMax, wMax, Nw).';

%% ------------------ effective coupling amplitudes (absorbing I into boundaries) ------------------
tauS = sqrt(max(par.gc1,0));
tauL = sqrt(max(par.gc2,0));
tauProd = tauS*tauL;

%% ------------------ solve z(w) and reconstruct ------------------
zSol   = complex(zeros(Nw,1));
uSol   = zeros(Nw,1);     % z = exp(-u + i v)
vSol   = zeros(Nw,1);
Dsol   = complex(zeros(Nw,1));
GSsol  = complex(zeros(Nw,1));
GLsol  = complex(zeros(Nw,1));

errFrob = zeros(Nw,1);
errMax  = zeros(Nw,1);

% warm start: optimize p=[log(u); v] with u=exp(p1), v free (NOT wrapped)
pPrev = [log(1.0); 0.0];    % u=1 => |z|=exp(-1)=0.367; v=0

opts = optimset('Display','off', 'TolX',1e-11, 'TolFun',1e-13, ...
                'MaxIter',3000, 'MaxFunEvals',9000);

for k = 1:Nw
    wk = w(k);

    S = S_TCMT_2x2(wk, par);
    r1 = S(1,1);
    r2 = S(2,2);
    t  = S(1,2);

    [z, u, v, D, GS, GL, pPrev] = solve_z_onefreq(r1, r2, t, tauProd, pPrev, opts);

    zSol(k)  = z;
    uSol(k)  = u;
    vSol(k)  = v;
    Dsol(k)  = D;
    GSsol(k) = GS;
    GLsol(k) = GL;

    Shat = S_WG_from_uv(u, v, D, GS, GL, tauProd);

    E = Shat - S;
    errFrob(k) = norm(E,'fro');
    errMax(k)  = max(abs(E(:)));
end

fprintf('\nReconstruction error summary over grid:\n');
fprintf('  max Frobenius error: %.3e\n', max(errFrob));
fprintf('  max entrywise error: %.3e\n', max(errMax));

%% ------------------ Solve z(w) and reconstruct ------------------
zSol   = complex(zeros(Nw,1));
uSol   = zeros(Nw,1);     % z = exp(-u + i v)
vSol   = zeros(Nw,1);
Dsol   = complex(zeros(Nw,1));
GSsol  = complex(zeros(Nw,1));
GLsol  = complex(zeros(Nw,1));

errFrob = zeros(Nw,1);
errMax  = zeros(Nw,1);

% warm start: optimize p=[log(u); v] with u=exp(p1), v free (NOT wrapped)
pPrev = [log(1.0); 0.0];    % u=1 => |z|=exp(-1)=0.367; v=0

opts = optimset('Display','off', 'TolX',1e-11, 'TolFun',1e-13, ...
                'MaxIter',3000, 'MaxFunEvals',9000);

for k = 1:Nw
    wk = w(k);

    S = S_TCMT_2x2(wk, par);
    r1 = S(1,1);
    r2 = S(2,2);
    t  = S(1,2);

    [z, u, v, D, GS, GL, pPrev] = solve_z_onefreq(r1, r2, t, tauProd, pPrev, opts);

    zSol(k)  = z;
    uSol(k)  = u;
    vSol(k)  = v;
    Dsol(k)  = D;
    GSsol(k) = GS;
    GLsol(k) = GL;

    Shat = S_WG_from_uv(u, v, D, GS, GL, tauProd);

    E = Shat - S;
    errFrob(k) = norm(E,'fro');
    errMax(k)  = max(abs(E(:)));
end

% Now extract the parameters at omega = 0 (first element of w)
omega0_index = find(w == 0);
fprintf('\nAt omega=0 (w = 0):\n');
fprintf('u(0) = %.6f, v(0) = %.6f, D(0) = %.6f\n', uSol(omega0_index), vSol(omega0_index), Dsol(omega0_index));
fprintf('GS(0) = %.6f, GL(0) = %.6f, tauProd = %.6f\n', GSsol(omega0_index), GLsol(omega0_index), tauProd);

fprintf('\nReconstruction error summary over grid:\n');
fprintf('  max Frobenius error: %.3e\n', max(errFrob));
fprintf('  max entrywise error: %.3e\n', max(errMax));



%% ------------------ continuous K_eff from (u,v) ------------------
% z = exp(-u + i v) => K_eff = -0.5 log z
% Choose continuous branch by unwrapping v across omega:
vUn = unwrap(vSol);
Keff = (uSol/2) - 1i*(vUn/2);

%% ------------------ optional: CPA input vector consistency at omega=0 ------------------
% TCMT CPA vector (up to global phase): [± i sqrt(gc1); sqrt(gc2)]
% WG CPA vector from v_min (right singular vector of S_WG at omega0)
omega0 = 0;
S0 = S_TCMT_2x2(omega0, par);
[rU, ~, rV] = svd(S0);
[~, idxMin] = min(diag(rU)); %#ok<ASGLU>  % not used; keep to show idea
[~, Sig0, V0] = svd(S0); svals0 = diag(Sig0);
[~,kmin] = min(svals0);
vmin0 = V0(:,kmin); vmin0 = vmin0/norm(vmin0);

aCPA_tcmt = [1i*sqrt(par.gc1); sqrt(par.gc2)];
aCPA_tcmt = aCPA_tcmt/norm(aCPA_tcmt);

% compare phase-insensitive overlap |<vmin0, aCPA_tcmt>|
ov = abs(vmin0' * aCPA_tcmt);
fprintf('\nCPA input check at omega=0:\n');
fprintf('  |< v_min_TCMT(0) , a_CPA_TCMT >| = %.6f (1 means identical up to global phase)\n', ov);
fprintf('  v_min_TCMT ratio a1/a2 = %.4g%+.4gi\n', real(vmin0(1)/vmin0(2)), imag(vmin0(1)/vmin0(2)));
fprintf('  a_CPA_TCMT ratio a1/a2 = %.4g%+.4gi\n', real(aCPA_tcmt(1)/aCPA_tcmt(2)), imag(aCPA_tcmt(1)/aCPA_tcmt(2)));
% 
% Initialize arrays to store singular values
singular_values_S = zeros(Nw, 2); % For the TCMT scattering matrix
singular_values_Shat = zeros(Nw, 2); % For the waveguide model (Shat)

for k = 1:Nw
    wk = w(k);
    
    % Compute the scattering matrix for TCMT at frequency wk
    S = S_TCMT_2x2(wk, par);
    
    % Compute the singular values of the TCMT scattering matrix S
    [~, Sig_S, ~] = svd(S);
    singular_values_S(k, :) = diag(Sig_S);  % Store the singular values
    
    % Compute the reconstructed scattering matrix Shat from the waveguide model
    Shat = S_WG_from_uv(uSol(k), vSol(k), Dsol(k), GSsol(k), GLsol(k), tauProd);
    
    % Compute the singular values of the reconstructed scattering matrix Shat
    [~, Sig_Shat, ~] = svd(Shat);
    singular_values_Shat(k, :) = diag(Sig_Shat);  % Store the singular values
end

% 提取 omega = 0 对应的索引（假设 w 包含精确 0）
omega0_index = find(abs(w) < 1e-10, 1);   % 更鲁棒写法

fprintf('\n=== ω = 0 处主要变量 (第一张图对应点) ===\n');
fprintf('u(0)          = %.10f\n', uSol(omega0_index));
fprintf('v(0)          = %.10f rad\n', vSol(omega0_index));
fprintf('z(0)          = %.10f %+10.10fi\n', real(zSol(omega0_index)), imag(zSol(omega0_index)));
fprintf('D(0)          = %.10f %+10.10fi\n', real(Dsol(omega0_index)), imag(Dsol(omega0_index)));
fprintf('Γ_S^eff(0)    = %.10f %+10.10fi\n', real(GSsol(omega0_index)), imag(GSsol(omega0_index)));
fprintf('Γ_L^eff(0)    = %.10f %+10.10fi\n', real(GLsol(omega0_index)), imag(GLsol(omega0_index)));
fprintf('tauProd       = %.6f\n', tauProd);
fprintf('K_eff(0)      = %.10f %+10.10fi\n', real(Keff(omega0_index)), imag(Keff(omega0_index)));

% CPA 输入一致性（如果已计算）
fprintf('CPA overlap   = %.10f\n', ov);
fprintf('v_min ratio   = %.6f %+6.6fi\n', real(vmin0(1)/vmin0(2)), imag(vmin0(1)/vmin0(2)));
fprintf('理论 CPA ratio= %.6f %+6.6fi\n', real(aCPA_tcmt(1)/aCPA_tcmt(2)), imag(aCPA_tcmt(1)/aCPA_tcmt(2)));

%% Plot Singular Values for TCMT and Waveguide Model
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 18]);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
nTiles = 6;
letters = char('a' + (0:(nTiles-1)));  % a..l
% (a) Singular values of the TCMT model and waveguide model comparison
ax1=nexttile; hold on; box on; grid on;
plot(w, singular_values_S(:, 1), 'LineWidth', 2); % Singular value 1 for TCMT model
plot(w, singular_values_Shat(:, 1), '--', 'LineWidth', 2); % Singular value 1 for Waveguide model
plot(w, singular_values_S(:, 2), 'LineWidth', 2); % Singular value 2 for TCMT model
plot(w, singular_values_Shat(:, 2), '--', 'LineWidth', 2); % Singular value 2 for Waveguide model
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
xlabel('$\omega$');
ylabel('Singular Values');
% title('(a) Singular Values Comparison: S vs Shat');
legend1=legend({'$\sigma_{max}$ (TCMT)', '$\sigma_{max}$ (Waveguide)', '$\sigma_{min}$ (TCMT)', '$\sigma_{min}$ (Waveguide)'}, 'Location', 'sw');
legend boxoff
set(legend1,...
    'Position',[0.304086156302872 0.742980215835156 0.180408592156842 0.111356211643593]);
panel_letter_local(ax1, letters(1));
% Optional: You can add more subplots to visualize the singular value behavior in more detail


% (b) Logarithmic magnitude and phase of z (dual y-axes)
ax2=nexttile; hold on; box on; grid on;
[ax, h1, h2] = plotyy(w, log10(max(abs(zSol), 1e-300)), w, unwrap(angle(zSol)));
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
set(h1, 'LineWidth', 1.8);
set(h2, 'LineWidth', 1.8);
xlabel('$\omega$');
ylabel(ax(1), '$\log_{10}|z|$');
ylabel(ax(2), '$\arg(z)$');
% title('(b) Logarithmic magnitude and phase of $z(\omega)$');
legend({'$\log_{10}|z|$', 'Phase of $z$'}, 'Location', 'best');
legend boxoff
panel_letter_local(ax2, letters(2));
% (c) Real and imaginary parts of K_eff (dual y-axes)
ax3=nexttile; hold on; box on; grid on;
[ax, h1, h2] = plotyy(w, real(Keff), w, imag(Keff));
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
set(h1, 'LineWidth', 1.8);
set(h2, 'LineWidth', 1.8);
xlabel('$\omega$');
ylabel(ax(1), '$\mathrm{Re}(K_{\mathrm{eff}})$');
ylabel(ax(2), '$\mathrm{Im}(K_{\mathrm{eff}})$');
% title('(c) Real and Imaginary parts of $K_{\mathrm{eff}}$');
legend({'$\mathrm{Re}(K_{\mathrm{eff}})$', '$\mathrm{Im}(K_{\mathrm{eff}})$'}, 'Location', 'best');
legend boxoff
panel_letter_local(ax3, letters(3));
% (d) Effective boundary magnitude and phase for $\Gamma_S^{eff}$ (dual y-axes)
ax4=nexttile; hold on; box on; grid on;
[ax, h1, h2] = plotyy(w, abs(GSsol), w, unwrap(angle(GSsol)));
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
set(h1, 'LineWidth', 1.8);
set(h2, 'LineWidth', 1.8);
xlabel('$\omega$');
ylabel(ax(1), '$|\Gamma_S^{\mathrm{eff}}|$');
ylabel(ax(2), '$\arg(\Gamma_S^{\mathrm{eff}})$');
% title('(d) Effective boundary for $\Gamma_S^{\mathrm{eff}}$');
legend({'$|\Gamma_S^{\mathrm{eff}}|$', '$\arg(\Gamma_S^{\mathrm{eff}})$'}, 'Location', 'best');
legend boxoff
panel_letter_local(ax4, letters(4));
% (e) Effective boundary magnitude and phase for $\Gamma_L^{eff}$ (dual y-axes)
ax5=nexttile; hold on; box on; grid on;
[ax, h1, h2] = plotyy(w, abs(GLsol), w, unwrap(angle(GLsol)));
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
set(h1, 'LineWidth', 1.8);
set(h2, 'LineWidth', 1.8);
xlabel('$\omega$');
ylabel(ax(1), '$|\Gamma_L^{\mathrm{eff}}|$');
ylabel(ax(2), '$\arg(\Gamma_L^{\mathrm{eff}})$');
% title('(e) Effective boundary for $\Gamma_L^{\mathrm{eff}}$');
legend({'$|\Gamma_L^{\mathrm{eff}}|$', '$\arg(\Gamma_L^{\mathrm{eff}})$'}, 'Location', 'best');
legend boxoff
panel_letter_local(ax5, letters(5));
% (f) Effective coupling $\kappa_{\mathrm{eff}}$ and phase difference (dual y-axes)
ax6=nexttile; hold on; box on; grid on;
kappaEff = abs(GSsol .* GLsol) .* exp(-2 * real(Keff)); % Compute kappa_eff
phaseDiff = unwrap(angle(GSsol .* GLsol)) - 2 * imag(Keff); % Compute phase difference
[ax, h1, h2] = plotyy(w, kappaEff, w, phaseDiff);
xline(0, '--k', 'LineWidth', 1.5,'HandleVisibility','off'); % Add vertical dashed line at omega=0
set(h1, 'LineWidth', 1.8);
set(h2, 'LineWidth', 1.8);
xlabel('$\omega$');
ylabel(ax(1), '$\kappa_{\mathrm{eff}}$');
ylabel(ax(2), '$\arg(\Gamma_S^{\mathrm{eff}} \cdot \Gamma_L^{\mathrm{eff}}) - 2 \mathrm{Im}(K_{\mathrm{eff}})$');
% title('(f) Effective coupling $\kappa_{\mathrm{eff}}$ and phase difference');
legend({'$\kappa_{\mathrm{eff}}$', 'Phase difference'}, 'Location', 'south');
legend boxoff
panel_letter_local(ax6, letters(6));

%Save the figure as PNG
print(fig, '-dpng', '-r600', 'Figure S12');
disp('Saved: Figure S12');
end
%% ------------------ Helper function for residuals ------------------
function residuals = compute_residuals(params, w, R_data, Phi_data)
    % params = [δ_R, δ_G, ω₀]
    delta_R = params(1);
    delta_G = params(2);
    omega0 = params(3);
    
    % Avoid division by zero
    if omega0 == 0
        omega0 = 1e-10;
    end
    
    % Compute normalized frequency
    Omega = w / omega0;
    
    % Compute K_fit(Ω) = √((δ_R + iΩ)(δ_G + iΩ))
    K_fit = sqrt((delta_R + 1i*Omega) .* (delta_G + 1i*Omega));
    R_fit = real(K_fit);
    Phi_fit = imag(K_fit);
    
    % Compute residuals for both real and imaginary parts
    % We weight them equally
    residuals = [R_fit - R_data; Phi_fit - Phi_data];
    
    % Optional: Add regularization to ensure parameters are positive
    % (lsqnonlin already handles bounds, but this helps with stability)
    if delta_R < 0 || delta_G < 0 || omega0 <= 0
        residuals = [residuals; 1e6 * max(0, -delta_R); 1e6 * max(0, -delta_G); 1e6 * max(0, -omega0)];
    end
end
%% ==============================================================
function S = S_TCMT_2x2(omega, par)
D1 = (omega - par.w01) + 1i*(par.g1 + par.gc1)/2;
D2 = (omega - par.w02) + 1i*(par.g2 + par.gc2)/2;
den = (D1.*D2 - par.kappa^2);
M = [ par.gc1.*D2, sqrt(par.gc1*par.gc2).*par.kappa; ...
      sqrt(par.gc1*par.gc2).*par.kappa, par.gc2.*D1 ];
S = eye(2) - 1i*(M./den);
end

%% ==============================================================
function [z, u, v, D, GS, GL, pOut] = solve_z_onefreq(r1, r2, t, tauProd, pInit, opts)
% z = exp(-u + i v), u>=0, v real (unbounded; continuity via warm start)

tFloor = 1e-14;
if abs(t) < tFloor
    t = t + tFloor;
end

obj = @(p) objective_p(p, r1, r2, t, tauProd);

pOut = fminsearch(obj, pInit, opts);

[u, v] = decode_p(pOut);
z = exp(-u + 1i*v);

sqrtz = exp(-u/2 + 1i*v/2);
D = (sqrtz * tauProd) / t;

[GS, GL] = solve_G_from_linear(r1, r2, D, z);

end

function f = objective_p(p, r1, r2, t, tauProd)
[u, v] = decode_p(p);
z = exp(-u + 1i*v);
sqrtz = exp(-u/2 + 1i*v/2);
D = (sqrtz * tauProd) / t;

[GS, GL] = solve_G_from_linear(r1, r2, D, z);

res = D - (1 + GS*GL*z);

f = real(res*conj(res));

% gentle regularization to avoid z too close to ±1 where 1-z^2 small
f = f + 1e-8*abs(1 - z)^2 + 1e-8*abs(1 + z)^2;
end

function [u, v] = decode_p(p)
% u>=0 enforced by u=exp(p1); v free (NO WRAP!)
u = exp(p(1));
v = p(2);
end

function [GS, GL] = solve_G_from_linear(r1, r2, D, z)
den = (1 - z^2);
if abs(den) < 1e-14
    den = den + 1e-14;
end
GS = D*(r1 - z*r2)/den;
GL = D*(r2 - z*r1)/den;
end

%% ==============================================================
function S = S_WG_from_uv(u, v, D, GS, GL, tauProd)
% Consistent reconstruction using sqrtz = exp(-u/2 + i v/2)

z = exp(-u + 1i*v);
sqrtz = exp(-u/2 + 1i*v/2);

r1 = (GS + GL*z)/D;
r2 = (GL + GS*z)/D;
t  = (sqrtz * tauProd)/D;

S = [r1, t; t, r2];
end

function panel_letter_local(ax, letter)
text(ax, -0.12, 0.98, letter, 'Units','normalized', 'FontSize',14, 'FontWeight','bold', ...
    'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');
end
