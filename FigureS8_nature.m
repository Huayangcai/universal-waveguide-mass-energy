function FigureS8_nature()
%% ============================================================
% Figure S15 (symbol-consistent with SI):
% Four fundamental laws for waveguide power absorption & emission
%
% SI mappings:
%   K(Ω)=R+iΦ,  E(Ω)=exp(-2K)                 (S2-3)
%   Γ̄_g(Ω) = (Γ_L E + Γ_S)/(1+Γ_S Γ_L E)     (Möbius composition)
%   Ab = I - S^† S ,  Em = I - S S^†          (S6)
%   α_p = ε_p = 1 - σ_p^2                     (Law 1)
%   Law 2/4 require correct in/out SVD bases: input~V, output~U.
%% ============================================================

clear; close all; clc;
rng(7);

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',15);
set(groot,'DefaultTextFontSize',15);
set(groot,'DefaultLineLineWidth',2.6);

%% ---------------- Frequency grid ----------------
Omega = linspace(0,2,801);      % normalized frequency Ω
nW = numel(Omega);

Omega_star = 1.05;
[~,idx_star] = min(abs(Omega - Omega_star));

%% ---------------- Multi-eigenchannel construction ----------------
Nch = 6; % number of eigenchannels (SI: intrinsic channels p)

% For each channel p, build Γ̄_{g,p}(Ω)
GammaBar_all = zeros(Nch,nW);
designs = make_multiport_designs(Omega, Nch);
for p = 1:Nch
    [GammaBar_g,~,~,~] = compute_all_curves_symbolic(Omega, designs(p));
    GammaBar_all(p,:) = GammaBar_g;
end

% Construct a reciprocal (complex-symmetric) S(Ω).
% To make Law-1 interpretation clean, we set singular values σ_p(Ω)=|Γ̄_{g,p}(Ω)|.
% A Takagi-like synthesis: S = Q * diag(σ) * Q^T  (complex symmetric)
Q = rand_unitary(Nch);
S_all = zeros(Nch,Nch,nW);
for k = 1:nW
    sigma_k = abs(GammaBar_all(:,k));
    sigma_k = min(sigma_k, 0.999999);   % passive safety |σ|<1
    S_all(:,:,k) = Q * diag(sigma_k) * (Q.'); % reciprocal: S = S^T
end

%% ---------------- Pick reference frequency Ω* ----------------
S0  = S_all(:,:,idx_star);

% SI operators:
Ab  = eye(Nch) - (S0')*S0;   % Ab = I - S^† S  (input-side absorption operator)
Em  = eye(Nch) - S0*(S0');   % Em = I - S S^†  (output-side emission operator)

% SVD: S0 = U * Sig * V^†
[U,Sig,V] = svd(S0,'econ');
sigma = diag(Sig);
alpha_eig = 1 - sigma.^2;    % Law 1: α_p = ε_p = 1 - σ_p^2

%% ---------------- Monte-Carlo tests for Law 2 & Law 4 ----------------
M = 350;
alpha_list      = zeros(M,1);
eps_list_phase  = zeros(M,1);
eps_list_pc     = zeros(M,1);

for m = 1:M
    % Random coefficients in eigenchannel space
    c = randn(Nch,1) + 1i*randn(Nch,1);
    c = c / norm(c);

    % -------- Input state |i> in H_in (expand on V) --------
    i_state = V * c;
    i_state = i_state / norm(i_state);
    alpha_list(m) = real(i_state' * Ab * i_state);  % α(i)=<i|Ab|i>

    % -------- Law 2: same eigenchannel power splitting (|d_p|^2=|c_p|^2) --------
    % Build output |o> in H_out (expand on U) with SAME magnitudes |c| but random phases
    theta = 2*pi*rand(Nch,1);
    d = abs(c) .* exp(1i*theta);
    o_phase = U * d;
    o_phase = o_phase / norm(o_phase);
    eps_list_phase(m) = real(o_phase' * Em * o_phase); % ε(o)=<o|Em|o>

    % -------- Law 4: reciprocal phase-conjugate pairing --------
    % Matched emission channel corresponds to phase-conjugated coefficient vector.
    o_pc = U * conj(c);
    o_pc = o_pc / norm(o_pc);
    eps_list_pc(m) = real(o_pc' * Em * o_pc);
end

%% ---------------- Law 3: basis-independent sum rule ----------------
Ktr = 80;
err_in = zeros(Ktr,1);
trAb = real(trace(Ab));

for t = 1:Ktr
    Bin = rand_unitary(Nch); % a random orthonormal input basis
    sumAlpha = 0;
    for n = 1:Nch
        en = Bin(:,n);
        sumAlpha = sumAlpha + real(en' * Ab * en);
    end
    err_in(t) = sumAlpha - trAb; % should be ~0
end

%% -------------------- Figure layout --------------------
width_cm = 22.5;
height_cm = 18;
fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

% (a) Law 1
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
bar(ax1, 1:Nch, alpha_eig, 0.72);
plot(ax1, 1:Nch, alpha_eig, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
xlabel(ax1, 'eigenchannel $p$');
ylabel(ax1, '$\alpha_p=\varepsilon_p=1-\sigma_p^2$');

% (b) Law 2
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, alpha_list, eps_list_phase, 'ko', 'MarkerSize',4);
plot(ax2, [0 1], [0 1], 'k--', 'LineWidth',1.2);
xlabel(ax2, '$\alpha(\,|i\rangle\,)=\langle i|\mathbf{A_b}|i\rangle$');
ylabel(ax2, '$\varepsilon(\,|o\rangle\,)=\langle o|\mathbf{E_m}|o\rangle$');

% (c) Law 3
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
histogram(ax3, err_in, 20);
xlabel(ax3, '$\sum_n\alpha(|e_n\rangle)-\mathbf{Tr}(\mathbf{A_b})$');
ylabel(ax3, 'count');

% (d) Law 4
ax4 = nexttile(tlo,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
plot(ax4, alpha_list, eps_list_pc, 'ko', 'MarkerSize',4);
plot(ax4, [0 1], [0 1], 'k--', 'LineWidth',1.2);
xlabel(ax4, '$\alpha(\,|i\rangle\,)$');
ylabel(ax4, '$\varepsilon(\,|o_{\mathrm{pc}}\rangle\,)$');

% panel labels
add_panel_label(ax1,'a'); add_panel_label(ax2,'b');
add_panel_label(ax3,'c'); add_panel_label(ax4,'d');

%% ---------------- Diagnostics ----------------
fprintf('=== Diagnostics at Omega*=%.4f ===\n', Omega(idx_star));
fprintf('Tr(Ab)=%.6f, Tr(Em)=%.6f, |diff|=%.3e\n', real(trace(Ab)), real(trace(Em)), abs(real(trace(Ab))-real(trace(Em))));
fprintf('Law2 RMS(ε-α)=%.3e\n', rms(eps_list_phase-alpha_list));
fprintf('Law4 RMS(ε_pc-α)=%.3e\n', rms(eps_list_pc-alpha_list));
fprintf('Law3 RMS sum error=%.3e\n', rms(err_in));

filename = 'Figure S8';
print(gcf, '-dpng', '-r600', [filename, '.png']);
disp(['Figure saved as ', filename, '.png']);
end

%% ============================================================
% Symbol-consistent 1D waveguide reduction:
% Γ̄_g(Ω) from Möbius composition with K(Ω)=R+iΦ and E(Ω)=exp(-2K).
%% ============================================================
function [GammaBar_g, kappa, Theta, alpha_1d] = compute_all_curves_symbolic(Omega, d)

% Terminals: Γ_S(Ω), Γ_L(Ω)
if isscalar(d.rS), rS = d.rS*ones(size(Omega)); else, rS = d.rS; end
if isscalar(d.rL), rL = d.rL*ones(size(Omega)); else, rL = d.rL; end
phiS = d.phiS; phiL = d.phiL;

GammaS = rS .* exp(1i*phiS);
GammaL = rL .* exp(1i*phiL);

% Propagation: K(Ω)=R+iΦ,  E(Ω)=exp(-2K)
R   = d.R0 + d.R2*(Omega.^2);
Phi = d.Phi1 .* Omega;
K   = R + 1i*Phi;
Ert = exp(-2*K);  % round-trip factor E(Ω)

% Möbius (bilinear) boundary composition: Γ̄_g(Ω)
GammaBar_g = (GammaL.*Ert + GammaS) ./ (1 + GammaS.*GammaL.*Ert);

% Passive safety: enforce |Γ̄_g|<1
mag = abs(GammaBar_g);
mag = min(mag, 0.999999);
GammaBar_g = mag .* exp(1i*angle(GammaBar_g));

% Feedback descriptors: κ(Ω), Θ(Ω) from ΓS ΓL E
kappa = abs(GammaS.*GammaL) .* exp(-2*R);
Theta = wrapToPi(angle(GammaS.*GammaL) - 2*Phi);

% 1D absorptivity/emissivity proxy for the effective single channel
alpha_1d = 1 - abs(GammaBar_g).^2;
alpha_1d = max(0, min(1, real(alpha_1d)));
end

%% ---------------- Designs ----------------
function ds = make_multiport_designs(Omega, Nch)
ds = repmat(struct(), Nch, 1);
for p = 1:Nch
    ds(p).rS   = 0.15 + 0.35*rand;
    ds(p).rL   = 0.25 + 0.55*rand;
    ds(p).R0   = 0.04 + 0.14*rand;
    ds(p).R2   = 0.01 + 0.06*rand;
    ds(p).Phi1 = (0.9 + 0.7*rand)*pi;

    ds(p).phiS = 2*pi*rand + (0.4+0.4*rand)*(2*pi)*Omega;
    ds(p).phiL = 2*pi*rand + (0.4+0.4*rand)*(2*pi)*Omega;
end
end

%% ---------------- Helpers ----------------
function Q = rand_unitary(N)
Z = randn(N,N) + 1i*randn(N,N);
[Q,~] = qr(Z,0);
Q = Q / exp(1i*angle(det(Q))/N);
end

function add_panel_label(ax, ch)
% top-left in axes normalized coordinates
text(ax, -0.12, 0.98, ch, 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');
end
