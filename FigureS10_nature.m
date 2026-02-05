function FigureS10_nature()
% ===============================================================
% Classic resonators (2x2) + selected waveguide laws (v2)
% Fixes:
%   - Law 1: consistent sorting of alpha_M and eig(Ab)
%   - Law 2: equal-weight pairing using SVD (U,V) like tunneling script
% Visual:
%   - emphasize corridor Lambda only; lighter kappa
%   - remove minor grid; cleaner labels; info textbox w/ white background
%   - inset color matches resonator; inset auto-placed to avoid peak occlusion
% Scaling:
%   - (a,b,d) Lambda on log scale
%   - (c) Helmholtz Lambda on linear scale (tight ylim)
% Export: FigureXX_ClassicResonators_Laws_Main_v2.png (600 dpi)
% ===============================================================

clear; close all; clc;

%% ---------------- Global style ----------------
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',13);
set(groot,'DefaultTextFontSize',13);
set(groot,'DefaultLineLineWidth',2.2);

%% ---------------- Four classic resonators ----------------
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

%% ---------------- Dispersion parameters (demo) ----------------
delta_R = 0.20;
delta_G = 0.10;

%% ---------------- Frequency grid ----------------
Omega_min = 0.05;
Omega_max = 10.0;
N_Omega   = 2200;
Omega = linspace(Omega_min, Omega_max, N_Omega);

[R_Om, Phi_Om] = Ksplit(Omega, delta_R, delta_G);
K  = R_Om + 1i*Phi_Om;
a1 = exp(-K);        % one-way amplitude
a2 = exp(-2*K);      % round-trip amplitude

wrapPi = @(x) mod(x + pi, 2*pi) - pi;

%% ---------------- Resonance corridor ----------------
ThetaSigma = pi/12;
minSep_Om  = 0.40;
maxMarkers = 1;      % one main resonance per resonator
eps_mag    = 1e-12;

%% ---------------- Precompute per resonator ----------------
RES = struct([]);
for i = 1:nCfg
    name   = cfg{i,1};
    GamS0  = cfg{i,2};
    GamL0  = cfg{i,3};

    % --- F(Ω), D(Ω) convention (like Figure24)
    F = GamS0 .* GamL0 .* a2;        % F = ΓS ΓL e^{-2K}
    D = 1 + F;

    kap   = abs(F);
    theta = angle(F);
    det   = abs(wrapPi(theta - pi));
    Lam   = 1 ./ max(abs(D), eps_mag);      % Λ = 1/|1+F|

    % --- pick one resonance within corridor: max Λ in corridor
    [Ores, ires] = pick_resonances(Omega, det, Lam, maxMarkers, minSep_Om, 25, ThetaSigma);
    if isempty(ires)
        [~, idx] = max(Lam);
        ires = idx; Ores = Omega(idx);
    end
    idx_res = ires(1);
    Om_res  = Omega(idx_res);

    % --- Build passive reciprocal 2-port cavity S(Ω) consistent with D=1+F
    % FP: den = 1 - r1*r2*a2. Want den = 1 + ΓS ΓL a2 => set r1=ΓS, r2=-ΓL.
    r1 = clamp_reflection(GamS0);
    r2 = -clamp_reflection(GamL0);

    t1 = unitary_t_from_r(r1);
    t2 = unitary_t_from_r(r2);

    den = 1 - r1*r2.*a2;     % equals 1 + ΓS ΓL a2
    S21 = (t1*t2) .* a1 ./ den;
    S12 = S21;
    S11 = r1 + (t1^2) .* r2 .* a2 ./ den;
    S22 = r2 + (t2^2) .* r1 .* a2 ./ den;

    % Evaluate S at resonance
    Sres = [S11(idx_res), S12(idx_res); S21(idx_res), S22(idx_res)];
    Ab = eye(2) - (Sres')*Sres;
    Em = eye(2) - Sres*(Sres');
    Ab = (Ab + Ab')/2;
    Em = (Em + Em')/2;

    [U, Ssvd, V] = svd(Sres, 'econ');
    sig = diag(Ssvd);
    alpha_M = max(min(1 - sig.^2, 1), 0); % Law 1 RHS

    TrAb = real(trace(Ab));
    TrEm = real(trace(Em));

    % Passivity check
    eSS = eig((Sres')*Sres);
    maxEigSHS = max(real(eSS));

    % --- Samples for Law 2 (equal-weight pairing) and Law 4 (reciprocity pairing)
    nPairs = 180;  % fewer points -> cleaner inset
    [alpha_i, eps_o] = law2_equal_weight_samples(Ab, Em, U, V, nPairs);
    [alpha_t, eps_t] = law4_reciprocity_samples(Ab, Em, nPairs);

    RES(i).name = name;
    RES(i).GamS = GamS0;
    RES(i).GamL = GamL0;
    RES(i).F    = F;
    RES(i).D    = D;
    RES(i).Lam  = Lam;
    RES(i).kap  = kap;
    RES(i).det  = det;
    RES(i).idx_res = idx_res;
    RES(i).Om_res  = Om_res;
    RES(i).Sres = Sres;
    RES(i).Ab = Ab;
    RES(i).Em = Em;
    RES(i).U  = U;
    RES(i).V  = V;
    RES(i).alpha_M = alpha_M;
    RES(i).TrAb = TrAb;
    RES(i).TrEm = TrEm;
    RES(i).maxEigSHS = maxEigSHS;
    RES(i).alpha_i = alpha_i;
    RES(i).eps_o   = eps_o;
    RES(i).alpha_t = alpha_t;
    RES(i).eps_t   = eps_t;
end

%% ---------------- Figure layout ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 18]);
tlo = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for i = 1:nCfg
    ax = nexttile(tlo, i);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on'); % no minor grid

    Lam  = RES(i).Lam;
    kap  = RES(i).kap;
    det  = RES(i).det;
    Omr  = RES(i).Om_res;
    idxr = RES(i).idx_res;

    % Corridor mask
    maskC = det < ThetaSigma;

    % ------ Left axis: Lambda
    yyaxis(ax,'left');
    ax.YColor = [0 0 0];

    LamC = Lam;
    LamC(~maskC) = NaN; % show corridor only
    LamC = max(LamC, 1e-9);

    plot(ax, Omega, LamC, '-', 'Color', colors(i,:), 'LineWidth', 2.9);

    if i == 3
        % Helmholtz: linear Lambda, tighten y-range for visibility
        set(ax,'YScale','linear');
        ylo = max(1, min(Lam(maskC))*0.995);
        yhi = max(Lam(maskC))*1.02;
        if ~isfinite(ylo) || ~isfinite(yhi) || yhi <= ylo
            ylo = 1; yhi = max(Lam)*1.05;
        end
        ylim(ax,[ylo yhi]);
        ylabel(ax,'$\Lambda(\Omega)=1/|1+F|$');
    else
        % a,b,d: log Lambda
        set(ax,'YScale','log');
        ylim(ax,[1 max(Lam(maskC))*1.35 + 1e-9]);
        ylabel(ax,'$\Lambda(\Omega)=1/|1+F|$ (log)');
    end

    xline(ax, Omr, 'k:', 'LineWidth', 1.05, 'HandleVisibility','off');

    % ------ Right axis: kappa (lighter)
    yyaxis(ax,'right');
    ax.YColor = [0.25 0.25 0.25];
    plot(ax, Omega, kap, '-', 'Color', lighten_color(colors(i,:),0.72), 'LineWidth', 3.5);
    ylim(ax,[0 1]);
    yticks(ax,0:0.2:1);

    % Only right column keeps right-axis ylabel (reduce clutter)
    if ismember(i,[2 4])
        ylabel(ax,'$\kappa(\Omega)=|F(\Omega)|$');
    else
        ylabel(ax,'');
    end

    % ------ Common x
    xlabel(ax,'$\Omega$');
    xlim(ax,[Omega_min Omega_max]);

    % ------ Title
    yyaxis(ax,'left');
    title(ax, sprintf('%s', RES(i).name), 'FontWeight','bold');

    % ------ Info box (avoid overlapping axes/data)
    infoStr = sprintf(['$\\Omega_{\\rm res}=%.2f$\\quad $\\kappa=%.2f$\\quad ' ...
                       '$\\max\\ (\\mathbf{S^\\dagger S})=%.3f$'], ...
                       Omr, kap(idxr), RES(i).maxEigSHS);
    text(ax, 0.02, 0.06, infoStr, 'Units','normalized', 'Interpreter','latex', ...
        'FontSize', 12, 'Color', [0 0 0], ...
        'BackgroundColor','w', 'Margin',2, 'EdgeColor','none');

    % ------ Panel label
    panelChar = char('a' + (i-1));
    text(ax, -0.13, 0.98, panelChar, 'Units','normalized', ...
        'FontSize',16,'FontWeight','bold','FontName','Helvetica', ...
        'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','k');

    % ------ Law inset
    add_law_inset(fig, ax, RES(i), i, colors(i,:), Omega_min, Omega_max);

    hold(ax,'off');
end

% Export
% set(findall(gcf,'-property','FontSize'),'FontSize',13);
% print(fig, '-dpng', '-r600', 'FigureXX_ClassicResonators_Laws_Main_v2.png');
% fprintf('Saved: FigureXX_ClassicResonators_Laws_Main_v2.png\n');
    filename = 'FigureS10';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);

end

%% ============================================================
% Inset: selected law diagnostics (auto placement)
function add_law_inset(fig, axParent, R, idxPanel, cMain, Omega_min, Omega_max)
pos = get(axParent,'Position');

% auto place inset to avoid resonance line: if res on right half -> inset left
midOmega = 0.5*(Omega_min + Omega_max);
placeLeft = (R.Om_res > midOmega);

w = pos(3)*0.34;
h = pos(4)*0.34;
if placeLeft
    x = pos(1) + pos(3)*0.08;
else
    x = pos(1) + pos(3)*0.58;
end
y = pos(2) + pos(4)*0.60;

axI = axes('Parent',fig,'Position',[x y w h]); %#ok<LAXES>
hold(axI,'on'); box(axI,'on');
set(axI,'FontSize',9);

switch idxPanel
    case 1
        % Quarter-wave: Law 1 (alpha_Mp vs eig(Ab)) - with consistent sorting
        alphaM = R.alpha_M(:);
        eigAb  = real(eig(R.Ab));

        alphaM = sort(alphaM,'descend');
        eigAb  = sort(eigAb,'descend');

        p = (1:2)';
        bar(axI, p, alphaM, 0.55, 'FaceAlpha',0.85, 'EdgeColor','none', 'FaceColor', cMain);
        plot(axI, p, eigAb, 'k.', 'MarkerSize', 12);

        xlim(axI,[0.5 2.5]); xticks(axI,[1 2]);
        ylim(axI,[0 1]);
        title(axI,'Law 1', 'FontWeight','bold');
        xlabel(axI,'$p$'); ylabel(axI,'$\alpha_{Mp}$');

        maxDiff = max(abs(alphaM - eigAb));
        % text(axI,0.05,0.88, sprintf('$\\max|\\Delta|=%.1e$', maxDiff), ...
        %     'Units','normalized','Interpreter','latex','FontSize',8);

    case 2
        % Fabry-Perot: Law 3 (trace invariance)
        TrA = R.TrAb; TrE = R.TrEm;
        vals = [TrA, TrE];
        bar(axI, [1 2], vals, 0.55, 'FaceAlpha',0.85, 'EdgeColor','none', 'FaceColor', cMain);
        xticks(axI,[1 2]); xticklabels(axI,{'Tr$(\mathbf{A_b})$','Tr$(\mathbf{E_m})$'});
        ylabel(axI,'Trace');
        title(axI,'Law 3', 'FontWeight','bold');
        ylim(axI,[0 max(vals)*1.25 + 1e-9]);

        % text(axI,0.05,0.88, sprintf('$|\\Delta|=%.2e$', abs(TrA-TrE)), ...
        %     'Units','normalized','Interpreter','latex','FontSize',8);

    case 3
        % Helmholtz: Law 2 scatter alpha_i vs epsilon_o (equal-weight pairing)
        a = max(R.alpha_i(:),0);
        e = max(R.eps_o(:),0);
        m = max([a; e; 1e-12]);

        scatter(axI, a, e, 14, 'filled', 'MarkerFaceAlpha',0.55, 'MarkerEdgeAlpha',0.0, ...
            'MarkerFaceColor', cMain);
        plot(axI, [0 m], [0 m], 'k--', 'LineWidth', 1.0);
        axis(axI,[0 m 0 m]); axis(axI,'square');
        title(axI,'Law 2', 'FontWeight','bold');
        xlabel(axI,'$\alpha_i$'); ylabel(axI,'$\epsilon_o$');

        rmsErr = sqrt(mean((a-e).^2));
        % text(axI,0.05,0.88, sprintf('RMS=%.2e', rmsErr), ...
        %     'Units','normalized','Interpreter','none','FontSize',8);

    case 4
        % Tunnelling: Law 4 scatter alpha(i) vs epsilon(i*)
        a = max(R.alpha_t(:),0);
        e = max(R.eps_t(:),0);
        m = max([a; e; 1e-12]);

        scatter(axI, a, e, 14, 'filled', 'MarkerFaceAlpha',0.55, 'MarkerEdgeAlpha',0.0, ...
            'MarkerFaceColor', cMain);
        plot(axI, [0 m], [0 m], 'k--', 'LineWidth', 1.0);
        axis(axI,[0 m 0 m]); axis(axI,'square');
        title(axI,'Law 4', 'FontWeight','bold');
        xlabel(axI,'$\alpha_i$'); ylabel(axI,'$\epsilon_{i^*}$');

        rmsErr = sqrt(mean((a-e).^2));
        % text(axI,0.05,0.88, sprintf('RMS=%.2e', rmsErr), ...
        %     'Units','normalized','Interpreter','none','FontSize',8);
end

hold(axI,'off');
end

%% ============================================================
% Ksplit: same as Figure24
function [R, Phi] = Ksplit(Omega, delta_R, delta_G)
P = delta_R*delta_G - Omega.^2;
Q = Omega.*(delta_R + delta_G);
Zabs = sqrt(P.^2 + Q.^2);
R   = sqrt((Zabs + P)/2);
Phi = sqrt((Zabs - P)/2);
R = real(R); Phi = real(Phi);
end

%% ============================================================
% pick resonances: local maxima of Lam within det corridor
function [Ores, ires] = pick_resonances(Omega, detAbs, Lam, maxNum, minSep, win, corridor)
Ores = []; ires = [];

cand = find_localmax(Lam);
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

function idx = find_localmax(x)
% Toolbox-free local maxima indices
x = x(:);
if numel(x) < 3
    idx = [];
    return;
end
mid = (x(2:end-1) > x(1:end-2)) & (x(2:end-1) >= x(3:end));
idx = find(mid) + 1;
end

%% ============================================================
% Unitary interface transmission for a given reflection r
function t = unitary_t_from_r(r)
% Lossless reciprocal interface: |r|^2 + |t|^2 = 1 and r t* + t r* = 0
% => t phase = arg(r) + pi/2
tmag = sqrt(max(0, 1 - abs(r)^2));
t = tmag * exp(1i*(angle(r) + pi/2));
end

%% ============================================================
% Clamp any reflection magnitude >= 1 to < 1 (passivity & avoid t=0 degeneracy)
function r2 = clamp_reflection(r)
mag = abs(r);
if mag >= 0.999
    r2 = (0.999) * r / max(mag, 1e-12);
else
    r2 = r;
end
end

%% ============================================================
% Law 2: equal-weight pairing via SVD (like Figure05 tunneling script)
function [alpha_i, eps_o] = law2_equal_weight_samples(Ab, Em, U, V, nPairs)
alpha_i = zeros(nPairs,1);
eps_o   = zeros(nPairs,1);
for k = 1:nPairs
    c = randn(2,1)+1i*randn(2,1); c = c/norm(c);
    mags = abs(c);
    phs  = 2*pi*rand(2,1);
    d    = mags .* exp(1i*phs);
    d    = d / norm(d);

    i_vec = V*c;
    o_vec = U*d;

    alpha_i(k) = real(i_vec' * Ab * i_vec);
    eps_o(k)   = real(o_vec' * Em * o_vec);
end
end

%% ============================================================
% Law 4: reciprocity pairing alpha(i) = epsilon(i*)
function [alpha_t, eps_t] = law4_reciprocity_samples(Ab, Em, nPairs)
alpha_t = zeros(nPairs,1);
eps_t   = zeros(nPairs,1);
for k = 1:nPairs
    iVec = randn(2,1)+1i*randn(2,1); iVec = iVec / norm(iVec);
    istar = conj(iVec);
    alpha_t(k) = real(iVec'  * Ab * iVec);
    eps_t(k)   = real(istar' * Em * istar);
end
end

%% ============================================================
function c2 = lighten_color(c, f)
c2 = (1-f)*c + f*[1 1 1];
c2 = max(0, min(1, c2));
end
