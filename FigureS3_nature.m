function FigureS3_nature()
% ======================================================================
% Figure S3:Universal error scaling of the low-mass approximation
% Universal error scaling of the low-mass approximation
%
% Key physics checks:
%   - ln D^+ is real/physical only if |beta|<1
%   - For 2nd-order series, also require U2 > |S| (=> U2>0 and |beta2|<1)
%     Otherwise beta/gamma/lnD become non-physical and should be masked.
%
% Adds:
%   - "low-mass" label at x = x_lowmass
%   - Slope guides and labels for BOTH 1st and 2nd orders in low-mass regime:
%       1st:  U,beta,gamma ~ \chi^{-4},   |Δ ln D^+| ~ \chi^{-2}
%       2nd:  U,beta,gamma ~ \chi^{-6},   |Δ ln D^+| ~ \chi^{-4}
%   - Improved annotation placement (auto anchor from valid data)
% ======================================================================

clear; close all; clc;

%% ---------------- Global style ----------------
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultTextFontName','Helvetica');
set(0,'DefaultAxesFontSize',13);
set(0,'DefaultLineLineWidth',1.8);

%% ---------------- User controls ----------------
x_min = 1e-1;
x_max = 1e2;
N     = 2000;

show_second_order = true;
show_two_branches = false;

x_lowmass = 10;                 % low-mass shading + label

tol_beta = 1e-12;               % lnD physical cutoff: |beta| < 1 - tol_beta

%% ---------------- Normalized variable ----------------
x = logspace(log10(x_min), log10(x_max), N);
g = 1.0;

Spos =  g * x;
Sneg = -g * x;

clip   = @(z) max(z, 1e-14);
relpct = @(a,b) 100 * abs((a-b) ./ clip(b));

%% ---------------- Helpers ----------------
safe_lnD = @(beta) local_safe_atanh(beta, tol_beta);

% choose nearest valid y at x0_guess; if invalid, search within a band
pick_y_at = @(curve, xvec, x0_guess) local_pick_y(curve, xvec, x0_guess);

%% ---------------- Exact quantities ----------------
U_exact_pos = sqrt(Spos.^2 + g^2);
U_exact_neg = sqrt(Sneg.^2 + g^2);

beta_exact_pos  = Spos ./ U_exact_pos;
beta_exact_neg  = Sneg ./ U_exact_neg;

gamma_exact_pos = U_exact_pos / g;
gamma_exact_neg = U_exact_neg / g;

lnD_exact_pos = safe_lnD(beta_exact_pos);
lnD_exact_neg = safe_lnD(beta_exact_neg);

%% ---------------- 1st-order approximation ----------------
U1_pos = abs(Spos) + g^2 ./ (2*abs(Spos));
U1_neg = abs(Sneg) + g^2 ./ (2*abs(Sneg));

beta1_pos  = Spos ./ U1_pos;
beta1_neg  = Sneg ./ U1_neg;

gamma1_pos = U1_pos / g;
gamma1_neg = U1_neg / g;

lnD1_pos = safe_lnD(beta1_pos);
lnD1_neg = safe_lnD(beta1_neg);

%% ---------------- 2nd-order approximation ----------------
if show_second_order
    U2_pos = abs(Spos) + g^2./(2*abs(Spos)) - g^4./(8*abs(Spos).^3);
    U2_neg = abs(Sneg) + g^2./(2*abs(Sneg)) - g^4./(8*abs(Sneg).^3);

    beta2_pos  = Spos ./ U2_pos;
    beta2_neg  = Sneg ./ U2_neg;

    gamma2_pos = U2_pos / g;
    gamma2_neg = U2_neg / g;

    lnD2_pos = safe_lnD(beta2_pos);
    lnD2_neg = safe_lnD(beta2_neg);

    % ---- Physics-first validity for ALL 2nd-order quantities ----
    % require U2 >= |S|  (guarantees U2>0 and |beta2|<=1)
    valid2_pos = isfinite(U2_pos) & (U2_pos > abs(Spos)); % strict >
    valid2_neg = isfinite(U2_neg) & (U2_neg > abs(Sneg)); % same since abs(Sneg)=abs(Spos)

    % mask ALL 2nd-order curves in non-physical region
    U2_pos(~valid2_pos)       = NaN;   beta2_pos(~valid2_pos)    = NaN;
    gamma2_pos(~valid2_pos)   = NaN;   lnD2_pos(~valid2_pos)     = NaN;

    U2_neg(~valid2_neg)       = NaN;   beta2_neg(~valid2_neg)    = NaN;
    gamma2_neg(~valid2_neg)   = NaN;   lnD2_neg(~valid2_neg)     = NaN;
end

%% ---------------- Errors ----------------
% 1st
epsU1_pos   = relpct(U1_pos,     U_exact_pos);
epsB1_pos   = relpct(beta1_pos,  beta_exact_pos);
epsG1_pos   = relpct(gamma1_pos, gamma_exact_pos);
epsLnD1_pos = abs(lnD1_pos - lnD_exact_pos);
epsLnD1_pos(~isfinite(lnD1_pos) | ~isfinite(lnD_exact_pos)) = NaN;

% 2nd
if show_second_order
    epsU2_pos   = relpct(U2_pos,     U_exact_pos);
    epsB2_pos   = relpct(beta2_pos,  beta_exact_pos);
    epsG2_pos   = relpct(gamma2_pos, gamma_exact_pos);
    epsLnD2_pos = abs(lnD2_pos - lnD_exact_pos);
    epsLnD2_pos(~isfinite(lnD2_pos) | ~isfinite(lnD_exact_pos)) = NaN;
end

%% ---------------- Figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 12]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log','YScale','log');

% Shaded low-mass regime
yl_temp = [1e-8, 1e6];
patch_x = [x_lowmass, x_max, x_max, x_lowmass];
patch_y = [yl_temp(1), yl_temp(1), yl_temp(2), yl_temp(2)];
p = patch('XData',patch_x,'YData',patch_y, ...
    'FaceColor',[0.92 0.92 0.92],'EdgeColor','none', ...
    'FaceAlpha',0.6,'HandleVisibility','off');
uistack(p,'bottom');

% 1st order (solid)
loglog(x, epsU1_pos,   '-',  'LineWidth',2.4, 'DisplayName','$\varepsilon_\mathcal{U}$ (1st)');
loglog(x, epsB1_pos,   '-',  'LineWidth',2.4, 'DisplayName','$\varepsilon_{\beta_{w}}$ (1st)');
loglog(x, epsG1_pos,   '-',  'LineWidth',2.4, 'DisplayName','$\varepsilon_{\gamma_{w}}$ (1st)');
loglog(x, epsLnD1_pos, '-',  'LineWidth',2.4, 'DisplayName','$|\Delta \ln D^+|$ (1st)');

% 2nd order (dashed)
if show_second_order
    loglog(x, epsU2_pos,   '--', 'LineWidth',1.7, 'DisplayName','$\varepsilon_\mathcal{U}$ (2nd)');
    loglog(x, epsB2_pos,   '--', 'LineWidth',1.7, 'DisplayName','$\varepsilon_{\beta_{w}}$ (2nd)');
    loglog(x, epsG2_pos,   '--', 'LineWidth',1.7, 'DisplayName','$\varepsilon_{\gamma_{w}}$ (2nd)');
    loglog(x, epsLnD2_pos, '--', 'LineWidth',1.7, 'DisplayName','$|\Delta \ln D^+|$ (2nd)');
end

xlabel('$\chi = |\mathcal{S}|/|\Gamma_g|$');
ylabel('Error ( \% for $\mathcal{U},\beta_{w},\gamma_{w}$ ; absolute for $|\Delta\ln D^+|$ )');
% title('Universal error scaling of the low-mass approximation');
xlim([x_min, x_max]);

% Robust ylim
allY = [epsU1_pos(:); epsB1_pos(:); epsG1_pos(:); epsLnD1_pos(:)];
if show_second_order
    allY = [allY; epsU2_pos(:); epsB2_pos(:); epsG2_pos(:); epsLnD2_pos(:)];
end
allY = allY(isfinite(allY) & allY>0);
ymin = max(min(allY)*0.6, 1e-8);
ymax = min(max(allY)*1.8, 1e6);
ylim([ymin, ymax]);

% Update shading to current ylim
yl = ylim(ax);
set(p,'YData',[yl(1), yl(1), yl(2), yl(2)]);

%% ---------------- low-mass line + label ----------------
xline(x_lowmass,'k--','LineWidth',1.2,'HandleVisibility','off');

y_label = 10^(0.92*log10(yl(2)) + 0.08*log10(yl(1)));
text(x_lowmass*1.03, y_label-95, '\textbf{low-mass}', ...
    'Interpreter','latex','FontSize',12,'Rotation',90, ...
    'HorizontalAlignment','left','VerticalAlignment','top');

%% ---------------- Slope guides & labels (1st + 2nd) ----------------
% Pick anchors inside low-mass region where curves are finite
x0_1 = max(8, 0.8*x_lowmass);          % for 1st guides
x0_2 = max(20, 2.0*x_lowmass);         % for 2nd guides (deeper in low-mass)

y0_U1  = pick_y_at(epsU1_pos,   x, x0_1);
y0_LD1 = pick_y_at(epsLnD1_pos, x, x0_1);

% guide span (keep inside x_max)
xref1 = [x0_1, min(x0_1*8, x_max)];
if isfinite(y0_U1) && isfinite(y0_LD1)
    % 1st: U/beta/gamma ~ \chi^-4
    yref_1_u = y0_U1 * (xref1/x0_1).^(-4);
    loglog(xref1, yref_1_u, 'k-.', 'LineWidth',1.2,'HandleVisibility','off');

    % 1st: lnD ~ \chi^-2
    yref_1_d = y0_LD1 * (xref1/x0_1).^(-2);
    loglog(xref1, yref_1_d, 'k:', 'LineWidth',1.2,'HandleVisibility','off');

    % labels (slight offsets)
    text(xref1(end), 0.0005, '$\propto \chi^{-4}$ (1st)', ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','right');
    text(xref1(end)+10, yref_1_d(end)*0.80, '$\propto \chi^{-2}$ (1st)', ...
        'Interpreter','latex','FontSize',11,'HorizontalAlignment','right');
end

if show_second_order
    y0_U2  = pick_y_at(epsU2_pos,   x, x0_2);
    y0_LD2 = pick_y_at(epsLnD2_pos, x, x0_2);

    xref2 = [x0_2, min(x0_2*4, x_max)];
    if isfinite(y0_U2) && isfinite(y0_LD2)
        % 2nd: U/beta/gamma ~ \chi^-6
        yref_2_u = y0_U2 * (xref2/x0_2).^(-6);
        loglog(xref2, yref_2_u, 'k--', 'LineWidth',1.1,'HandleVisibility','off');

        % 2nd: lnD ~ \chi^-4
        yref_2_d = y0_LD2 * (xref2/x0_2).^(-4);
        loglog(xref2, yref_2_d, 'k-', 'LineWidth',0.9,'HandleVisibility','off');

        % labels (offsets to avoid overlaps)
        text(xref2(end), 1e-7, '$\propto \chi^{-6}$ (2nd)', ...
            'Interpreter','latex','FontSize',11,'HorizontalAlignment','right');
        text(30, 3e-8, '$\propto \chi^{-4}$ (2nd)', ...
            'Interpreter','latex','FontSize',11,'HorizontalAlignment','right');
    end
end

legend('Location','southwest','FontSize',10,'Box','on');

%% ---------------- Export ----------------
filename = 'Figure S3';
print(fig,'-dpng','-r600',[filename,'.png']);
fprintf('Saved: %s.png\n', filename);

end

% ===== local helpers =====
function lnD = local_safe_atanh(beta, tol_beta)
    lnD = NaN(size(beta));
    valid = isfinite(beta) & (abs(beta) < 1 - tol_beta);
    lnD(valid) = atanh(beta(valid));
end

function y0 = local_pick_y(curve, xvec, x0_guess)
    % pick nearest finite y around x0_guess; if not found, search in neighborhood
    y0 = NaN;
    finiteMask = isfinite(curve) & (curve > 0) & isfinite(xvec) & (xvec > 0);
    if ~any(finiteMask), return; end

    % nearest on log scale (more appropriate for log axes)
    [~, idx0] = min(abs(log10(xvec) - log10(x0_guess)));
    if finiteMask(idx0)
        y0 = curve(idx0);
        return;
    end

    % otherwise search outward
    idxList = find(finiteMask);
    [~, k] = min(abs(log10(xvec(idxList)) - log10(x0_guess)));
    y0 = curve(idxList(k));
end
