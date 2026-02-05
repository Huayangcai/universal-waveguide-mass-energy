function FigureS1_nature()
%% ============================================================
% Figure S1:Dispersion classification from the dimensionless RLGC electrical length
% - Draw crossover lines for ALL cases (including Ωc=0) on log-x
% - Crossover line color matches the corresponding case curve (ROBUST: read from line handles)
% - Add a thick red dashed reference line at Ω = 1
%% ============================================================
clear; close all; clc;

%% ---------- Frequency grid ----------
Omega = logspace(-4, 3, 1400);
Omega_min = min(Omega);
Omega_max = max(Omega);

%% ---------- 3 representative cases ----------
cases = [ ...
    1e-2, 1e-2; ...
    1e-2, 5e-1; ...
    5e-1, 0    ...
];

%% ---------- auto labels from cases ----------
labels = cell(size(cases,1),1);
for k = 1:size(cases,1)
    dR = cases(k,1); dG = cases(k,2);
    if dR==0 && dG==0
        tag = 'Lossless';
    elseif dR==dG
        tag = 'Symmetric';
    elseif dG==0
        tag = 'Series-only';
    elseif dR==0
        tag = 'Shunt-only';
    else
        tag = 'Asymmetric';
    end
    labels{k} = sprintf('%s: $\\delta_R=%s,\\ \\delta_G=%s$', tag, sciTex(dR), sciTex(dG));
end

Nc   = size(cases,1);
cols = lines(Nc);
lw   = 2.4;

%% ---------- Compute Γ and Ωc ----------
GR_all = zeros(Nc, numel(Omega));
GI_all = zeros(Nc, numel(Omega));
Oc_all = zeros(Nc,1);

for k = 1:Nc
    dR = cases(k,1); dG = cases(k,2);
    [GR, GI] = localGammaSplit(Omega, dR, dG);
    GR_all(k,:) = GR;
    GI_all(k,:) = GI;
    Oc_all(k) = sqrt(max(dR*dG, 0));
end

%% ---------- Plot (2 rows, 1 col) ----------
% fig = figure('Color','w','Position',[90 60 1120 860]);
% % fig = figure('Color','w','Units','centimeters','Position',[2 2 8 8]);
% tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

fig = figure('Color','w','Units','centimeters','Position',[2 2 20 16]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
% ===================== Top: alphaL =====================
ax1 = nexttile(1); hold(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log','YScale','log','Layer','top');  % Layer top so lines are visible
grid(ax1,'on'); ax1.GridAlpha=0.2; ax1.MinorGridAlpha=0.1;
ax1.XMinorGrid='on'; ax1.YMinorGrid='on';

h1 = gobjects(Nc,1);
for k = 1:Nc
    h1(k) = plot(ax1, Omega, GR_all(k,:), 'LineWidth', lw, 'Color', cols(k,:));
end

% draw crossover lines AFTER curves (so they sit on top)
drawCrossoverLines(ax1, Oc_all, h1, Omega_min, Omega_max);

% ===== Add reference line at Ω = 1 (thick red dashed) =====
xline(ax1, 1, '--', 'Color', [1 0 0], 'LineWidth', 3.0, 'HandleVisibility','off');

% xlabel(ax1,'$\Omega$','Interpreter','latex');
ylabel(ax1,'$\alpha L=\Re[K(\Omega)]$','Interpreter','latex');
% title(ax1,'Real part: $\alpha L$','Interpreter','latex');
xlim(ax1,[Omega_min Omega_max]);

% ---- optional annotations (NOTE: these are in axis data units) ----
text(ax1, 0.08, 0.1, '$\Omega=\Omega_c$', 'Interpreter','latex');
text(ax1, 0.002,0.1, '$\Omega=\Omega_c$', 'Interpreter','latex');
text(ax1, 1,    0.1, '$\Omega=1$', 'Interpreter','latex');
% text(ax1,0.00015,0.2,'(a)')
% ===================== Legend =====================
% Use the TOP curve handles for legend so labels match cases
lg = legend(ax1, h1, labels, 'Interpreter','latex', ...
    'Location','east','Orientation','vertical');
lg.Box = 'on'; lg.FontSize = 11;

% ===================== Bottom: betaL =====================
ax2 = nexttile(2); hold(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log','YScale','log','Layer','top');
grid(ax2,'on'); ax2.GridAlpha=0.2; ax2.MinorGridAlpha=0.1;
ax2.XMinorGrid='on'; ax2.YMinorGrid='on';

h2 = gobjects(Nc,1);
for k = 1:Nc
    h2(k) = plot(ax2, Omega, GI_all(k,:), 'LineWidth', lw, 'Color', cols(k,:));
end

drawCrossoverLines(ax2, Oc_all, h2, Omega_min, Omega_max);

% ===== Add reference line at Ω = 1 (thick red dashed) =====
xline(ax2, 1, '--', 'Color', [1 0 0], 'LineWidth', 3.0, 'HandleVisibility','off');

% optional reference: lossless beta = Omega
plot(ax2, Omega, Omega, 'k:', 'LineWidth', 1.1, 'HandleVisibility','off');

xlabel(ax2,'$\Omega$','Interpreter','latex');
ylabel(ax2,'$\beta L=\Im[K(\Omega)]$','Interpreter','latex');
% title(ax2,'Imag part: $\beta L$','Interpreter','latex');
xlim(ax2,[Omega_min Omega_max]);
% text(ax2,0.00015,4000,'(b)')
   text(ax1, -0.08, 0.98, 'a', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');
   text(ax2, -0.08, 0.98, 'b', 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color','k');
%% -------------------- Save figure --------------------
% Set ALL font sizes right before output
set(findall(gcf,'-property','FontSize'),'FontSize',14);

filename = 'Figure S1';
print(fig,'-dpng','-r600',[filename,'.png']);
fprintf('Saved: %s.png\n', filename);

end
%% ===================== Local functions =====================
function drawCrossoverLines(ax, Oc_all, hcurves, Omega_min, Omega_max)
% Draw Ωc lines on a log-x axis; color taken from curve handles, then sanitized to valid RGB.

    for kk = 1:numel(Oc_all)
        Oc = Oc_all(kk);

        % Clamp Ωc for log-x visibility
        if Oc <= Omega_min
            Oc_plot = Omega_min * 1.02;
        elseif Oc >= Omega_max
            Oc_plot = Omega_max / 1.02;
        else
            Oc_plot = Oc;
        end

        % --- Robust color read + sanitize ---
        c = get(hcurves(kk), 'Color');     % may be 1x3, sometimes not ideal
        c = sanitizeRGB(c);               % ensure 1x3 double in [0,1]

        cl = xline(ax, Oc_plot, '--', 'LineWidth', 1.6, ...
            'Color', c, 'HandleVisibility','off');
        try, cl.Layer = 'top'; catch, end
    end
end

function c = sanitizeRGB(c)
% Ensure c is a valid RGB triplet in [0,1] for ConstantLine.Color

    % Convert strings like 'b'/'r' if ever occurs
    if ischar(c) || isstring(c)
        switch char(c)
            case {'y'}, c = [1 1 0];
            case {'m'}, c = [1 0 1];
            case {'c'}, c = [0 1 1];
            case {'r'}, c = [1 0 0];
            case {'g'}, c = [0 1 0];
            case {'b'}, c = [0 0 1];
            case {'w'}, c = [1 1 1];
            otherwise,  c = [0 0 0];
        end
        return;
    end

    % Numeric: make it 1x3 real double
    c = double(c);
    c = real(c(:).');         % row
    if numel(c) >= 3
        c = c(1:3);           % drop alpha etc.
    else
        c = [0 0 0];
        return;
    end

    % Fix NaN/Inf
    if any(~isfinite(c))
        c = [0 0 0];
        return;
    end

    % If looks like 0..255 scale, convert
    if max(c) > 1
        if max(c) <= 255
            c = c / 255;
        else
            c = c / max(c);   % fallback normalize
        end
    end

    % Clamp to [0,1]
    c = max(min(c, 1), 0);
end



function [Gamma_real, Gamma_imag] = localGammaSplit(Omega, dR, dG)
% Passive-branch split for Γ = sqrt[(dR+iΩ)(dG+iΩ)]
    P = dR.*dG - Omega.^2;
    Q = Omega.*(dR + dG);
    absZ = sqrt(P.^2 + Q.^2);
    Gamma_real = sqrt(max((absZ + P)/2, 0));
    Gamma_imag = sqrt(max((absZ - P)/2, 0));
    Gamma_imag = sign(Q + 0).*Gamma_imag;
end

function s = sciTex(x)
% LaTeX scientific notation formatter: 10^{e} or m×10^{e}
    if x==0, s='0'; return; end
    e = floor(log10(abs(x)));
    m = x/10^e;
    if abs(m-1) < 1e-12
        s = sprintf('10^{%d}', e);
    else
        s = sprintf('%.2g\\times10^{%d}', m, e);
    end
end


