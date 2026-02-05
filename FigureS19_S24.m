clc; clear; close all;
warning('off','all');

%% ---------------- Constant Definitions ----------------
T = 298.0; F = 96485.33212; Rg = 8.314462618;
f = F/(Rg*T);
nernst = log(10)/f;   % 2.303RT/F ~ 0.0591 V at 298K

%% ---------------- Control Parameters ----------------
maxExp = 70;

% Stage 1 parameters (no ohmic drop fitting)
st1 = struct('maxIter', 900, 'maxFeval', 40000, 'tolFun', 1e-6, ...
             'tolX', 1e-6, 'nStarts', 4);

% Stage 2 parameters (full fitting)
st2 = struct('newtonMaxIter', 10, 'newtonTol', 5e-11, 'lsqMaxIter', 700, ...
             'lsqMaxFeval', 60000, 'tolFun', 1e-6, 'tolX', 1e-6, 'nStarts', 4);

% Model selection (AIC criterion)
sel = struct('enableAIC', false, 'AIC_deltaMin', 2.0, 'forceAB', true);

% Stage 2 weak regularization parameters
st2.wRohm = 0.01;
st2.Rohm0 = 0.020;
st2.RohmScale = 0.80;

% Penalty term to prevent j* from approaching 0
st2.wCollapse = 0.01;
st2.jstarFloor = 1e-10;

%% ---------------- Diagnostic Parameters ----------------
prune = struct();
prune.enable = true;
prune.fracTol  = 0.03;   % Area contribution <3% considered weak
prune.peakTol  = 0.05;
prune.nGrid    = 500;
prune.etaPctRange = [5 95];
prune.cathodicOnly = true;  % Only for cathodic process
prune.jRelMin  = 1e-3;

%% ---------------- Read Data ----------------
csvFile = 'NC 2025-SI FigS9b_sourcedata.csv';
if ~exist(csvFile,'file')
    here = fileparts(mfilename('fullpath'));
    alt  = fullfile(here, csvFile);
    if exist(alt,'file'), csvFile = alt; else, error('Cannot find Fig3b_sourcedata.csv'); end
end

TBL = readtable(csvFile);
v = lower(string(TBL.Properties.VariableNames));
col_pH = find(contains(v,'ph'),1,'first');
col_U  = find(contains(v,'u'),1,'first');
col_J  = find(contains(v,'current') | contains(v,'j'),1,'first');

if isempty(col_pH) || isempty(col_U) || isempty(col_J)
    error('CSV columns not recognized: %s', strjoin(string(TBL.Properties.VariableNames),', '));
end

pHcol = TBL{:,col_pH};
Uraw  = TBL{:,col_U};
j_raw = TBL{:,col_J};

mask = isfinite(Uraw) & isfinite(j_raw);
Uraw = Uraw(mask); j_raw = j_raw(mask); pHcol = pHcol(mask);

pHnum  = parse_pH_column(pHcol);
pHvals = sort(unique(pHnum(:)));
nG = numel(pHvals);

fprintf('Loaded %d points across %d pH conditions. (Assume U is vs SHE)\n', numel(Uraw), nG);
    
%% ---------------- Fit each pH curve ----------------
RES = repmat(struct('pH', [], 'P', [], 'P1', [], 'Iref', [], ...
                    'U', [], 'eta', [], 'j', [], 'modelTag', ""), nG, 1);
    
fprintf('\n=== Fitting progress (AIC selection: AB / A / B) ===\n');
    
for g = 1:nG
    ph = pHvals(g);
    idx = (pHnum == ph);
        
    Ug = Uraw(idx);
    jg = j_raw(idx);
        
    % HER overpotential vs SHE
    etag = Ug + nernst * ph;
        
    % Sort by η
    [etag, ord] = sort(etag(:));
    Ug = Ug(ord); jg = jg(ord);
        
    Iref = median(abs(jg(abs(jg) > 0)));
    if ~isfinite(Iref) || Iref <= 0, Iref = max(abs(jg)) + 1e-12; end
        
    fprintf('\n--- pH=%.3g (%d points) ---\n', ph, numel(etag));
        
    % ===================== Fit AB mode =====================
    fprintf('  Stage 1(AB): Initial fitting without ohmic drop...\n');
    P1_AB = fit_stage1_noohm(etag, jg, Iref, f, maxExp, st1);
    
    fprintf('  Stage 2(AB): Full implicit optimization...\n');
    P2_AB = fit_stage2_full(etag, jg, Iref, f, maxExp, st2, P1_AB);
    
    % Check and ensure mode A jstar >= mode B jstar
    P2_AB = ensure_modeA_primary(P2_AB);
        
    if prune.enable
        [~, infoAB] = decide_single_mode_by_contrib(P2_AB, etag, maxExp, ...
            st2.newtonMaxIter, st2.newtonTol, prune);
        fprintf('  [Diagnostic] Mode A contribution=%.3g, Mode B contribution=%.3g (valid points=%d)\n', ...
            infoAB.fracA, infoAB.fracB, infoAB.N);
    end
        
    % ===================== Fit single mode A/B =====================
    fprintf('  Stage 1(A): Initial fitting without ohmic drop...\n');
    P1_A = fit_stage1_single_noohm(etag, jg, Iref, f, maxExp, st1, 'A', P2_AB);
    fprintf('  Stage 2(A): Full implicit optimization...\n');
    P2_A = fit_stage2_single_full(etag, jg, Iref, f, maxExp, st2, P1_A, 'A');
    
    fprintf('  Stage 1(B): Initial fitting without ohmic drop...\n');
    P1_B = fit_stage1_single_noohm(etag, jg, Iref, f, maxExp, st1, 'B', P2_AB);
    fprintf('  Stage 2(B): Full implicit optimization...\n');
    P2_B = fit_stage2_single_full(etag, jg, Iref, f, maxExp, st2, P1_B, 'B');
    
    % ===================== AIC model selection =====================
    if sel.enableAIC
        rssAB = rss_asinh_model(P2_AB, etag, jg, Iref, maxExp, st2.newtonMaxIter, st2.newtonTol);
        rssA = rss_asinh_model(P2_A, etag, jg, Iref, maxExp, st2.newtonMaxIter, st2.newtonTol);
        rssB = rss_asinh_model(P2_B, etag, jg, Iref, maxExp, st2.newtonMaxIter, st2.newtonTol);
            
        n = numel(etag);
        AIC_AB = aic_from_rss(rssAB, n, 9);  % AB: 9 parameters
        AIC_A = aic_from_rss(rssA, n, 5);    % Single mode: 5 parameters
        AIC_B = aic_from_rss(rssB, n, 5);
            
        bestSingleAIC = min(AIC_A, AIC_B);
        [~, idxMin] = min([AIC_AB, AIC_A, AIC_B]); % 1=AB, 2=A, 3=B
            
        choose = idxMin;
            
        if ~sel.forceAB && idxMin == 1
            if (bestSingleAIC - AIC_AB) < sel.AIC_deltaMin
                choose = 2 + (AIC_B < AIC_A); % Select better single mode
            end
        end
            
        if choose == 1
            P1 = P1_AB; P2 = P2_AB; modelTag = "AB";
        elseif choose == 2
            P1 = P1_A; P2 = P2_A; modelTag = "A";
        else
            P1 = P1_B; P2 = P2_B; modelTag = "B";
        end
            
        fprintf('  [AIC] AB=%.2f, A=%.2f, B=%.2f -> Choose %s\n', ...
            AIC_AB, AIC_A, AIC_B, modelTag);
    else
        P1 = P1_AB; P2 = P2_AB; modelTag = "AB";
    end
        
    % Store results
    RES(g) = struct('pH', ph, 'P', P2, 'P1', P1, 'Iref', Iref, ...
                    'U', Ug, 'eta', etag, 'j', jg, 'modelTag', modelTag);
        
    fprintf('  Final(%s): j*A=%.2e ξA=%.3f λA=%.3f jlimA=%.3g | ', ...
        modelTag, P2.jstarA, P2.xiA, P2.lamA, P2.jlimA);
    fprintf('j*B=%.2e ξB=%.3f λB=%.3f jlimB=%.3g | Rohm=%.3g\n', ...
        P2.jstarB, P2.xiB, P2.lamB, P2.jlimB, P2.Rohm);
end

% Open file for writing parameters
fid = fopen('optimized_parameters_Au.txt', 'w');

% Write header
fprintf(fid, 'pH\tModel\tjstarA\txiA\tlamA\tjlimA\tjstarB\txiB\tlamB\tjlimB\tRohm\n');

% Iterate through all results and write to file
for g = 1:nG
    res = RES(g);
    
    % Check P structure
    if isstruct(res.P)
        % If structure, access fields directly
        fprintf(fid, '%.2f\t%s\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\n', ...
            res.pH, res.modelTag, ...
            res.P.jstarA, res.P.xiA, res.P.lamA, res.P.jlimA, ...
            res.P.jstarB, res.P.xiB, res.P.lamB, res.P.jlimB, ...
            res.P.Rohm);
    else
        % If P is array/matrix
        fprintf('Warning: pH=%.2f P is not a structure, using backup output\n', res.pH);
        
        % Backup: adjust indices based on actual data structure
        if numel(res.P) >= 9
            fprintf(fid, '%.2f\t%s\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\n', ...
                res.pH, res.modelTag, ...
                res.P(1), res.P(2), res.P(3), res.P(4), ...  % jstarA, xiA, lamA, jlimA
                res.P(5), res.P(6), res.P(7), res.P(8), ...  % jstarB, xiB, lamB, jlimB
                res.P(9));  % Rohm
        else
            % If array length insufficient, output NaN
            fprintf(fid, '%.2f\t%s\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n', ...
                res.pH, res.modelTag);
        end
    end
end

% Close file
fclose(fid);
fprintf('\nOptimized parameters saved to: optimized_parameters_Au.txt\n');

%% ---------------- Basic Visualization ----------------
plotOpt = struct();
plotOpt.nFine = 1200;
    
% Polarization curve decomposition plot
plot_polarization_decomp(RES, plotOpt, maxExp);
    
% Tafel plot analysis
TAF = plot_tafel_decomp_and_fit(RES, plotOpt, maxExp);
    
% Total fitting comparison plot
plotOpt2 = struct();
plotOpt2.titleTotal = 'Total fitting vs data (all pH)';
plotOpt2.figFolder = 'figures';
plotOpt2.saveFigures = true;
plot_polarization_comparison_multi_pH(RES, plotOpt2, maxExp, ...
    st2.newtonMaxIter, st2.newtonTol);
    
%% ---------------- Cai-Smith Plot Analysis ----------------
caiOpt = struct();
caiOpt.useEtaEff = true;
caiOpt.phaseWeight = 'cos';
caiOpt.showGrid = true;
caiOpt.capGammaTo1 = false;
caiOpt.title = 'Cai-Smith Plot (multi-pH overlay)';
caiOpt.plotStorageCurve = true;
caiOpt.storageUseEtaEff = true;
caiOpt.storageYScale = 'log';
caiOpt.nEta = 900;
caiOpt.edgeFrac = 0.01;
caiOpt.etaPctRange = [1 99];

% Call plotting function
plot_waveguide_electrochem_summary_figure(RES, maxExp, st2.newtonMaxIter, st2.newtonTol);

opt = struct();
opt.colorMode = 'Pi_dens';      % or 'I_store' / 'rho' / 'theta'
opt.perTileColorbar = true;
opt.zoomToData = true;
opt.showUnitCircle = false;
opt.savePng = 'CS_electrochem.png';
plot_caismith_electrochem_waveguide(RES, opt);
    
%% ---------------- Waveguide Theory Advanced Analysis ----------------
fprintf('\n=== Performing waveguide theory advanced analysis ===\n');

% Resonance and attenuation dynamics analysis
analyze_resonance_dynamics(RES, maxExp, st2.newtonMaxIter, st2.newtonTol);
    
fprintf('\n=== Analysis completed ===\n');

%% =====================================================================
% Helper Functions
%% =====================================================================

function pHnum = parse_pH_column(pHcol)
    if isnumeric(pHcol)
        pHnum = double(pHcol(:));
        return;
    end
    s = string(pHcol(:));
    pHnum = nan(numel(s),1);
    for i = 1:numel(s)
        tok = regexp(s(i), '([0-9]*\.?[0-9]+)', 'tokens', 'once');
        if ~isempty(tok)
            pHnum(i) = str2double(tok{1});
        end
    end
    if any(~isfinite(pHnum))
        error('Some pH values cannot be parsed, please check pH column.');
    end
end

function [jstar2, xi2, lam2, jlim2] = seed_second_mode_from_residual(eta, dj, f, maxExp, jmax)
    % Estimate initial parameters for second mode from residual
    eta = eta(:); dj = dj(:);
    
    absr = abs(dj);
    thr = prctile_fallback(absr, 80);
    mask = isfinite(eta) & isfinite(dj) & (absr >= thr);
    if nnz(mask) < 10
        mask = isfinite(eta) & isfinite(dj);
    end
    
    e = eta(mask);
    d = dj(mask);
    
    xi2  = 0.35;
    lam2 = 0.40;
    
    a = xi2*lam2*f;
    b = -(1-xi2)*lam2*f;
    
    Ap = clip(a.*e, maxExp);
    Am = clip(b.*e, maxExp);
    phi = exp(Ap) - exp(Am);
    phi(abs(phi) < 1e-12) = 1e-12;
    
    jstar_raw = d ./ phi;
    jstar2 = median(abs(jstar_raw(isfinite(jstar_raw))));
    jstar2 = max(jstar2, 1e-10);
    jstar2 = min(jstar2, max(0.2*jmax, 1e-3));
    
    jlim2 = max(0.8*jmax, 1e-3);
end

%% =====================================================================
% Basic Visualization Functions
%% =====================================================================

function plot_polarization_decomp(RES, plotOpt, maxExp)
% Set global defaults
set(groot, 'DefaultAxesFontName', 'Helvetica');
set(groot, 'DefaultTextFontName', 'Helvetica');
set(groot, 'DefaultAxesFontSize', 12);          % tick labels
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.15);  % xlabel/ylabel
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.25);  % title
set(groot, 'DefaultLegendFontSize', 12);
set(groot, 'DefaultAxesLineWidth', 1.0);

    nG = numel(RES);
    C = lines(nG);
    
    allJ = vertcat(RES.j);
    jmag = abs(allJ(allJ~=0));
    if isempty(jmag), jmag = 1; end
    yL = max(prctile_fallback(jmag, 99.5), max(jmag));
    
    ncol = 4;
    nrow = ceil(nG/ncol);
    width_cm = 22.5;
    height_cm = 18;
    fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
    tlo = tiledlayout(fig,nrow, ncol,'Padding','compact','TileSpacing','compact');
    
    % Initialize array to store RMSE
    RMSE_values = zeros(nG, 1);
    
    for g = 1:nG
        nexttile;
        hold on; box on; grid on;
        
        P = RES(g).P;
        etag = RES(g).eta(:);
        jg   = RES(g).j(:);
        
        etaFine = linspace(min(etag), max(etag), plotOpt.nFine).';
        jTot = safe_two_mode_model(P, etaFine, maxExp, 10, 5e-11, []);
        [jA, jB] = two_mode_decompose_given_total(P, etaFine, jTot, maxExp);
        
        plot(etag, jg, 'ko', 'MarkerSize', 3.0, 'LineWidth', 1.0);
        plot(etaFine, jTot, 'r-', 'LineWidth', 2.5);
        plot(etaFine, jA, 'b--', 'LineWidth', 2.0);
        plot(etaFine, jB, 'g:', 'LineWidth', 2.0);
        
        % Calculate RMSE
        j_pred = safe_two_mode_model(P, etag, maxExp, 10, 5e-11, []);
        RMSE = sqrt(mean((j_pred - jg).^2));
        RMSE_values(g) = RMSE;
        
        % Format RMSE to 2 significant figures
        if RMSE >= 100
            RMSE_str = sprintf('%.2f', RMSE);
        elseif RMSE >= 10
            RMSE_str = sprintf('%.2f', RMSE);
        elseif RMSE >= 1
            RMSE_str = sprintf('%.2f', RMSE);
        elseif RMSE >= 0.1
            RMSE_str = sprintf('%.2f', RMSE);
        elseif RMSE >= 0.01
            RMSE_str = sprintf('%.2f', RMSE);
        else
            RMSE_str = sprintf('%.2f', RMSE);
        end
        
        xlabel('$\eta_\mathrm{eff}$ (V)', 'Interpreter','latex','FontSize',12);
        ylabel('$j$ (mA cm$^{-2}$)', 'Interpreter','latex','FontSize',12);
        
        ylim([-1.15*yL, 0.15*yL]);
        
        if g == 1
            legend1=legend({'Data','$j_\mathrm{tot}','$j_\mathrm{A}','$j_\mathrm{B}'}, 'Location',...
                'southeast','FontSize',12,'Interpreter','latex');
            legend boxoff
            set(legend1,...
             'Position',[0.129032822837497 0.564893662799677 0.187843137254902 0.11572712542964]);
        end
        
        % Add label (a, b, c, ...) to each subplot
        text(-0.15, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end
    
    % If subplot count < 12, show RMSE statistics in last subplot
    if nG < 12
        % Determine last subplot position
        last_tile = nG;
        
        % Get last subplot axis
        ax_last = nexttile(last_tile);
        
        % Add RMSE statistics to current axis
        hold(ax_last, 'on');
        
        % Calculate statistics
        mean_RMSE = mean(RMSE_values);
        min_RMSE = min(RMSE_values);
        max_RMSE = max(RMSE_values);
        std_RMSE = std(RMSE_values);
        
        % Format values to 2 significant figures
        formatValue = @(x) sprintf('%.2f', x);
        
        % Create statistics text
        stats_text = {};
        
        % Add RMSE for each pH
        for g = 1:nG
            stats_text{end+1} = sprintf('  pH=%.1f: %s mA/cm²', RES(g).pH, formatValue(RMSE_values(g)));
        end
        
        % Add text to right side of axis
        text(ax_last, 1.2, 0.98, stats_text, ...
            'Units','normalized', 'FontSize', 12, 'FontWeight','normal', ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', ...
            'BackgroundColor', [1 1 1 0.9], 'EdgeColor', 'k', 'Margin', 3);
        
    else
        % If 8 or more subplots, show statistics in subplot 8
        stats_tile = min(8, nG);  % Use subplot 8, or last if nG<8
        
        % Get subplot axis
        ax_stats = nexttile(stats_tile);
        
        % Clear current axis
        cla(ax_stats);
        set(ax_stats, 'Visible', 'off');
        
        % Calculate statistics
        mean_RMSE = mean(RMSE_values);
        min_RMSE = min(RMSE_values);
        max_RMSE = max(RMSE_values);
        std_RMSE = std(RMSE_values);
        
        % Format values to 2 significant figures
        formatValue = @(x) sprintf('%.2g', x);
        
        % Create statistics text
        stats_text = {
            'RMSE Statistics Summary';
            '';
            sprintf('Mean: %s mA/cm²', formatValue(mean_RMSE));
            sprintf('Min: %s mA/cm²', formatValue(min_RMSE));
            sprintf('Max: %s mA/cm²', formatValue(max_RMSE));
            sprintf('Std: %s mA/cm²', formatValue(std_RMSE));
            '';
            'RMSE per pH:'
        };
        
        % Add RMSE for each pH
        for g = 1:nG
            stats_text{end+1} = sprintf('  pH=%.0f: %s mA/cm²', RES(g).pH, formatValue(RMSE_values(g)));
        end
        
        % Add label to statistics subplot
        text(ax_stats, -0.20, 0.98, char('a' + stats_tile - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end
    
    % Adjust all subplot proportions
    set(tlo, 'TileSpacing', 'compact', 'Padding', 'compact');
    filename = 'FigureS20';
    print(gcf, '-dpng', '-r600', [filename, '.png']);
    disp(['Figure saved as ', filename, '.png']);
end

function plot_polarization_comparison_multi_pH(RES, plotOpt, maxExp, newtonMaxIter, newtonTol)
    if nargin < 2 || isempty(plotOpt), plotOpt = struct(); end
    plotOpt = set_default(plotOpt, 'nFine', 700);
    plotOpt = set_default(plotOpt, 'saveFigures', true);
    plotOpt = set_default(plotOpt, 'figFolder', 'figures');
    plotOpt = set_default(plotOpt, 'titleTotal', 'Total fitting vs data (all pH)');
    
    nG = numel(RES);
    C  = lines(nG);
    
    % Total overlay plot
    fig1 = figure('Color','w', 'Name','Total fitting overlay', 'Position',[90 70 980 700]);
    ax = axes(fig1);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    
    for g = 1:nG
        eta = RES(g).eta(:);
        j = RES(g).j(:);
        P = RES(g).P;
        
        [eta, ord] = sort(eta);
        j = j(ord);
        etaFine = linspace(min(eta), max(eta), plotOpt.nFine).';
        jTot = safe_two_mode_model(P, etaFine, maxExp, newtonMaxIter, newtonTol, []);
        
        scatter(ax, eta, j, 18, 'MarkerEdgeColor', C(g,:), ...
            'MarkerFaceColor', C(g,:), 'MarkerFaceAlpha', 0.55, ...
            'HandleVisibility', 'off');
        plot(ax, etaFine, jTot, '-', 'Color', C(g,:), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('pH=%.3g (%s)', RES(g).pH, RES(g).modelTag));
    end
    
    xlabel(ax, '$\eta$ (V)', 'Interpreter','latex');
    ylabel(ax, '$j$ (mA cm$^{-2}$)', 'Interpreter','latex');
    title(ax, plotOpt.titleTotal, 'Interpreter','none');
    legend(ax, 'Location', 'eastoutside');
    yline(ax, 0, 'k-', 'HandleVisibility','off');
    xline(ax, 0, 'k-', 'HandleVisibility','off');
    
    if plotOpt.saveFigures
        save_figure(fig1, 'polarization_total_multi_pH', plotOpt.figFolder);
    end
end

function TAF = plot_tafel_decomp_and_fit(RES, plotOpt, maxExp)
    plotOpt = set_default(plotOpt, 'nFine', 900);
    plotOpt = set_default(plotOpt, 'newtonMaxIter', 12);
    plotOpt = set_default(plotOpt, 'newtonTol', 5e-11);
    plotOpt = set_default(plotOpt, 'useEtaEff', true);
    
    cfg = struct();
    cfg.f_lo       = 0.03;
    cfg.f_hi       = 0.30;
    cfg.minN       = 12;
    cfg.minDecades = 0.6;
    cfg.R2_min     = 0.90;
    cfg.branch     = 'cathodic';
    
    nG = numel(RES);
    C  = lines(nG);
    
    pH_vals = nan(nG,1);
    betaA   = nan(nG,1); R2A = nan(nG,1); NA = nan(nG,1);
    betaB   = nan(nG,1); R2B = nan(nG,1); NB = nan(nG,1);
    
    ncol = 4;
    nrow = ceil(nG/ncol);
    width_cm = 22.5;
    height_cm = 18;
    fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
    tlo = tiledlayout(fig,nrow, ncol,'Padding','compact','TileSpacing','compact');
    
    for g = 1:nG
        nexttile;
        hold on; box on; grid on;
        
        etag = RES(g).eta(:);
        [etag, ~] = sort(etag);
        etaFine = linspace(min(etag), max(etag), plotOpt.nFine).';
        
        P = RES(g).P;
        jTot = safe_two_mode_model(P, etaFine, maxExp, plotOpt.newtonMaxIter, plotOpt.newtonTol, []);
        [jA, jB] = two_mode_decompose_given_total(P, etaFine, jTot, maxExp);
        Rohm = P.Rohm;
        
        if plotOpt.useEtaEff
            xTot = etaFine - Rohm*(jTot./1000);
            xA   = etaFine - Rohm*(jA./1000);
            xB   = etaFine - Rohm*(jB./1000);
            xlab = '$\eta_{\mathrm{eff}}$ (V)';
        else
            xTot = etaFine; xA = etaFine; xB = etaFine;
            xlab = '$\eta$ (V)';
        end
        
        yTot = safe_log10abs(jTot);
        yA   = safe_log10abs(jA);
        yB   = safe_log10abs(jB);
        
        plot(xTot, yTot, 'r-', 'LineWidth', 2.2);
        plot(xA, yA, 'k--', 'LineWidth', 1.8);
        plot(xB, yB, 'b:', 'LineWidth', 1.8);
        
        FA = fit_tafel_relative_threshold(xA, jA, cfg);
        FB = fit_tafel_relative_threshold(xB, jB, cfg);
        
        text(-0.2, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
        
        % Show Tafel slope for Mode A/B
        showR2N = false;   % Set false to show only slope
        
        if FA.ok
            plot(FA.xLine, FA.yLine, 'g--', 'LineWidth', 2.5);
            if showR2N
                txtA = sprintf('A: %.0f mV/dec  (R^2=%.2f, N=%d)', FA.beta_mVdec, FA.R2, FA.N);
            else
                txtA = sprintf('A: %.0f mV/dec', FA.beta_mVdec);
            end
        else
            txtA = 'A: n/a';
        end

        if FB.ok
            plot(FB.xLine, FB.yLine, 'g:', 'LineWidth', 2.5);
            if showR2N
                txtB = sprintf('B: %.0f mV/dec  (R^2=%.2f, N=%d)', FB.beta_mVdec, FB.R2, FB.N);
            else
                txtB = sprintf('B: %.0f mV/dec', FB.beta_mVdec);
            end
        else
            txtB = 'B: n/a';
        end

        % Place in top left, two lines
        x0 = 0.02; y0 = 0.20; dy = 0.08;  % Normalized coordinates
        text(x0, y0, txtA, 'Units','normalized', 'HorizontalAlignment','left', ...
            'VerticalAlignment','top', 'Color',[0.85 0 0], 'FontSize',10, 'FontWeight','bold');
        text(x0, y0-dy, txtB, 'Units','normalized', 'HorizontalAlignment','left', ...
            'VerticalAlignment','top', 'Color',[0 0.25 0.9], 'FontSize',10, 'FontWeight','bold');
        
        xlabel(xlab, 'Interpreter','latex');
        ylabel('$\log_{10}|j|$', 'Interpreter','latex');
        ylim([-8 2]);
        
        if g == 1
            legend1=legend({'$j_\mathrm{tot}$','$j_\mathrm{A}$','$j_\mathrm{B}$',...
                'Tafel $j_\mathrm{A}','Tafel $j_\mathrm{B}'}, 'Location','southwest',...
                'Interpreter','latex','FontSize',12);
            legend boxoff
            set(legend1,...
                'Position',[0.629472135195725 0.112193774957868 0.079556626539964 0.15634765625]);
        end
        
        pH_vals(g) = RES(g).pH;
    end
    
    TAF = table(pH_vals, betaA, R2A, NA, betaB, R2B, NB, ...
        'VariableNames', {'pH','betaA_mVdec','R2A','NA','betaB_mVdec','R2B','NB'});
    filename = 'FigureS21';
    print(gcf, '-dpng', '-r600', [filename, '.png']);
    disp(['Figure saved as ', filename, '.png']);
end

%% =====================================================================
% Cai-Smith Plot Analysis (Improved)
%% =====================================================================

function CaiTAB = plot_caismith_multi_pH_overlay_improved(RES, caiOpt, maxExp, newtonMaxIter, newtonTol)
    % Set default parameters
    caiOpt = set_default(caiOpt, 'useEtaEff', true);
    caiOpt = set_default(caiOpt, 'phaseWeight', 'cos');
    caiOpt = set_default(caiOpt, 'showGrid', true);
    caiOpt = set_default(caiOpt, 'capGammaTo1', false);
    caiOpt = set_default(caiOpt, 'title', 'Cai-Smith Plot (multi-pH overlay)');
    caiOpt = set_default(caiOpt, 'plotStorageCurve', true);
    caiOpt = set_default(caiOpt, 'storageUseEtaEff', true);
    caiOpt = set_default(caiOpt, 'storageYScale', 'log');
    caiOpt = set_default(caiOpt, 'nEta', 900);
    caiOpt = set_default(caiOpt, 'edgeFrac', 0.01);
    caiOpt = set_default(caiOpt, 'etaPctRange', [1 99]);
    
    nG = numel(RES);
    C  = lines(nG);
    LS = {'-','--',':','-.'};
    
    % Update CUR structure to include all new parameters
    CUR = repmat(struct('pH',[], 'eta',[], 'eta_eff',[], 'Gmag',[], ...
        'Gph',[], 'I_store',[], 'eta_trans',[], 'U',[], 'S_rel',[], ...
        'beta_w',[], 'phi',[], 'D_plus',[], 'D_minus',[], ...
        'I_store_max',[], 'eta_trans_max',[], 'imax',[], ...
        'I_store_int',[], 'eta_trans_int',[]), nG, 1);
    
    % Calculate Cai-Smith parameters for each pH
    for g = 1:nG
        eta_raw = RES(g).eta(:);
        lo = prctile_fallback(eta_raw, caiOpt.etaPctRange(1));
        hi = prctile_fallback(eta_raw, caiOpt.etaPctRange(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        etaGrid = linspace(lo, hi, caiOpt.nEta).';
        
        P = RES(g).P;
        
        % Use improved function to calculate all parameters
        [Gmag, Gph, I_store, eta_eff, U, S_rel, beta_w, phi, D_plus, D_minus, eta_trans] = ...
            calc_caismith_from_params_improved(...
                P, etaGrid, caiOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, caiOpt.phaseWeight);
        
        if caiOpt.capGammaTo1
            Gmag = min(Gmag, 1);
        end
        
        % Use I_store as storage intensity indicator
        imax = pick_peak_interior(I_store, caiOpt.edgeFrac);
        
        % Store all parameters in CUR structure
        CUR(g).pH = RES(g).pH;
        CUR(g).eta = etaGrid;
        CUR(g).eta_eff = eta_eff;
        CUR(g).Gmag = Gmag;
        CUR(g).Gph = Gph;
        CUR(g).I_store = I_store;
        CUR(g).eta_trans = eta_trans;
        CUR(g).U = U;
        CUR(g).S_rel = S_rel;
        CUR(g).beta_w = beta_w;
        CUR(g).phi = phi;
        CUR(g).D_plus = D_plus;
        CUR(g).D_minus = D_minus;
        CUR(g).imax = imax;
        CUR(g).I_store_max = I_store(imax);
        CUR(g).eta_trans_max = eta_trans(imax);
        CUR(g).I_store_int = trapz(eta_eff, I_store);
        CUR(g).eta_trans_int = trapz(eta_eff, eta_trans);
    end
    
    % Determine angle range
    all_deg = [];
    for g = 1:nG
        deg = CUR(g).Gph * 180/pi;
        deg = mod(deg + 180, 360) - 180;
        all_deg = [all_deg; deg(:)];
    end
    all_deg = all_deg(isfinite(all_deg));
    if isempty(all_deg)
        thMin = -90; thMax = 90;
    else
        thMin = min(all_deg);
        thMax = max(all_deg);
        pad = 0.12 * (thMax - thMin + 1e-6);
        thMin = max(thMin - pad, -90);
        thMax = min(thMax + pad, 90);
    end
    
    % ---------- Figure A: Cai-Smith disk overlay ----------
    figure('Color','w', 'Name','Cai-Smith plot overlay', 'Position',[120 60 980 820]);
    pax = polaraxes;
    hold(pax, 'on');
    
    if caiOpt.showGrid
        for r = [0.2 0.4 0.6 0.8 1.0]
            th = linspace(0, 2*pi, 240);
            h = polarplot(pax, th, r*ones(size(th)), 'k:', 'LineWidth', 0.5);
            h.HandleVisibility = 'off';
        end
        for thd = [-60 -45 -30 -15 0 15 30 45 60]
            th = thd * pi/180;
            rr = linspace(0, 1, 60);
            h = polarplot(pax, th*ones(size(rr)), rr, 'k:', 'LineWidth', 0.35);
            h.HandleVisibility = 'off';
            if thd ~= 0
                text(pax, th, 1.06, sprintf('%d°', thd), 'FontSize', 12, ...
                    'HorizontalAlignment', 'center', 'Interpreter', 'none');
            end
        end
    end
    
    for g = 1:nG
        th = CUR(g).Gph;
        rr = CUR(g).Gmag;
        ls = LS{mod(g-1, numel(LS)) + 1};
        
        polarplot(pax, th, rr, 'LineStyle', ls, 'Color', C(g,:), ...
            'LineWidth', 2.2, 'HandleVisibility', 'off');
        
        im = CUR(g).imax;
        polarplot(pax, th(im), rr(im), 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', C(g,:), 'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.0, 'DisplayName', ...
            sprintf('pH=%.3g (I_store_max=%.2e, η*=%.3g)', CUR(g).pH, CUR(g).I_store_max, CUR(g).eta_eff(im)));
    end
    
    pax.ThetaLim = [thMin thMax];
    rlim(pax, [0 1.08]);
    grid(pax, 'on');
    title(pax, caiOpt.title, 'Interpreter', 'none');
    legend(pax, 'Location', 'southoutside', 'Orientation', 'vertical', 'FontSize', 10);
    
    % ---------- Figure B: Storage intensity curve ----------
    if caiOpt.plotStorageCurve
        useEff = caiOpt.storageUseEtaEff;
        
        figure('Color','w', 'Name','Storage intensity vs overpotential', 'Position',[140 80 1060 760]);
        tl = tiledlayout(2, 1, 'Padding','compact', 'TileSpacing','compact');
        
        ax1 = nexttile(tl, 1);
        hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
        for g = 1:nG
            if useEff
                x = CUR(g).eta_eff;
            else
                x = CUR(g).eta;
            end
            I_store_g = CUR(g).I_store;
            I_store_norm = I_store_g / (max(I_store_g) + 1e-300);
            
            plot(ax1, x, I_store_norm, 'LineWidth', 2.0, 'Color', C(g,:), ...
                'LineStyle', LS{mod(g-1, numel(LS))+1}, ...
                'DisplayName', sprintf('pH=%.3g', CUR(g).pH));
            
            im = CUR(g).imax;
            plot(ax1, x(im), I_store_norm(im), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', C(g,:), 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        if useEff
            xlab = '$\eta_{\mathrm{eff}}$ (V)';
        else
            xlab = '$\eta$ (V)';
        end
        xlabel(ax1, xlab, 'Interpreter', 'latex');
        ylabel(ax1, '$S_{rel}$', 'Interpreter', 'latex');
        title(ax1, 'Normalized integrated storage intensity vs overpotential (normalized by pH)', 'Interpreter', 'none');
        legend(ax1, 'Location', 'eastoutside');
        
        ax2 = nexttile(tl, 2);
        hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
        for g = 1:nG
            if useEff
                x = CUR(g).eta_eff;
            else
                x = CUR(g).eta;
            end
            I_store_g= CUR(g).I_store;
            
            if strcmpi(caiOpt.storageYScale, 'log')
                semilogy(ax2, x, max(I_store_g, 1e-300), 'LineWidth', 2.0, ...
                    'Color', C(g,:), 'LineStyle', LS{mod(g-1, numel(LS))+1});
            else
                plot(ax2, x, I_store_g, 'LineWidth', 2.0, 'Color', C(g,:), ...
                    'LineStyle', LS{mod(g-1, numel(LS))+1});
            end
            
            im = CUR(g).imax;
            ymark = I_store_g(im);
            if strcmpi(caiOpt.storageYScale, 'log')
                ymark = max(ymark, 1e-300);
            end
            plot(ax2, x(im), ymark, 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', C(g,:), 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        if useEff
            xlab = '$\eta_{\mathrm{eff}}$ (V)';
        else
            xlab = '$\eta$ (V)';
        end
        xlabel(ax2, xlab, 'Interpreter', 'latex');
        
        if strcmpi(caiOpt.storageYScale, 'log')
            ylabel(ax2, '$I_{\mathrm{store}}$ (absolute value, log)', 'Interpreter', 'latex');
            title(ax2, 'Absolute storage intensity vs overpotential (semi-log)', 'Interpreter', 'none');
        else
            ylabel(ax2, '$I_{\mathrm{store}}$ (absolute value)', 'Interpreter', 'latex');
            title(ax2, 'Absolute storage intensity vs overpotential', 'Interpreter', 'none');
        end
        
        sgtitle(tl, 'Storage intensity comparison across different pH', 'Interpreter', 'none');
    end
    
    % Output table - updated to include complete parameters
    pH = nan(nG,1);
    I_store_max = nan(nG,1);
    I_store_int = nan(nG,1);
    eta_trans_max = nan(nG,1);
    eta_trans_int = nan(nG,1);
    etaEff_at_I_store_max = nan(nG,1);
    GammaMag_at_I_store_max = nan(nG,1);
    GammaPhase_deg_at_I_store_max = nan(nG,1);
    mean_beta_w = nan(nG,1);
    mean_phi = nan(nG,1);
    
    for g = 1:nG
        pH(g) = CUR(g).pH;
        I_store_max(g) = CUR(g).I_store_max;
        I_store_int(g) = CUR(g).I_store_int;
        eta_trans_max(g) = CUR(g).eta_trans_max;
        eta_trans_int(g) = CUR(g).eta_trans_int;
        im = CUR(g).imax;
        etaEff_at_I_store_max(g) = CUR(g).eta_eff(im);
        GammaMag_at_I_store_max(g) = CUR(g).Gmag(im);
        deg = CUR(g).Gph(im) * 180/pi;
        GammaPhase_deg_at_I_store_max(g) = mod(deg + 180, 360) - 180;
        mean_beta_w(g) = mean(CUR(g).beta_w(isfinite(CUR(g).beta_w)));
        mean_phi(g) = mean(CUR(g).phi(isfinite(CUR(g).phi)));
    end
    
    % Create complete Cai-Smith analysis results table
    CaiTAB = table(pH, I_store_max, I_store_int, eta_trans_max, eta_trans_int, ...
        etaEff_at_I_store_max, GammaMag_at_I_store_max, GammaPhase_deg_at_I_store_max, ...
        mean_beta_w, mean_phi, ...
        'VariableNames', {'pH', 'I_store_max', 'I_store_int', 'eta_trans_max', ...
        'eta_trans_int', 'etaEff_at_I_store_max', 'GammaMag_at_I_store_max', ...
        'GammaPhase_deg_at_I_store_max', 'mean_beta_w', 'mean_phi'});
    
    % Backward compatibility: also output old format table
    if nargout > 0
        % Keep old variable names for compatibility
        Smax_abs = I_store_max;
        S_int_abs = I_store_int;
        etaEff_at_Smax = etaEff_at_I_store_max;
        GammaMag_at_Smax = GammaMag_at_I_store_max;
        GammaPhase_deg_at_Smax = GammaPhase_deg_at_I_store_max;
        
        % Output backward compatible table
        CaiTAB_legacy = table(pH, Smax_abs, S_int_abs, etaEff_at_Smax, ...
            GammaMag_at_Smax, GammaPhase_deg_at_Smax, ...
            'VariableNames', {'pH', 'Smax_abs', 'S_int_abs', ...
            'etaEff_at_Smax', 'GammaMag_at_Smax', 'GammaPhase_deg_at_Smax'});
        
        fprintf('\n=== Cai-Smith analysis results (new format) ===\n');
        disp(CaiTAB);
        
        fprintf('\n=== Cai-Smith analysis results (old format - backward compatible) ===\n');
        disp(CaiTAB_legacy);
    end
end

function [Gmag, Gph, I_store, eta_eff, U, S_rel, beta_w, rapidity, D_plus, D_minus, eta_trans] = ...
    calc_caismith_from_params_improved(P, etaGrid, useEtaEff, maxExp, newtonMaxIter, newtonTol, phaseWeight)
% Waveguide-state diagnostics from polarization only (updated to Sec.2.5 & Methods)
% Convention: b_p -> cathodic reduction branch, a_p -> anodic oxidation branch.
% Use UNSATURATED directional fluxes for Gamma_g phase/amplitude; saturation only affects net j via fitting.

    if nargin < 7 || isempty(phaseWeight), phaseWeight = 'none'; end

    etaGrid = etaGrid(:);
    eps0 = 1e-30;

    % ---- (1) eta_eff from fitted implicit j(eta) ----
    if useEtaEff
        jTot = safe_two_mode_model(P, etaGrid, maxExp, newtonMaxIter, newtonTol, []);
        eta_eff = etaGrid - P.Rohm .* (jTot./1000);   % jTot: mA/cm^2 -> A/cm^2
    else
        eta_eff = etaGrid;
    end

    % ---- (2) Directional fluxes from UNSATURATED kinetics ----
    % J_ox,p = jstar_p * exp(a_p * eta_eff)   (anodic oxidation)
    % J_red,p = jstar_p * exp(b_p * eta_eff)  (cathodic reduction)
    % NOTE: do NOT apply wg_mass_energy_map_optimized() to these components.
    JoxA  = P.jstarA .* exp( clip(P.aA .* eta_eff, maxExp) );
    JredA = P.jstarA .* exp( clip(P.bA .* eta_eff, maxExp) );

    JoxB  = P.jstarB .* exp( clip(P.aB .* eta_eff, maxExp) );
    JredB = P.jstarB .* exp( clip(P.bB .* eta_eff, maxExp) );

    Jox  = JoxA  + JoxB;
    Jred = JredA + JredB;

    % ---- (3) Power-flow & exchange proxies ----
    U = Jox + Jred;          % exchange proxy  (>=0)
    S_rel = Jox - Jred;      % transfer proxy  (signed)
    beta_w = S_rel ./ (U + eps0);     % in [-1,1]

    % numerical safety for atanh
    beta_w = max(min(beta_w, 1-1e-12), -1+1e-12);

    rapidity = atanh(beta_w);
    D_plus  = exp(rapidity);
    D_minus = exp(-rapidity);

    % storage intensity (normalized): m/U = 2*sqrt(Jox*Jred)/(Jox+Jred)
    I_store = 2 * sqrt(max(Jox,0).*max(Jred,0)) ./ (U + eps0);
    I_store = max(min(I_store, 1), 0);   % [0,1]
    eta_trans = abs(beta_w);             % [0,1]

    % ---- (4) Polarization-derived coherence angle and Gamma_g ----
    chi = (Jred - Jox) ./ (Jred + Jox + eps0);     % in [-1,1]
    chi = max(min(chi, 1), -1);

    varphi = acos(chi);                 % [0, pi]
    theta  = varphi - pi/2;             % [-pi/2, +pi/2] for Cai-Smith disk display

    rho = sqrt( min(Jred, Jox) ./ (max(Jred, Jox) + eps0) );  % [0,1]
    rho(~isfinite(rho)) = 0;

    % complex Gamma_g
    Gamma_complex = rho .* exp(1i*theta);

    Gmag = abs(Gamma_complex);          % = rho
    Gph  = angle(Gamma_complex);        % = theta in [-pi/2, pi/2]

    % cleanup
    Gmag(~isfinite(Gmag)) = 0;
    Gph(~isfinite(Gph))   = 0;
    U(~isfinite(U)) = eps0;
    S_rel(~isfinite(S_rel)) = 0;
    beta_w(~isfinite(beta_w)) = 0;
    rapidity(~isfinite(rapidity)) = 0;
    D_plus(~isfinite(D_plus)) = 1;
    D_minus(~isfinite(D_minus)) = 1;
    I_store(~isfinite(I_store)) = 0;
    eta_trans(~isfinite(eta_trans)) = 0;
end

function imax = pick_peak_interior(S, frac)
    S = S(:);
    n = numel(S);
    if n < 5
        [~, imax] = max(S);
        return;
    end
    k0 = max(2, round(frac*n));
    k1 = min(n-1, round((1-frac)*n));
    if k1 <= k0
        [~, imax] = max(S);
        return;
    end
    [~, ii] = max(S(k0:k1));
    imax = ii + k0 - 1;
end

%% =====================================================================
% Waveguide Theory Advanced Analysis Functions
%% =====================================================================

function analyze_resonance_dynamics(RES, maxExp, newtonMaxIter, newtonTol)
% analyze_resonance_dynamics (TWO-FIG VERSION, with legend in tile-8)
% ------------------------------------------------------------
% Fig-1: feedback & retention (rho + tau_ec), legend in tile 8
% Fig-2: storage, sensitivity, optimum, legend in tile 8
%   Right axis now: |d chi / d eta_eff| (instead of |d arg Gamma / d eta_eff|)

    fprintf('\n=== Resonance & attenuation diagnostics (Gamma_g-based) ===\n');

    if nargin < 2 || isempty(maxExp),        maxExp = 70; end
    if nargin < 3 || isempty(newtonMaxIter), newtonMaxIter = 10; end
    if nargin < 4 || isempty(newtonTol),     newtonTol = 5e-11; end

    eta0 = 0.05;
    fixedEtaWindow = [-1.5, 1.0];
    nEta = 700;

    nG = numel(RES);
    nShow = min(10, nG);

    OUT = repmat(struct('pH',NaN,'eta_eff_s',[],'rho_s',[],'tau_ec',[], ...
                        'Istore_s',[],'dchi_abs',[],'eta_eff2s',[],'Pi_dens',[], ...
                        'maskC',[],'etaStar',NaN), nShow, 1);

    for g = 1:nShow
        P = RES(g).P;

        etaGrid = linspace(fixedEtaWindow(1), fixedEtaWindow(2), nEta).';

        % Get rho, theta, I_store, eta_eff quickly
        [rho, theta, I_store, eta_eff] = ...
            calc_caismith_from_params_improved(P, etaGrid, true, maxExp, newtonMaxIter, newtonTol, 'none');

        rho     = max(min(rho(:),1),0);
        theta   = unwrap(theta(:));
        I_store = max(min(I_store(:),1),0);
        eta_eff = eta_eff(:);

        % sort by eta_eff
        [eta_eff_s, ord] = sort(eta_eff);
        rho_s    = rho(ord);
        Istore_s = I_store(ord);

        % attenuation proxy
        tau_ec = 1 ./ max(-log(rho_s + 1e-12), 1e-12);

        % ---------------- NEW: chi and |dchi/deta_eff| ----------------
        % compute chi on the SAME eta_eff grid using UNSATURATED directional fluxes
        chi = compute_chi_from_params_unsat(P, eta_eff, maxExp);
        chi_s = chi(ord);
        dchi = gradient(chi_s, eta_eff_s);
        dchi_abs = abs(dchi);

        % --- compute jTot for Pi_use/Pi_dens ---
        jTot = safe_two_mode_model(P, etaGrid, maxExp, newtonMaxIter, newtonTol, []);
        jTot = jTot(:);
        eta_eff2 = etaGrid - P.Rohm.*(jTot./1000);

        [eta_eff2s, ord2] = sort(eta_eff2(:));
        jTot_s = jTot(ord2);

        rho_interp = interp1(eta_eff_s, rho_s, eta_eff2s, 'linear', 'extrap');
        rho_interp = max(min(rho_interp,1),0);

        Pi_use  = max(0, -jTot_s) .* max(1 - rho_interp.^2, 0);
        Pi_dens = Pi_use ./ (abs(eta_eff2s) + eta0);

        maskC = (eta_eff2s < 0) & isfinite(Pi_dens) & isfinite(eta_eff2s);
        etaStar = NaN;
        if any(maskC)
            [~, im] = max(Pi_dens(maskC));
            idx = find(maskC);
            etaStar = eta_eff2s(idx(im));
        end

        OUT(g).pH = RES(g).pH;
        OUT(g).eta_eff_s = eta_eff_s;
        OUT(g).rho_s = rho_s;
        OUT(g).tau_ec = tau_ec;
        OUT(g).Istore_s = Istore_s;
        OUT(g).dchi_abs = dchi_abs;
        OUT(g).eta_eff2s = eta_eff2s;
        OUT(g).Pi_dens = Pi_dens;
        OUT(g).maskC = maskC;
        OUT(g).etaStar = etaStar;

        fprintf('pH=%.3g: eta_eff* (Pi_dens max, cathodic) = %.4f V\n', RES(g).pH, etaStar);
    end

    % ============================================================
    % FIGURE 1: feedback & retention
    % ============================================================
    ncol = 4;
    nrow = 3;
    fig1 = figure('Color','w','Units','centimeters','Position',[2 2 22.5 18]);
    tlo1 = tiledlayout(fig1, nrow, ncol, 'Padding','compact', 'TileSpacing','compact');

    for g = 1:nShow
        ax = nexttile(tlo1, g);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        eta_eff_s = OUT(g).eta_eff_s;
        rho_s     = OUT(g).rho_s;
        tau_ec    = OUT(g).tau_ec;
        etaStar   = OUT(g).etaStar;

        yyaxis(ax,'left');
        plot(ax, eta_eff_s, rho_s, 'b-', 'LineWidth', 2);
        ylabel(ax,'$\rho=|\Gamma_g|$','Interpreter','latex');
        ylim(ax,[0 1.05]);

        yyaxis(ax,'right');
        semilogy(ax, eta_eff_s, tau_ec, 'r-', 'LineWidth', 2);
        ylabel(ax,'$\tilde{\tau}_{\rm ec}=1/(-\ln\rho)$','Interpreter','latex');

        xlabel(ax,'$\eta_{\mathrm{eff}}$ (V)','Interpreter','latex');

        text(ax, -0.25, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end

    % ---- tile 11 legend (Fig 1) ----
    axL1 = nexttile(tlo1, 11);
    axis(axL1,'off'); hold(axL1,'on');

    h1 = plot(axL1, NaN, NaN, 'b-', 'LineWidth', 2);
    h2 = plot(axL1, NaN, NaN, 'r-', 'LineWidth', 2);

    lg1 = legend(axL1, [h1 h2], ...
        {'$\rho=|\Gamma_g|$', '$\tilde{\tau}_{\rm ec}=1/(-\ln\rho)$'}, ...
        'Interpreter','latex', 'Location','best');
    lg1.Box = 'off';
    lg1.FontSize = 12;

    filename = 'FigureS23';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);

    % ============================================================
    % FIGURE 2: storage, sensitivity, optimum
    % ============================================================
    fig2 = figure('Color','w','Units','centimeters','Position',[2 2 22.5 18]);
    tlo2 = tiledlayout(fig2, nrow, ncol, 'Padding','compact', 'TileSpacing','compact');

    for g = 1:nShow
        ax = nexttile(tlo2, g);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        eta_eff_s = OUT(g).eta_eff_s;
        Istore_s  = OUT(g).Istore_s;
        dchi_abs  = OUT(g).dchi_abs;

        eta_eff2s = OUT(g).eta_eff2s;
        Pi_dens   = OUT(g).Pi_dens;
        maskC     = OUT(g).maskC;
        etaStar   = OUT(g).etaStar;

        yyaxis(ax,'left');
        plot(ax, eta_eff_s, Istore_s, 'g-', 'LineWidth', 2);
        ylabel(ax,'$I_{\rm store}$','Interpreter','latex');
        ylim(ax,[0 1.05]);

        yyaxis(ax,'right');
        semilogy(ax, eta_eff_s, max(dchi_abs, 1e-12), 'm-', 'LineWidth', 2);
        ylabel(ax,'$|\mathrm{d}\beta_{w}/\mathrm{d}\eta_{\rm eff}|$','Interpreter','latex');

        % overlay normalized Pi_dens on left axis
        yyaxis(ax,'left');
        y2 = Pi_dens; y2(~maskC) = NaN;
        if any(isfinite(y2))
            y2n = y2 ./ max(y2(isfinite(y2)));
            plot(ax, eta_eff2s, y2n, 'r--', 'LineWidth', 1.5, 'HandleVisibility','off');
        end

        xlabel(ax,'$\eta_{\mathrm{eff}}$ (V)','Interpreter','latex');

        if isfinite(etaStar)
            xline(ax, etaStar, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off');
        end

        text(ax, -0.25, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end

    % ---- tile 11 legend (Fig 2) ----
    axL2 = nexttile(tlo2, 11);
    axis(axL2,'off'); hold(axL2,'on');

    h4 = plot(axL2, NaN, NaN, 'g-', 'LineWidth', 2);
    h5 = plot(axL2, NaN, NaN, 'm-', 'LineWidth', 2);
    h6 = plot(axL2, NaN, NaN, 'r--', 'LineWidth', 1.5);
    h7 = plot(axL2, NaN, NaN, 'k--', 'LineWidth', 2.0);

    lg2 = legend(axL2, [h4 h5 h6 h7], ...
        {'$I_{\rm store}$', '$|\mathrm{d}\beta_{w}/\mathrm{d}\eta_{\rm eff}|$', '$\Pi_{\rm dens}$ (norm.)', '$\eta_{\rm eff}^*$'}, ...
        'Interpreter','latex', 'Location','best');
    lg2.Box = 'off';
    lg2.FontSize = 12;

    filename = 'FigureS24';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);
end

% =====================================================================
% helper: chi from UNSATURATED directional fluxes
% =====================================================================
function chi = compute_chi_from_params_unsat(P, eta_eff, maxExp)
    eta_eff = eta_eff(:);
    eps0 = 1e-30;

    % UNSATURATED directional fluxes
    JoxA  = P.jstarA .* exp( clip(P.aA.*eta_eff, maxExp) );
    JredA = P.jstarA .* exp( clip(P.bA.*eta_eff, maxExp) );
    JoxB  = P.jstarB .* exp( clip(P.aB.*eta_eff, maxExp) );
    JredB = P.jstarB .* exp( clip(P.bB.*eta_eff, maxExp) );

    Jox  = JoxA + JoxB;
    Jred = JredA + JredB;

    chi = (Jred - Jox) ./ (Jred + Jox + eps0);
    chi = max(min(chi, 1), -1);   % numerical safety
    chi(~isfinite(chi)) = 0;
end

%% =====================================================================
% Fitting Functions
%% =====================================================================

function P1 = fit_stage1_noohm(eta, jObs, Iref, f, maxExp, st1)
    % Stage 1: Fitting without ohmic drop
    eta = eta(:); jObs = jObs(:);
    jabs = abs(jObs(isfinite(jObs)));
    jmax = max(jabs); if ~isfinite(jmax) || jmax<=0, jmax=1e-3; end
    
    % Use residual-driven initial value estimation
    useResidualSeed = true;
    if useResidualSeed
        try
            P1A = fit_stage1_single_noohm(eta, jObs, Iref, f, maxExp, st1, 'A', []);
            P1B = fit_stage1_single_noohm(eta, jObs, Iref, f, maxExp, st1, 'B', []);
            
            jA0 = one_mode_explicit_noohm(P1A, eta, maxExp, 'A');
            jB0 = one_mode_explicit_noohm(P1B, eta, maxExp, 'B');
            
            rA = asinh(jObs./Iref) - asinh(jA0./Iref);
            rB = asinh(jObs./Iref) - asinh(jB0./Iref);
            rssA = sum(rA(isfinite(rA)).^2);
            rssB = sum(rB(isfinite(rB)).^2);
            
            if rssA <= rssB
                baseMode = 'A';
                jbase = jA0;
                Pbase = P1A;
            else
                baseMode = 'B';
                jbase = jB0;
                Pbase = P1B;
            end
            
            dj = (jObs - jbase);
            [jstar2, xi2, lam2, jlim2] = seed_second_mode_from_residual(eta, dj, f, maxExp, jmax);
            
            if baseMode == 'A'
                xiA0 = Pbase.xiA; lamA0 = Pbase.lamA; jstarA0 = Pbase.jstarA; jlimA0 = Pbase.jlimA;
                xiB0 = xi2; lamB0 = lam2; jstarB0 = jstar2; jlimB0 = jlim2;
            else
                xiB0 = Pbase.xiB; lamB0 = Pbase.lamB; jstarB0 = Pbase.jstarB; jlimB0 = Pbase.jlimB;
                xiA0 = xi2; lamA0 = lam2; jstarA0 = jstar2; jlimA0 = jlim2;
            end
            
        catch
            useResidualSeed = false;
        end
    end
    
    if ~useResidualSeed
        % Backup heuristic initial values
        xiA0 = 0.25; lamA0 = 0.25; jstarA0 = max(0.08*median(jabs(jabs>0)), 1e-10);
        xiB0 = 0.65; lamB0 = 0.18; jstarB0 = max(0.03*median(jabs(jabs>0)), 1e-10);
        jlim0 = max(prctile_fallback(jabs, 85), 1e-3);
        jlimA0 = 1.4*jlim0;
        jlimB0 = 1.2*jlim0;
    end
    
    p0 = [log(jstarA0), inv_sigmoid(xiA0), inv_softplus(lamA0), ...
          log(jstarB0), inv_sigmoid(xiB0), inv_softplus(lamB0), ...
          inv_softplus(jlimA0), inv_softplus(jlimB0)];
    
    lb = [-30, -12, -15, -30, -12, -15, ...
          inv_softplus(max(0.2*jmax,1e-6)), inv_softplus(max(0.2*jmax,1e-6))];
    ub = [30, 12, 15, 30, 12, 15, ...
          inv_softplus(max(30*jmax,1e-3)), inv_softplus(max(30*jmax,1e-3))];
    
    % Weak regularization
    w_xi = 0.02; L0 = 2.3; sx = 4.0;
    w_lim = 0.02; ratio_max = 8.0; sl = 1.5;
    
    resfun = @(p) residual_stage1_explicit(p, eta, jObs, Iref, f, maxExp, ...
        w_xi, L0, sx, w_lim, ratio_max, sl, jlimA0, jlimB0);
    
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxIterations', st1.maxIter, 'MaxFunctionEvaluations', st1.maxFeval, ...
        'TolFun', st1.tolFun, 'TolX', st1.tolX);
    
    best_p = p0; best_cost = inf;
    
    for s = 1:st1.nStarts
        if s == 1
            pS = p0;
        else
            pS = clamp_vec(p0 + 0.20*randn(size(p0)), lb, ub);
        end
        [p_try, resnorm, ~, exitflag] = lsqnonlin(resfun, pS, lb, ub, opts);
        if exitflag > 0 && resnorm < best_cost
            best_cost = resnorm; best_p = p_try;
        end
    end
    
    P1 = decode_stage1(best_p, f);
end

function r = residual_stage1_explicit(p, eta, jObs, Iref, f, maxExp, ...
    w_xi, L0, sx, w_lim, ratio_max, sl, jlimA0, jlimB0)
    
    P = decode_stage1(p, f);
    jpred = two_mode_explicit_noohm(P, eta, maxExp);
    
    r_data = asinh(jObs./Iref) - asinh(jpred./Iref);
    r_data(~isfinite(r_data)) = 1e6;
    
    excess_xi = max(0, abs(p(2))-L0) + max(0, abs(p(5))-L0);
    r_xi = sqrt(w_xi) * (excess_xi/sx);
    
    jlimA = softplus(p(7));
    jlimB = softplus(p(8));
    excess_lim = max(0, log(jlimA/(ratio_max*jlimA0))) + max(0, log(jlimB/(ratio_max*jlimB0)));
    r_lim = sqrt(w_lim) * (excess_lim/sl);
    
    r = [r_data; r_xi; r_lim];
end

function P = decode_stage1(p, f)
    P = struct();
    P.jstarA = exp(p(1));
    P.xiA = min(max(sigmoid(p(2)), 1e-6), 1-1e-6);
    P.lamA = softplus(p(3));
    
    P.jstarB = exp(p(4));
    P.xiB = min(max(sigmoid(p(5)), 1e-6), 1-1e-6);
    P.lamB = softplus(p(6));
    
    P.jlimA = max(softplus(p(7)), 1e-10);
    P.jlimB = max(softplus(p(8)), 1e-10);
    P.Rohm = 0;
    
    % Calculate a,b parameters
    P.aA = P.xiA * P.lamA * f;
    P.bA = -(1 - P.xiA) * P.lamA * f;
    P.aB = P.xiB * P.lamB * f;
    P.bB = -(1 - P.xiB) * P.lamB * f;
    
    % Ensure mode A is primary (jstarA >= jstarB)
    if P.jstarA < P.jstarB
        % Swap all parameters
        [P.jstarA, P.jstarB] = deal(P.jstarB, P.jstarA);
        [P.xiA, P.xiB] = deal(P.xiB, P.xiA);
        [P.lamA, P.lamB] = deal(P.lamB, P.lamA);
        [P.jlimA, P.jlimB] = deal(P.jlimB, P.jlimA);
        [P.aA, P.aB] = deal(P.aB, P.aA);
        [P.bA, P.bB] = deal(P.bB, P.bA);
        
        % Record swap flag
        P.modesSwapped = true;
    else
        P.modesSwapped = false;
    end
end

function j = two_mode_explicit_noohm(P, eta, maxExp)
    eta = eta(:);
    
    % Mode A
    ApA = clip(P.aA.*eta, maxExp);
    AmA = clip(P.bA.*eta, maxExp);
    xA = P.jstarA .* (exp(ApA) - exp(AmA));
    jA = wg_mass_energy_map_optimized(xA, P.jlimA);
    
    % Mode B
    ApB = clip(P.aB.*eta, maxExp);
    AmB = clip(P.bB.*eta, maxExp);
    xB = P.jstarB .* (exp(ApB) - exp(AmB));
    jB = wg_mass_energy_map_optimized(xB, P.jlimB);
    
    j = jA + jB;
end

function P2 = fit_stage2_full(eta, jObs, Iref, f, maxExp, st2, P1)
    % Stage 2: Full implicit fitting
    eta = eta(:); jObs = jObs(:);
    jabs = abs(jObs(isfinite(jObs)));
    jmax = max(jabs); if ~isfinite(jmax) || jmax<=0, jmax=1e-3; end
    
    % Calculate Rohm upper bound (data-driven)
    eta95 = prctile_fallback(abs(eta), 95);
    j95 = prctile_fallback(abs(jObs), 95);
    Rcap = max(0.8*eta95/(j95/1000+1e-12), 1e-12);
    
    p0 = encode_stage2_from_P1(P1, Rcap);
    
    lb = [-30, -12, -15, -30, -12, -15, inv_softplus(0), ...
          inv_softplus(max(0.2*jmax,1e-6)), inv_softplus(max(0.2*jmax,1e-6))];
    ub = [30, 12, 15, 30, 12, 15, inv_softplus(Rcap), ...
          inv_softplus(max(30*jmax,1e-3)), inv_softplus(max(30*jmax,1e-3))];
    
    j_init0 = two_mode_explicit_noohm(P1, eta, maxExp);
    
    w_xi = 0.02; L0 = 2.3; sx = 4.0;
    w_lim = 0.02; ratio_max = 8.0; sl = 1.5;
    jlimA0 = P1.jlimA; jlimB0 = P1.jlimB;
    
    resfun = make_residual_stage2(eta, jObs, Iref, f, maxExp, st2, ...
        w_xi, L0, sx, w_lim, ratio_max, sl, jlimA0, jlimB0, j_init0, Rcap);
    
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxIterations', st2.lsqMaxIter, 'MaxFunctionEvaluations', st2.lsqMaxFeval, ...
        'TolFun', st2.tolFun, 'TolX', st2.tolX);
    
    best_p = p0; best_cost = inf;
    
    for s = 1:st2.nStarts
        if s == 1
            pS = p0;
        else
            pS = clamp_vec(p0 + 0.16*randn(size(p0)), lb, ub);
        end
        [p_try, resnorm, ~, exitflag] = lsqnonlin(resfun, pS, lb, ub, opts);
        if exitflag > 0 && resnorm < best_cost
            best_cost = resnorm; best_p = p_try;
        end
    end
    
    P2 = decode_stage2(best_p, f, Rcap);
end

function p0 = encode_stage2_from_P1(P1, Rcap)
    Rohm0 = min(max(0.1, 0), Rcap);
    p0 = [log(P1.jstarA), inv_sigmoid(P1.xiA), inv_softplus(P1.lamA), ...
          log(P1.jstarB), inv_sigmoid(P1.xiB), inv_softplus(P1.lamB), ...
          inv_softplus(Rohm0), inv_softplus(P1.jlimA), inv_softplus(P1.jlimB)];
end

function resfun = make_residual_stage2(etaGrid, jObs, Iref, f, maxExp, st2, ...
    w_xi, L0, sx, w_lim, ratio_max, sl, jlimA0, jlimB0, j_init0, Rcap)
    
    p_last = [];
    j_last = [];
    
    resfun = @eval_res;
    
    function r = eval_res(p)
        P = decode_stage2(p, f, Rcap);
        
        if ~isempty(p_last) && norm(p-p_last) < 0.22 && numel(j_last) == numel(etaGrid)
            j0 = j_last;
        else
            j0 = j_init0;
        end
        
        jpred = safe_two_mode_model(P, etaGrid, maxExp, st2.newtonMaxIter, st2.newtonTol, j0);
        
        p_last = p;
        j_last = jpred;
        
        r_data = asinh(jObs./Iref) - asinh(jpred./Iref);
        r_data(~isfinite(r_data)) = 1e6;
        
        excess_xi = max(0, abs(p(2))-L0) + max(0, abs(p(5))-L0);
        r_xi = sqrt(w_xi) * (excess_xi/sx);
        
        jlimA = softplus(p(8));
        jlimB = softplus(p(9));
        excess_lim = max(0, log(jlimA/(ratio_max*jlimA0))) + max(0, log(jlimB/(ratio_max*jlimB0)));
        r_lim = sqrt(w_lim) * (excess_lim/sl);
        
        % Rohm weak prior
        r_R = [];
        if isfield(st2, 'wRohm') && st2.wRohm > 0
            R = softplus(p(7));
            r_R = sqrt(st2.wRohm) * ((R - st2.Rohm0) / (st2.RohmScale + 1e-12));
        end
        
        % Penalty to prevent j* approaching 0
        r_C = [];
        if isfield(st2, 'wCollapse') && st2.wCollapse > 0
            jA = exp(p(1));
            jB = exp(p(4));
            jf = st2.jstarFloor;
            r_C = sqrt(st2.wCollapse) * max(0, log(jf/(min(jA,jB)+1e-300)));
        end
        
        r = [r_data; r_xi; r_lim; r_R; r_C];
    end
end

function P = decode_stage2(p, f, Rcap)
    P = struct();
    
    P.jstarA = exp(p(1));
    P.xiA = min(max(sigmoid(p(2)), 1e-6), 1-1e-6);
    P.lamA = softplus(p(3));
    
    P.jstarB = exp(p(4));
    P.xiB = min(max(sigmoid(p(5)), 1e-6), 1-1e-6);
    P.lamB = softplus(p(6));
    
    P.Rohm = min(max(softplus(p(7)), 0), Rcap);
    P.jlimA = max(softplus(p(8)), 1e-10);
    P.jlimB = max(softplus(p(9)), 1e-10);
    
    % Calculate a,b parameters
    P.aA = P.xiA * P.lamA * f;
    P.bA = -(1 - P.xiA) * P.lamA * f;
    P.aB = P.xiB * P.lamB * f;
    P.bB = -(1 - P.xiB) * P.lamB * f;
    
    % Ensure mode A is primary (jstarA >= jstarB)
    if P.jstarA < P.jstarB
        % Swap all parameters
        [P.jstarA, P.jstarB] = deal(P.jstarB, P.jstarA);
        [P.xiA, P.xiB] = deal(P.xiB, P.xiA);
        [P.lamA, P.lamB] = deal(P.lamB, P.lamA);
        [P.jlimA, P.jlimB] = deal(P.jlimB, P.jlimA);
        [P.aA, P.aB] = deal(P.aB, P.aA);
        [P.bA, P.bB] = deal(P.bB, P.bA);
        
        % Record swap flag
        P.modesSwapped = true;
    else
        P.modesSwapped = false;
    end
    
    P.Rcap = Rcap;
end

%% =====================================================================
% Model Evaluation Functions
%% =====================================================================

function jpred = safe_two_mode_model(P, eta, maxExp, newtonMaxIter, newtonTol, j0)
    % Enhanced error-handling model evaluation function
    try
        jpred = two_mode_model_vecfast(P, eta, maxExp, newtonMaxIter, newtonTol, j0);
        
        % Check result validity
        if any(~isfinite(jpred))
            warning('Model returned non-finite values, using backup method');
            jpred = two_mode_explicit_noohm(P, eta, maxExp);
        end
        
        % Check numerical stability
        if max(abs(jpred)) > 1e10
            warning('Model output abnormally large, clipping');
            jpred = sign(jpred) .* min(abs(jpred), 1e10);
        end
        
    catch ME
        warning('Model evaluation failed: %s', ME.message);
        jpred = zeros(size(eta));
    end
end

function j = two_mode_model_vecfast(P, eta, maxExp, maxIter, tol, j0)
    % Vectorized Newton method for solving implicit equation
    eta = eta(:);
    n = numel(eta);
    
    if nargin < 6 || isempty(j0) || numel(j0) ~= n
        j = zeros(n,1);
    else
        j = j0(:);
        j(~isfinite(j)) = 0;
    end
    
    for it = 1:maxIter
        [Fv, dF] = eval_F_dF(P, eta, j, maxExp);
        
        rel = max(abs(Fv) ./ (1+abs(j)));
        if rel <= tol, return; end
        
        step = Fv ./ dF;
        step(~isfinite(step)) = 0;
        step = sign(step) .* min(abs(step), 0.8*(1+abs(j)));
        
        j1 = j - step;
        
        if it >= 3
            F0 = abs(Fv);
            [F1, ~] = eval_F_dF(P, eta, j1, maxExp);
            bad = ~(isfinite(F1) & abs(F1) < 0.90*F0);
            if any(bad)
                j2 = j - 0.5*step;
                [F2, ~] = eval_F_dF(P, eta, j2, maxExp);
                use2 = bad & isfinite(F2) & (abs(F2) < abs(F1));
                j1(use2) = j2(use2);
            end
        end
        
        j = j1;
    end
end

function [Fv, dF] = eval_F_dF(P, eta, j, maxExp)
    % Calculate implicit equation and its derivative
    eta_eff = eta - P.Rohm .* (j./1000);
    
    % Mode A
    ApA = clip(P.aA .* eta_eff, maxExp);
    AmA = clip(P.bA .* eta_eff, maxExp);
    ePA = exp(ApA); eMA = exp(AmA);
    
    xA = P.jstarA .* (ePA - eMA);
    dcoreA = P.jstarA .* (P.aA.*ePA - P.bA.*eMA);
    
    [jA, dsatA] = wg_mass_energy_map_optimized(xA, P.jlimA);
    
    % Mode B
    ApB = clip(P.aB .* eta_eff, maxExp);
    AmB = clip(P.bB .* eta_eff, maxExp);
    ePB = exp(ApB); eMB = exp(AmB);
    
    xB = P.jstarB .* (ePB - eMB);
    dcoreB = P.jstarB .* (P.aB.*ePB - P.bB.*eMB);
    
    [jB, dsatB] = wg_mass_energy_map_optimized(xB, P.jlimB);
    
    RHS = jA + jB;
    Fv = j - RHS;
    
    dRHS_detaeff = dsatA .* dcoreA + dsatB .* dcoreB;
    dRHS_dj = dRHS_detaeff .* (-P.Rohm./1000);
    dF = 1 - dRHS_dj;
    
    dF(~isfinite(dF) | abs(dF) < 1e-18) = 1;
    Fv(~isfinite(Fv)) = 0;
end

function [jA, jB] = two_mode_decompose_given_total(P, eta, jTot, maxExp)
    % Decompose based on mode swap status
    eta = eta(:); jTot = jTot(:);
    eta_eff = eta - P.Rohm .* (jTot./1000);
    
    % Check if mode swap flag exists
    if isfield(P, 'modesSwapped') && P.modesSwapped
        fprintf('  Note: Using swapped mode parameters for decomposition\n');
    end
    
    % Mode A calculation
    ApA = clip(P.aA .* eta_eff, maxExp);
    AmA = clip(P.bA .* eta_eff, maxExp);
    xA = P.jstarA .* (exp(ApA) - exp(AmA));
    jA = wg_mass_energy_map_optimized(xA, P.jlimA);
    
    % Mode B calculation
    ApB = clip(P.aB .* eta_eff, maxExp);
    AmB = clip(P.bB .* eta_eff, maxExp);
    xB = P.jstarB .* (exp(ApB) - exp(AmB));
    jB = wg_mass_energy_map_optimized(xB, P.jlimB);
end

%% =====================================================================
% Single Mode Fitting Functions
%% =====================================================================

function P1 = fit_stage1_single_noohm(eta, jObs, Iref, f, maxExp, st1, whichMode, Pseed)
    eta = eta(:); jObs = jObs(:);
    jabs = abs(jObs(isfinite(jObs)));
    jmax = max(jabs); if ~isfinite(jmax) || jmax<=0, jmax=1e-3; end
    
    if nargin >= 8 && isstruct(Pseed)
        jstar0 = max(get_jstar_mode(Pseed, whichMode), 1e-12);
        xi0 = get_xi_mode(Pseed, whichMode);
        lam0 = get_lam_mode(Pseed, whichMode);
        jlim0 = max(get_jlim_mode(Pseed, whichMode), 1e-6);
    else
        xi0 = 0.25; lam0 = 0.25;
        jstar0 = max(0.08*median(jabs(jabs>0)), 1e-10);
        jlim0 = max(prctile_fallback(jabs,85), 1e-3);
    end
    
    p0 = [log(jstar0), inv_sigmoid(xi0), inv_softplus(lam0), inv_softplus(jlim0)];
    lb = [-30, -12, -15, inv_softplus(max(0.2*jmax,1e-6))];
    ub = [30, 12, 15, inv_softplus(max(30*jmax,1e-3))];
    
    resfun = @(p) residual_stage1_single(p, eta, jObs, Iref, f, maxExp, whichMode);
    
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxIterations', st1.maxIter, 'MaxFunctionEvaluations', st1.maxFeval, ...
        'TolFun', st1.tolFun, 'TolX', st1.tolX);
    
    best_p = p0; best_cost = inf;
    
    for s = 1:st1.nStarts
        if s == 1
            pS = p0;
        else
            pS = clamp_vec(p0 + 0.20*randn(size(p0)), lb, ub);
        end
        [p_try, resnorm, ~, exitflag] = lsqnonlin(resfun, pS, lb, ub, opts);
        if exitflag > 0 && resnorm < best_cost
            best_cost = resnorm; best_p = p_try;
        end
    end
    
    P1 = decode_stage1_single(best_p, f, whichMode);
end

function r = residual_stage1_single(p, eta, jObs, Iref, f, maxExp, whichMode)
    P = decode_stage1_single(p, f, whichMode);
    jpred = one_mode_explicit_noohm(P, eta, maxExp, whichMode);
    r = asinh(jObs./Iref) - asinh(jpred./Iref);
    r(~isfinite(r)) = 1e6;
end

function P = decode_stage1_single(p, f, whichMode)
    P = struct();
    jstar = exp(p(1));
    xi = min(max(sigmoid(p(2)), 1e-6), 1-1e-6);
    lam = softplus(p(3));
    jlim = max(softplus(p(4)), 1e-10);
    
    xi_off = 0.5; lam_off = 0.2; jlim_off = 1; jstar_off = 0;
    
    if whichMode == 'A'
        P.jstarA = jstar; P.xiA = xi; P.lamA = lam; P.jlimA = jlim;
        P.jstarB = jstar_off; P.xiB = xi_off; P.lamB = lam_off; P.jlimB = jlim_off;
    else
        P.jstarB = jstar; P.xiB = xi; P.lamB = lam; P.jlimB = jlim;
        P.jstarA = jstar_off; P.xiA = xi_off; P.lamA = lam_off; P.jlimA = jlim_off;
    end
    P.Rohm = 0;
    
    P.aA = P.xiA * P.lamA * f;
    P.bA = -(1 - P.xiA) * P.lamA * f;
    P.aB = P.xiB * P.lamB * f;
    P.bB = -(1 - P.xiB) * P.lamB * f;
end

function j = one_mode_explicit_noohm(P, eta, maxExp, whichMode)
    eta = eta(:);
    if whichMode == 'A'
        Ap = clip(P.aA.*eta, maxExp);
        Am = clip(P.bA.*eta, maxExp);
        x = P.jstarA .* (exp(Ap) - exp(Am));
        j = wg_mass_energy_map_optimized(x, P.jlimA);
    else
        Ap = clip(P.aB.*eta, maxExp);
        Am = clip(P.bB.*eta, maxExp);
        x = P.jstarB .* (exp(Ap) - exp(Am));
        j = wg_mass_energy_map_optimized(x, P.jlimB);
    end
end

function P2 = fit_stage2_single_full(eta, jObs, Iref, f, maxExp, st2, P1, whichMode)
    eta = eta(:); jObs = jObs(:);
    jabs = abs(jObs(isfinite(jObs)));
    jmax = max(jabs); if ~isfinite(jmax) || jmax<=0, jmax=1e-3; end
    
    eta95 = prctile_fallback(abs(eta), 95);
    j95 = prctile_fallback(abs(jObs), 95);
    Rcap = max(0.8*eta95/(j95/1000 + 1e-12), 1e-12);
    
    jstar0 = get_jstar_mode(P1, whichMode);
    xi0 = get_xi_mode(P1, whichMode);
    lam0 = get_lam_mode(P1, whichMode);
    jlim0 = get_jlim_mode(P1, whichMode);
    Rohm0 = min(max(0.1,0), Rcap);
    
    p0 = [log(jstar0), inv_sigmoid(xi0), inv_softplus(lam0), ...
          inv_softplus(Rohm0), inv_softplus(jlim0)];
    
    lb = [-30, -12, -15, inv_softplus(0), inv_softplus(max(0.2*jmax,1e-6))];
    ub = [30, 12, 15, inv_softplus(Rcap), inv_softplus(max(30*jmax,1e-3))];
    
    j_init0 = one_mode_explicit_noohm(P1, eta, maxExp, whichMode);
    
    resfun = make_residual_stage2_single(eta, jObs, Iref, f, maxExp, st2, whichMode, j_init0, Rcap);
    
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxIterations', st2.lsqMaxIter, 'MaxFunctionEvaluations', st2.lsqMaxFeval, ...
        'TolFun', st2.tolFun, 'TolX', st2.tolX);
    
    best_p = p0; best_cost = inf;
    
    for s = 1:st2.nStarts
        if s == 1
            pS = p0;
        else
            pS = clamp_vec(p0 + 0.16*randn(size(p0)), lb, ub);
        end
        [p_try, resnorm, ~, exitflag] = lsqnonlin(resfun, pS, lb, ub, opts);
        if exitflag > 0 && resnorm < best_cost
            best_cost = resnorm; best_p = p_try;
        end
    end
    
    P2 = decode_stage2_single(best_p, f, Rcap, whichMode);
end

function resfun = make_residual_stage2_single(etaGrid, jObs, Iref, f, maxExp, st2, whichMode, j_init0, Rcap)
    p_last = [];
    j_last = [];
    
    resfun = @eval_res;
    
    function r = eval_res(p)
        P = decode_stage2_single(p, f, Rcap, whichMode);
        
        if ~isempty(p_last) && norm(p-p_last) < 0.22 && numel(j_last) == numel(etaGrid)
            j0 = j_last;
        else
            j0 = j_init0;
        end
        
        jpred = safe_two_mode_model(P, etaGrid, maxExp, st2.newtonMaxIter, st2.newtonTol, j0);
        
        p_last = p;
        j_last = jpred;
        
        r = asinh(jObs./Iref) - asinh(jpred./Iref);
        r(~isfinite(r)) = 1e6;
        
        if isfield(st2, 'wRohm') && st2.wRohm > 0
            R = softplus(p(4));
            rR = sqrt(st2.wRohm) * ((R - st2.Rohm0) / (st2.RohmScale + 1e-12));
            r = [r; rR];
        end
    end
end

function P = decode_stage2_single(p, f, Rcap, whichMode)
    P = struct();
    jstar = exp(p(1));
    xi = min(max(sigmoid(p(2)), 1e-6), 1-1e-6);
    lam = softplus(p(3));
    Rohm = min(max(softplus(p(4)), 0), Rcap);
    jlim = max(softplus(p(5)), 1e-10);
    
    xi_off = 0.5; lam_off = 0.2; jlim_off = 1; jstar_off = 0;
    
    if whichMode == 'A'
        P.jstarA = jstar; P.xiA = xi; P.lamA = lam; P.jlimA = jlim;
        P.jstarB = jstar_off; P.xiB = xi_off; P.lamB = lam_off; P.jlimB = jlim_off;
    else
        P.jstarB = jstar; P.xiB = xi; P.lamB = lam; P.jlimB = jlim;
        P.jstarA = jstar_off; P.xiA = xi_off; P.lamA = lam_off; P.jlimA = jlim_off;
    end
    P.Rohm = Rohm;
    
    P.aA = P.xiA * P.lamA * f;
    P.bA = -(1 - P.xiA) * P.lamA * f;
    P.aB = P.xiB * P.lamB * f;
    P.bB = -(1 - P.xiB) * P.lamB * f;
    
    P.Rcap = Rcap;
end

%% =====================================================================
% AIC and Diagnostic Functions
%% =====================================================================

function rss = rss_asinh_model(P, eta, jObs, Iref, maxExp, newtonMaxIter, newtonTol)
    jpred = safe_two_mode_model(P, eta(:), maxExp, newtonMaxIter, newtonTol, []);
    r = asinh(jObs(:)./Iref) - asinh(jpred(:)./Iref);
    r(~isfinite(r)) = 0;
    rss = sum(r.^2);
end

function AIC = aic_from_rss(rss, n, k)
    n = max(n,1);
    rss = max(rss, 1e-300);
    AIC = n*log(rss/n) + 2*k;
end

function [modeSel, info] = decide_single_mode_by_contrib(P, etaData, maxExp, newtonMaxIter, newtonTol, prune)
    etaData = etaData(:);
    lo = prctile_fallback(etaData, prune.etaPctRange(1));
    hi = prctile_fallback(etaData, prune.etaPctRange(2));
    if ~(isfinite(lo) && isfinite(hi) && hi>lo)
        lo = min(etaData);
        hi = max(etaData);
    end
    etaGrid = linspace(lo, hi, prune.nGrid).';
    
    jTot = safe_two_mode_model(P, etaGrid, maxExp, newtonMaxIter, newtonTol, []);
    [jA, jB] = two_mode_decompose_given_total(P, etaGrid, jTot, maxExp);
    
    jt = abs(jTot);
    jmax = max(jt(isfinite(jt)));
    if ~(isfinite(jmax) && jmax>0)
        modeSel = 'AB';
        info = struct('fracA', nan, 'fracB', nan, 'peakA', nan, 'peakB', nan, 'N', 0);
        return;
    end
    
    mask = isfinite(etaGrid) & isfinite(jTot) & isfinite(jA) & isfinite(jB) & (jt >= prune.jRelMin*jmax);
    if prune.cathodicOnly
        mask = mask & (etaGrid < 0);
    end
    
    N = nnz(mask);
    if N < 20
        modeSel = 'AB';
        info = struct('fracA', nan, 'fracB', nan, 'peakA', nan, 'peakB', nan, 'N', N);
        return;
    end
    
    x = etaGrid(mask);
    JA = abs(jA(mask));
    JB = abs(jB(mask));
    JT = abs(jTot(mask));
    
    areaT = trapz(x, JT);
    areaA = trapz(x, JA);
    areaB = trapz(x, JB);
    
    fracA = areaA/(areaT + 1e-300);
    fracB = areaB/(areaT + 1e-300);
    peakA = max(JA)/(max(JT) + 1e-300);
    peakB = max(JB)/(max(JT) + 1e-300);
    
    info = struct('fracA', fracA, 'fracB', fracB, 'peakA', peakA, 'peakB', peakB, 'N', N);
    
    weakA = (fracA < prune.fracTol) && (peakA < prune.peakTol);
    weakB = (fracB < prune.fracTol) && (peakB < prune.peakTol);
    
    if weakA && ~weakB
        modeSel = 'B';
    elseif weakB && ~weakA
        modeSel = 'A';
    else
        modeSel = 'AB';
    end
end

%% =====================================================================
% Tafel Fitting Helper Functions
%% =====================================================================

function F = fit_tafel_relative_threshold(x, j, cfg)
    F = struct('ok', false, 'beta_mVdec', nan, 'R2', nan, 'N', 0, ...
        'x1', nan, 'x2', nan, 'xLine', [], 'yLine', []);
    x = x(:); j = j(:);
    
    if strcmpi(cfg.branch, 'cathodic')
        mask = (j < 0);
    elseif strcmpi(cfg.branch, 'anodic')
        mask = (j > 0);
    else
        mask = true(size(j));
    end
    
    x = x(mask); j = j(mask);
    y = log10(abs(j));
    good = isfinite(x) & isfinite(y) & isfinite(j);
    x = x(good); j = j(good); y = y(good);
    
    if numel(x) < cfg.minN, return; end
    absj = abs(j);
    jmax = max(absj);
    if ~(isfinite(jmax) && jmax>0), return; end
    
    th_lo = cfg.f_lo * jmax;
    th_hi = cfg.f_hi * jmax;
    sel = (absj >= th_lo) & (absj <= th_hi);
    if nnz(sel) < cfg.minN, return; end
    
    xs = x(sel); ys = y(sel);
    if (max(ys) - min(ys)) < cfg.minDecades, return; end
    
    p = polyfit(xs, ys, 1);
    yhat = polyval(p, xs);
    SSr = sum((ys - yhat).^2);
    SSt = sum((ys - mean(ys)).^2);
    R2 = 1 - SSr/(SSt + 1e-12);
    if R2 < cfg.R2_min, return; end
    
    m = p(1);
    if ~isfinite(m) || abs(m) < 1e-12, return; end
    
    F.ok = true;
    F.beta_mVdec = 1000/abs(m);
    F.R2 = R2;
    F.N = numel(xs);
    F.x1 = min(xs);
    F.x2 = max(xs);
    
    xLine = linspace(F.x1, F.x2, 30).';
    yLine = polyval(p, xLine);
    F.xLine = xLine;
    F.yLine = yLine;
end

%% =====================================================================
% Mathematical and Utility Functions
%% =====================================================================

function y = sigmoid(x)
    y = 1./(1 + exp(-x));
end

function y = softplus(x)
    y = log1p(exp(-abs(x))) + max(x, 0);
end

function x = inv_softplus(y)
    x = log(max(exp(y) - 1, 1e-12));
end

function x = inv_sigmoid(y)
    y = min(max(y, 1e-12), 1-1e-12);
    x = log(y./(1-y));
end

function z = clip(x, maxExp)
    z = min(max(x, -maxExp), maxExp);
end

function v = clamp_vec(v, lo, hi)
    v = min(max(v, lo), hi);
end

function q = prctile_fallback(x, p)
    x = x(:); x = x(isfinite(x));
    if isempty(x), q = NaN; return; end
    x = sort(x); p = min(max(p, 0), 100);
    if numel(x) == 1, q = x; return; end
    r = 1 + (numel(x)-1)*(p/100);
    k = floor(r); a = r - k;
    k = max(1, min(k, numel(x)-1));
    q = (1-a)*x(k) + a*x(k+1);
end

function y = safe_log10abs(j)
    y = log10(abs(j));
    y(~isfinite(y)) = NaN;
end

function S = set_default(S, field, value)
    if ~isfield(S, field) || isempty(S.(field))
        S.(field) = value;
    end
end

function save_figure(fig_handle, filename, folder)
    if nargin < 3 || isempty(folder), folder = 'figures'; end
    folder = char(folder); filename = char(filename);
    
    if isfolder(folder)
        folderAbs = folder;
    else
        folderAbs = fullfile(pwd, folder);
    end
    
    [ok, msg] = mkdir(folderAbs);
    if ~(ok || isfolder(folderAbs))
        warning('Cannot create folder "%s": %s. Using temporary directory.', folderAbs, msg);
        folderAbs = fullfile(tempdir, 'figures');
        mkdir(folderAbs);
    end
    
    pngFile = fullfile(folderAbs, [filename, '.png']);
    figFile = fullfile(folderAbs, [filename, '.fig']);
    
    try
        exportgraphics(fig_handle, pngFile, 'Resolution', 300);
    catch
        print(fig_handle, pngFile, '-dpng', '-r300');
    end
    
    try
        savefig(fig_handle, figFile);
    catch
        saveas(fig_handle, figFile);
    end
    
    fprintf('Figures saved to:\n  %s\n  %s\n', pngFile, figFile);
end

%% =====================================================================
% Waveguide Mass-Energy Mapping Function (Optimized)
%% =====================================================================

function [j, djdx] = wg_mass_energy_map_optimized(x, jlim)
    % Optimized mass-energy mapping function
    jlim = max(jlim, 1e-12);
    z = x ./ jlim;
    
    % Use hybrid method for numerical stability
    small_mask = abs(z) < 0.1;
    large_mask = ~small_mask;
    
    j = zeros(size(x), 'like', x);
    
    % Small value approximation (Taylor expansion)
    if any(small_mask(:))
        z_small = z(small_mask);
        j(small_mask) = z_small .* jlim .* (1 - 0.5*z_small.^2);
    end
    
    % Large value exact calculation
    if any(large_mask(:))
        z_large = z(large_mask);
        den = sqrt(1 + z_large.^2);
        j(large_mask) = (z_large ./ den) .* jlim;
    end
    
    if nargout > 1
        djdx = zeros(size(x), 'like', x);
        if any(small_mask(:))
            djdx(small_mask) = 1 - 1.5*z_small.^2;
        end
        if any(large_mask(:))
            djdx(large_mask) = 1 ./ (den.^3);
        end
    end
end

%% =====================================================================
% Simple Get Functions
%% =====================================================================

function v = get_jstar_mode(P, which)
    if which == 'A'
        v = P.jstarA;
    else
        v = P.jstarB;
    end
end

function v = get_xi_mode(P, which)
    if which == 'A'
        v = P.xiA;
    else
        v = P.xiB;
    end
end

function v = get_lam_mode(P, which)
    if which == 'A'
        v = P.lamA;
    else
        v = P.lamB;
    end
end

function v = get_jlim_mode(P, which)
    if which == 'A'
        v = P.jlimA;
    else
        v = P.jlimB;
    end
end

function P = ensure_modeA_primary(P)
    % Ensure mode A jstar >= mode B jstar
    % If jstarA < jstarB, swap all parameters of the two modes
    
    if ~isfield(P, 'jstarA') || ~isfield(P, 'jstarB')
        return;  % If not dual mode, return directly
    end
    
    % Check if swap needed
    if P.jstarA >= P.jstarB
        return;  % Already mode A primary, no swap needed
    end
    
    fprintf('  Note: Swapping mode parameters to ensure mode A primary (jstarA=%.2e -> %.2e, jstarB=%.2e -> %.2e)\n', ...
        P.jstarA, P.jstarB, P.jstarB, P.jstarA);
    
    % Swap all parameters
    [P.jstarA, P.jstarB] = deal(P.jstarB, P.jstarA);
    [P.xiA, P.xiB] = deal(P.xiB, P.xiA);
    [P.lamA, P.lamB] = deal(P.lamB, P.lamA);
    [P.jlimA, P.jlimB] = deal(P.jlimB, P.jlimA);
    [P.aA, P.aB] = deal(P.aB, P.aA);
    [P.bA, P.bB] = deal(P.bB, P.bA);
    
    % Record swap history
    if ~isfield(P, 'swapHistory')
        P.swapHistory = struct();
    end
    P.swapHistory.modesSwapped = true;
    P.swapHistory.original_jstarA = P.jstarB;  % Note: after swap jstarA is original jstarB
    P.swapHistory.original_jstarB = P.jstarA;  % after swap jstarB is original jstarA
end

function plot_waveguide_electrochem_summary_figure(RES, maxExp, newtonMaxIter, newtonTol)
% Summary figure (2x2) with optimum defined by MAX useful-output density.
%   Pi_use(eta_eff)  = max(0,-jTot) * (1 - rho^2)
%   Pi_dens(eta_eff) = Pi_use / (|eta_eff| + eta0)
% Optimum in panel (c): argmax Pi_dens on cathodic side (eta_eff<0).

    if nargin < 2 || isempty(maxExp), maxExp = 70; end
    if nargin < 3 || isempty(newtonMaxIter), newtonMaxIter = 10; end
    if nargin < 4 || isempty(newtonTol), newtonTol = 5e-11; end

    nG = numel(RES);

    % === FIX: keep original 7 lines() colors + add more if needed ===
    C  = extended_line_colors(nG);

    % ---------- options ----------
    opt = struct();
    opt.useEtaEff = true;
    opt.nEta = 500;
    opt.etaPctRange = [5 95];
    opt.useFixedWindow = true;
    opt.fixedEtaWindow = [-1.5 1.0];
    opt.cathodicOnly = false;

    % useful-density regularization (V): prevents blow-up near 0
    opt.eta0 = 0.05;   % try 0.05~0.20

    % ---------- precompute waveguide diagnostics per pH ----------
    storage_data = cell(nG,1);
    for g = 1:nG
        eta_raw = RES(g).eta(:);

        if opt.useFixedWindow
            lo = opt.fixedEtaWindow(1);
            hi = opt.fixedEtaWindow(2);
        else
            lo = prctile_fallback(eta_raw, opt.etaPctRange(1));
            hi = prctile_fallback(eta_raw, opt.etaPctRange(2));
            if ~(isfinite(lo) && isfinite(hi) && hi > lo)
                lo = min(eta_raw); hi = max(eta_raw);
            end
        end

        etaGrid = linspace(lo, hi, opt.nEta).';

        [Gmag, Gph, ~, eta_eff_map, U, S_rel, beta_w] = ...
            calc_caismith_from_params_improved(RES(g).P, etaGrid, opt.useEtaEff, ...
                                               maxExp, newtonMaxIter, newtonTol, 'none');

        storage_data{g} = struct();
        storage_data{g}.etaGrid_raw = etaGrid;
        storage_data{g}.eta_eff_map = eta_eff_map(:);
        storage_data{g}.rho = Gmag(:);
        storage_data{g}.theta = Gph(:);
        storage_data{g}.U = U(:);
        storage_data{g}.S = S_rel(:);
        storage_data{g}.beta_w = beta_w(:);
        storage_data{g}.chi = -beta_w(:);
    end

    % ---------- figure ----------
    width_cm = 22.5;
    height_cm = 18;
    fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
    tlo = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

    % ==================== (a) polarization fits ====================
    ax1 = nexttile(tlo,1);
    hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

    all_j = vertcat(RES.j);
    yL = prctile_fallback(abs(all_j), 99);
    yL = max(yL, max(abs(all_j))*1.1);
    if ~(isfinite(yL) && yL>0), yL = 1; end

    for g = 1:nG
        eta = RES(g).eta(:);
        j   = RES(g).j(:);
        P   = RES(g).P;

        scatter(ax1, eta, j, 18, 'MarkerEdgeColor', C(g,:), ...
            'MarkerFaceColor', C(g,:), 'MarkerFaceAlpha', 0.45, ...
            'HandleVisibility','off');

        etaFine = linspace(min(eta), max(eta), 700).';
        jTot = safe_two_mode_model(P, etaFine, maxExp, newtonMaxIter, newtonTol, []);

        plot(ax1, etaFine, jTot, '-', 'Color', C(g,:), 'LineWidth', 2.0, ...
            'DisplayName', sprintf('pH=%.1f', RES(g).pH));
    end

    xlabel(ax1, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex');
    ylabel(ax1, '$j$ (mA cm$^{-2}$)', 'Interpreter','latex');
    yline(ax1, 0, 'k-', 'HandleVisibility','off');
    xline(ax1, 0, 'k-', 'HandleVisibility','off');
    ylim(ax1, [-yL, 0.15*yL]);
    if opt.useFixedWindow, xlim(ax1, opt.fixedEtaWindow); end
    lg1 = legend(ax1, 'Location','east', 'FontSize', 12);
    legend(ax1,'boxoff');

    % ==================== (b) rho and chi vs eta_eff ====================
    ax2 = nexttile(tlo,2);
    hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
    set(ax2,'FontSize',10,'Layer','top','Clipping','on');

    yyaxis(ax2,'left');
    for g = 1:nG
        d = storage_data{g};
        plot(ax2, d.eta_eff_map, d.rho, '-', 'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
    end
    ylabel(ax2, '$\rho=|\Gamma_g|$', 'Interpreter','latex', 'FontSize',12, 'Color','k');
    ylim(ax2, [0, 1.05]);
    ax2.YColor = 'k';

    xlabel(ax2, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex','FontSize',12);
    if opt.useFixedWindow, xlim(ax2, opt.fixedEtaWindow); end

    yyaxis(ax2,'right');
    for g = 1:nG
        d = storage_data{g};
        chi = max(min(d.chi,1),-1);
        plot(ax2, d.eta_eff_map, chi, '--', 'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
    end
    ylabel(ax2, '$\beta_{w}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.85 0.33 0.10]);
    ylim(ax2, [-1.05, 1.05]);
    ax2.YColor = [0.85 0.33 0.10];

    yyaxis(ax2,'left');  hR = plot(ax2, NaN, NaN, 'k-',  'LineWidth',2.0, 'DisplayName', '$\rho=|\Gamma_g|$');
    yyaxis(ax2,'right'); hC = plot(ax2, NaN, NaN, 'k--', 'LineWidth',2.0, 'DisplayName', '$\beta_{w}$');
    legend1=legend(ax2, [hR hC], 'Location','best', 'Interpreter','latex', 'FontSize',12, 'Box','off');
    set(legend1,...
        'Position',[0.621377893541425 0.717865915229878 0.114346409217984 0.0559640529108982],...
        'Interpreter','latex',...
        'FontSize',12);

    % ==================== (c) output density + optimum ====================
    ax3 = nexttile(tlo,3);
    hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');

    eta0 = opt.eta0;

    for g = 1:nG
        P = RES(g).P;
        d = storage_data{g};
        etaGrid = d.etaGrid_raw(:);

        jTot = safe_two_mode_model(P, etaGrid, maxExp, newtonMaxIter, newtonTol, []);
        eta_eff = etaGrid - P.Rohm.*(jTot./1000);
        rho = d.rho(:);

        n = min([numel(eta_eff), numel(rho), numel(jTot)]);
        eta_eff = eta_eff(1:n); rho = rho(1:n); jTot = jTot(1:n);

        Puse  = max(0, -jTot) .* (1 - rho.^2);
        Puse =abs(jTot) .* (1 - rho.^2);     % mA/cm^2
        Pdens = Puse ./ (abs(eta_eff) + eta0);

        if opt.cathodicOnly
            mask = (eta_eff < 0) & isfinite(eta_eff) & isfinite(Pdens);
        else
            mask = isfinite(eta_eff) & isfinite(Pdens);
        end
        if nnz(mask) < 8, continue; end

        x = eta_eff(mask);
        z = Pdens(mask);

        [x, ord] = sort(x);
        z = z(ord);

        plot(ax3, x, z, '-', 'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');

        [zmax, im] = max(z);
        xstar = x(im);

        plot(ax3, xstar, zmax, 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', C(g,:), 'MarkerEdgeColor','k', 'LineWidth', 1.0, ...
            'DisplayName', sprintf('$\\eta^*_{\\mathrm{eff}}=%.2f$ V', xstar));

        fprintf('pH=%.2f: density optimum (eta0=%.2f): eta_eff*=%.3f V, Pi_dens*=%.3g\n', ...
            RES(g).pH, eta0, xstar, zmax);
    end

    set(ax3,'YScale','log');
    ax3.YAxis.MinorTick = 'off';
    set(ax3, 'MinorGridLineStyle', 'none');
    xlabel(ax3, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex');
    ylabel(ax3, '$\Pi_{\mathrm{dens}}$ (mA cm$^{-2}$ V$^{-1}$)', 'Interpreter','latex');
    if opt.useFixedWindow, xlim(ax3, opt.fixedEtaWindow); end
    lg3 = legend(ax3, 'Location','southwest', 'FontSize', 12, 'Interpreter','latex');
    set(lg3,'NumColumns',2);
    legend(ax3,'boxoff');
    set(lg3,...
        'Position',[0.0858296896907499 0.0859403890721938 0.299076816334444 0.158660140692019],...
        'NumColumns',2,...
        'Interpreter','latex');

    % ==================== (d) normalized U and |S| ====================
    ax4 = nexttile(tlo,4);
    hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');

    for g = 1:nG
        d = storage_data{g};
        x = d.eta_eff_map;
        yU = max(d.U, 0);
        yS = max(abs(d.S), 0);
        plot(ax4, x, yU./max(yU), '-',  'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
        plot(ax4, x, yS./max(yS), '--', 'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
    end

    xlabel(ax4, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex');
    ylabel(ax4, 'Normalized $\mathcal{U}$ and $|\mathcal{S}|$', 'Interpreter','latex');
    if opt.useFixedWindow, xlim(ax4, opt.fixedEtaWindow); end

    hU = plot(ax4, NaN, NaN, 'k-',  'LineWidth', 2.0, 'DisplayName', '$\mathcal{U}$');
    hS = plot(ax4, NaN, NaN, 'k--', 'LineWidth', 2.0, 'DisplayName', '$|\mathcal{S}|$');
    lg4 = legend(ax4, [hU hS], 'Location','southwest', 'FontSize', 12, 'Interpreter','latex');
    legend(ax4,'boxoff');
    set(ax4,'YScale','log');

    % ---- panel labels ----
    text(ax1, -0.12, 0.98, 'a', 'Units','normalized', 'FontSize',16, 'FontWeight','bold', ...
        'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');
    text(ax2, -0.12, 0.98, 'b', 'Units','normalized', 'FontSize',16, 'FontWeight','bold', ...
        'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');
    text(ax3, -0.12, 0.98, 'c', 'Units','normalized', 'FontSize',16, 'FontWeight','bold', ...
        'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');
    text(ax4, -0.12, 0.98, 'd', 'Units','normalized', 'FontSize',16, 'FontWeight','bold', ...
        'FontName','Helvetica', 'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top');

    % ---- save ----
    filename = 'Figure19';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);
end

% ======================================================================
% Helper: keep lines(7) exactly, then add extra distinct colors if needed
% ======================================================================
function C = extended_line_colors(n)
    base = lines(7);
    if n <= 7
        C = base(1:n,:);
        return;
    end

    C = zeros(n,3);
    C(1:7,:) = base;

    k = n - 7;

    % Candidate pool in HSV (dense), then pick farthest colors from base (greedy max-min)
    Nh = 720;
    H = linspace(0,1,Nh+1)'; H(end)=[];
    S = 0.80*ones(size(H));
    V = 0.90*ones(size(H));
    cand = hsv2rgb([H S V]);

    % Remove candidates too close to white (optional; keeps print-friendly)
    cand = cand(sum((cand - 1).^2,2) > 0.15);

    minD2 = inf(size(cand,1),1);
    for j = 1:7
        d2 = sum((cand - base(j,:)).^2,2);
        minD2 = min(minD2, d2);
    end

    picked = zeros(k,3);
    for i = 1:k
        [~, idx] = max(minD2);
        picked(i,:) = cand(idx,:);

        % update min-distance considering new pick
        d2new = sum((cand - picked(i,:)).^2,2);
        minD2 = min(minD2, d2new);
        minD2(idx) = -inf; % don't pick again
    end

    C(8:end,:) = picked;
end

function OUT = plot_caismith_electrochem_waveguide(RES, opt)
% plot_caismith_electrochem_waveguide (FINAL)
% ------------------------------------------------------------
% Publication-style Cai–Smith disks for electrochemical waveguide mapping.
% Key features (FINAL):
%   (1) Color = normalized output density per pH:
%         Cnorm = Pi_dens / max(Pi_dens)  within each pH (on mask)
%       so Cnorm ∈ [0,1] for all pH.
%   (2) ONE global colorbar on the far right (no per-tile colorbars).
%   (3) No bar chart in tile 8; tile 8 left empty.
%   (4) Mark and annotate the Pi_dens maximum point on each disk
%       (eta_eff* and Pi_dens^max).

if nargin < 2, opt = struct(); end
opt = setdef(opt,'useEtaEff',true);
opt = setdef(opt,'nEta',700);
opt = setdef(opt,'etaPctRange',[5 95]);
opt = setdef(opt,'useFixedWindow',true);
opt = setdef(opt,'fixedEtaWindow',[-1.5 1.0]);
opt = setdef(opt,'cathodicOnly',true);
opt = setdef(opt,'eta0',0.05);
opt = setdef(opt,'zoomToData',true);
opt = setdef(opt,'zoomMargin',0.12);
opt = setdef(opt,'showUnitCircle',true);
opt = setdef(opt,'showGrid',false);
opt = setdef(opt,'maxExp',70);
opt = setdef(opt,'newtonMaxIter',10);
opt = setdef(opt,'newtonTol',5e-11);
opt = setdef(opt,'pointSize',14);
opt = setdef(opt,'savePng','');
opt = setdef(opt,'colormapName','turbo');

assert(isstruct(RES) && ~isempty(RES), 'RES must be a non-empty struct array.');
assert(isfinite(opt.eta0) && opt.eta0 > 0, 'opt.eta0 must be > 0 for stable Pi_dens normalization.');

% ---- global style (publication-ish) ----
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',12);
set(groot,'DefaultAxesLineWidth',0.9);
set(groot,'DefaultLineLineWidth',1.6);

nG = numel(RES);
OUT = repmat(struct('pH',NaN,'eta',[],'eta_eff',[],'jTot',[], ...
    'Jox',[],'Jred',[],'U',[],'S',[],'beta',[],'I_store',[], ...
    'rho',[],'theta',[],'Gamma',[],'Pi_use',[],'Pi_dens',[], ...
    'C',[],'mask',[],'clim',[0 1], 'idxStar',NaN, 'Pi_dens_max',NaN), nG, 1);

% ---------------- compute per pH ----------------
for g = 1:nG
    P = RES(g).P;
    P = ensure_ab_fields(P);  % ensure aA,bA,aB,bB exist

    % eta grid
    if opt.useFixedWindow
        lo = opt.fixedEtaWindow(1); hi = opt.fixedEtaWindow(2);
    else
        eta_raw = RES(g).eta(:);
        lo = prctile_fallback(eta_raw, opt.etaPctRange(1));
        hi = prctile_fallback(eta_raw, opt.etaPctRange(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw); hi = max(eta_raw);
        end
    end
    etaGrid = linspace(lo, hi, opt.nEta).';

    % implicit j(eta) -> eta_eff
    jTot = safe_two_mode_model(P, etaGrid, opt.maxExp, opt.newtonMaxIter, opt.newtonTol, []);
    if opt.useEtaEff
        eta_eff = etaGrid - P.Rohm.*(jTot./1000); % mA/cm2 -> A/cm2 for Rohm
    else
        eta_eff = etaGrid;
    end

    % UNSATURATED directional fluxes (used for Gamma only)
    JoxA  = P.jstarA .* exp( clip(P.aA.*eta_eff, opt.maxExp) );
    JredA = P.jstarA .* exp( clip(P.bA.*eta_eff, opt.maxExp) );
    JoxB  = P.jstarB .* exp( clip(P.aB.*eta_eff, opt.maxExp) );
    JredB = P.jstarB .* exp( clip(P.bB.*eta_eff, opt.maxExp) );
    Jox = JoxA + JoxB;
    Jred = JredA + JredB;

    eps0 = 1e-30;
    U = Jox + Jred;
    S = Jox - Jred;
    beta = S ./ (U + eps0);
    beta = max(min(beta,1-1e-12),-1+1e-12);

    I_store = 2*sqrt(max(Jox,0).*max(Jred,0)) ./ (U + eps0);
    I_store = max(min(I_store,1),0);

    % Cai–Smith Gamma
    chi = (Jred - Jox) ./ (Jred + Jox + eps0);
    chi = max(min(chi,1),-1);
    theta = acos(chi) - pi/2;                 % [-pi/2, pi/2]
    rho = sqrt( min(Jred,Jox) ./ (max(Jred,Jox)+eps0) );
    rho(~isfinite(rho)) = 0;
    rho = max(min(rho,1),0);
    Gamma = rho .* exp(1i*theta);

    % useful output & density (HER cathodic)
    Pi_use  = max(0, -jTot) .* max(1 - rho.^2, 0);
    Pi_dens = Pi_use ./ (abs(eta_eff) + opt.eta0);

    % mask
    mask = isfinite(real(Gamma)) & isfinite(imag(Gamma)) & isfinite(eta_eff) & isfinite(Pi_dens) & (Pi_dens > 0);
    if opt.cathodicOnly
        mask = mask & (eta_eff < 0);
    end

    % maximum Pi_dens on mask
    idxStar = NaN;
    Pi_dens_max = NaN;
    if any(mask)
        tmp = Pi_dens; tmp(~mask) = -Inf;
        [Pi_dens_max, idxStar] = max(tmp);
        if ~isfinite(Pi_dens_max) || Pi_dens_max <= 0
            Pi_dens_max = NaN; idxStar = NaN;
        end
    end

    % normalized color (0..1) within each pH
    r = Pi_dens ./ (Pi_dens_max + 1e-300);   % in [0,1] on mask
    r(~mask) = NaN;
    r = min(max(r,0),1);

    epsr = 1e-4;
    Cnorm = (log10(r + epsr) - log10(epsr)) / (log10(1 + epsr) - log10(epsr));
    Cnorm = min(max(Cnorm,0),1);
    OUT(g).C = Cnorm;

    % store
    OUT(g).pH = RES(g).pH;
    OUT(g).eta = etaGrid;
    OUT(g).eta_eff = eta_eff;
    OUT(g).jTot = jTot;
    OUT(g).Jox = Jox;
    OUT(g).Jred = Jred;
    OUT(g).U = U;
    OUT(g).S = S;
    OUT(g).beta = beta;
    OUT(g).I_store = I_store;
    OUT(g).rho = rho;
    OUT(g).theta = theta;
    OUT(g).Gamma = Gamma;
    OUT(g).Pi_use = Pi_use;
    OUT(g).Pi_dens = Pi_dens;
    OUT(g).C = Cnorm;
    OUT(g).mask = mask;
    OUT(g).idxStar = idxStar;
    OUT(g).Pi_dens_max = Pi_dens_max;
end

% ---------------- figure: 7 disks + empty tile 8 ----------------
nShow = min(nG,10);

ncol = 4;
nrow = 3;
width_cm = 24;
height_cm = 18;

fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(fig, nrow, ncol, 'Padding','compact', 'TileSpacing','compact');

% colormap
apply_colormap(fig, opt.colormapName);

axisCol = [0.65 0.65 0.65];
trajCol = [0.80 0.80 0.80];

for g = 1:nShow
    ax = nexttile(tlo, g);
    hold(ax,'on'); box(ax,'on');

    G = OUT(g).Gamma;
    xr = real(G); yi = imag(G);

    mk = OUT(g).mask & isfinite(OUT(g).C) & isfinite(xr) & isfinite(yi);
    okTraj = isfinite(xr) & isfinite(yi);

    % faint unit circle (recommended for publication)
    if opt.showUnitCircle
        th = linspace(0,2*pi,360);
        plot(ax, cos(th), sin(th), '-', 'Color',[0.87 0.87 0.87], 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    % optional grid (usually OFF)
    if opt.showGrid
        th = linspace(0,2*pi,360);
        for rr = [0.25 0.5 0.75 1.0]
            plot(ax, rr*cos(th), rr*sin(th), ':', 'Color',[0.75 0.75 0.75], 'LineWidth', 0.7, 'HandleVisibility','off');
        end
    end

    % trajectory light gray
    plot(ax, xr(okTraj), yi(okTraj), '-', 'Color', trajCol, 'LineWidth', 1.0, 'HandleVisibility','off');

    % colored points (normalized density)
    if any(mk)
        scatter(ax, xr(mk), yi(mk), opt.pointSize, OUT(g).C(mk), 'filled');
    else
        text(ax, 0.5, 0.5, 'no valid points', 'Units','normalized', ...
            'HorizontalAlignment','center', 'Color',[0.3 0.3 0.3]);
    end

    % fixed color scale for ALL tiles
    caxis(ax, [0 1]);

    % optimum marker + annotation (eta* and Pi_dens^max in absolute units)
    if isfinite(OUT(g).idxStar)
        k = OUT(g).idxStar;
        plot(ax, xr(k), yi(k), 'o', 'MarkerSize', 7.5, ...
            'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth', 1.2, ...
            'HandleVisibility','off');

        etaStar = OUT(g).eta_eff(k);
        PmaxAbs = OUT(g).Pi_dens_max;

        txt = sprintf('$\\eta^*_{\\rm eff}=%.2f\\,\\mathrm{V}$\n$\\Pi_{\\rm dens}^{\\max}=%.3g$', etaStar, PmaxAbs);
        text(ax, 0.3, 0.3, txt, 'Units','normalized', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Interpreter','latex', 'FontSize', 10, ...
            'BackgroundColor',[1 1 1 0.85], 'EdgeColor',[0.8 0.8 0.8], 'Margin', 3);
    end

    % zoom-to-data (optional)
    if opt.zoomToData && any(mk)
        xd = xr(mk); yd = yi(mk);
        xmin = min(xd); xmax = max(xd);
        ymin = min(yd); ymax = max(yd);
        cx = 0.5*(xmin+xmax); cy = 0.5*(ymin+ymax);
        span = max([xmax-xmin, ymax-ymin, 1e-3]);
        half = 0.5*span*(1 + 2*opt.zoomMargin);

        xL = [cx-half, cx+half];
        yL = [cy-half, cy+half];
        xL = max(min(xL, 1.05), -1.05);
        yL = max(min(yL, 1.05), -1.05);
        xlim(ax, xL); ylim(ax, yL);
    else
        xlim(ax, [-1.05 1.05]); ylim(ax, [-1.05 1.05]);
    end

    % orientation cross
    xl = xlim(ax); yl = ylim(ax);
    plot(ax, [0 0], yl, '-', 'Color', axisCol, 'LineWidth', 0.8, 'HandleVisibility','off');
    plot(ax, xl, [0 0], '-', 'Color', axisCol, 'LineWidth', 0.8, 'HandleVisibility','off');

    xlabel(ax,'Re(\Gamma_g)','Interpreter','tex');
    ylabel(ax,'Im(\Gamma_g)','Interpreter','tex');
    
    % Add label (a, b, c, ...) to each subplot
    text(-0.25, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
end

% tile 11: EMPTY (keep layout)
ax11 = nexttile(tlo, 11);
axis(ax11,'off');

% -------- single global colorbar on the far right --------
axCB = axes(fig, 'Position',[0.80 0.10 0.02 0.18], 'Visible','off');
apply_colormap(fig, opt.colormapName);
caxis(axCB, [0 1]);

cb = colorbar(axCB, 'Location','eastoutside');
cb.TickDirection = 'out';
cb.LineWidth = 0.8;
cb.FontSize = 12;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\Pi_{\mathrm{dens}}/\max(\Pi_{\mathrm{dens}})$';
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'0','0.5','1'};

% export
if ~isempty(opt.savePng)
    set(fig,'PaperPositionMode','auto');
    try
        exportgraphics(fig, opt.savePng, 'Resolution', 600);
    catch
        print(fig, opt.savePng, '-dpng', '-r600');
    end
    fprintf('Saved: %s\n', opt.savePng);
end
print(fig,'-dpng','-r600','Figure22.png');
fprintf('Saved: FigureS22.png\n');
end

% ===================== helpers =====================

function P = ensure_ab_fields(P)
% Ensure P has aA,bA,aB,bB. If missing, infer from (xi,lam) and f if present.
if ~isfield(P,'aA') || ~isfield(P,'bA') || ~isfield(P,'aB') || ~isfield(P,'bB')
    if isfield(P,'xiA') && isfield(P,'lamA') && isfield(P,'xiB') && isfield(P,'lamB') && isfield(P,'f')
        f = P.f;
    else
        if isfield(P,'T') && isfield(P,'F') && isfield(P,'Rg')
            f = P.F/(P.Rg*P.T);
        else
            error('P must contain aA,bA,aB,bB OR contain xi/lam and f (or T,F,Rg).');
        end
    end
    P.aA = P.xiA * P.lamA * f;
    P.bA = -(1 - P.xiA) * P.lamA * f;
    P.aB = P.xiB * P.lamB * f;
    P.bB = -(1 - P.xiB) * P.lamB * f;
end
end

function S = setdef(S, f, v)
if ~isfield(S,f) || isempty(S.(f)), S.(f) = v; end
end

function apply_colormap(fig, name)
try
    if exist(name,'file') == 2
        colormap(fig, feval(name, 256));
    else
        % fallbacks
        if strcmpi(name,'turbo') && exist('turbo','file')==2
            colormap(fig, turbo(256));
        else
            colormap(fig, parula(256));
        end
    end
catch
    colormap(fig, parula(256));
end
end