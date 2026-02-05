clc; clear; close all;
warning('off','all');

%% ---------------- Constant Definitions ----------------
T = 298.0; F = 96485.33212; Rg = 8.314462618;
f = F/(Rg*T);
nernst = log(10)/f;   % 2.303RT/F ~ 0.0591 V at 298K

%% ---------------- Control Parameters ----------------
maxExp = 70;

% Stage 1 parameters (fitting without ohmic drop)
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
prune.cathodicOnly = true;  % Only for cathodic processes
prune.jRelMin  = 1e-3;

%% ---------------- read data ----------------
csvFile = 'Fig3b_sourcedata.csv';
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
    
    % ===================== Fit AB model =====================
    fprintf('  Stage 1 (AB): Initial fitting without ohmic drop...\n');
    P1_AB = fit_stage1_noohm(etag, jg, Iref, f, maxExp, st1);
    
    fprintf('  Stage 2 (AB): Full implicit optimization...\n');
    P2_AB = fit_stage2_full(etag, jg, Iref, f, maxExp, st2, P1_AB);
    
    % Check and ensure mode A's jstar >= mode B's jstar
    P2_AB = ensure_modeA_primary(P2_AB);
    
    if prune.enable
        [~, infoAB] = decide_single_mode_by_contrib(P2_AB, etag, maxExp, ...
            st2.newtonMaxIter, st2.newtonTol, prune);
        fprintf('  [Diagnostic] Mode A contribution=%.3g, Mode B contribution=%.3g (valid points=%d)\n', ...
            infoAB.fracA, infoAB.fracB, infoAB.N);
    end
    
    % ===================== Fit single mode A/B =====================
    fprintf('  Stage 1 (A): Initial fitting without ohmic drop...\n');
    P1_A = fit_stage1_single_noohm(etag, jg, Iref, f, maxExp, st1, 'A', P2_AB);
    fprintf('  Stage 2 (A): Full implicit optimization...\n');
    P2_A = fit_stage2_single_full(etag, jg, Iref, f, maxExp, st2, P1_A, 'A');
    
    fprintf('  Stage 1 (B): Initial fitting without ohmic drop...\n');
    P1_B = fit_stage1_single_noohm(etag, jg, Iref, f, maxExp, st1, 'B', P2_AB);
    fprintf('  Stage 2 (B): Full implicit optimization...\n');
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
                choose = 2 + (AIC_B < AIC_A); % Select the better single mode
            end
        end
        
        if choose == 1
            P1 = P1_AB; P2 = P2_AB; modelTag = "AB";
        elseif choose == 2
            P1 = P1_A; P2 = P2_A; modelTag = "A";
        else
            P1 = P1_B; P2 = P2_B; modelTag = "B";
        end
        
        fprintf('  [AIC] AB=%.2f, A=%.2f, B=%.2f -> Selected %s\n', ...
            AIC_AB, AIC_A, AIC_B, modelTag);
    else
        P1 = P1_AB; P2 = P2_AB; modelTag = "AB";
    end
    
    % Store results
    RES(g) = struct('pH', ph, 'P', P2, 'P1', P1, 'Iref', Iref, ...
                    'U', Ug, 'eta', etag, 'j', jg, 'modelTag', modelTag);
    
    fprintf('  Final (%s): j*A=%.2e ξA=%.3f λA=%.3f jlimA=%.3g | ', ...
        modelTag, P2.jstarA, P2.xiA, P2.lamA, P2.jlimA);
    fprintf('j*B=%.2e ξB=%.3f λB=%.3f jlimB=%.3g | Rohm=%.3g\n', ...
        P2.jstarB, P2.xiB, P2.lamB, P2.jlimB, P2.Rohm);
end

% Open file to write parameters
fid = fopen('optimized_parameters_Au.txt', 'w');

% Write header
fprintf(fid, 'pH\tModel\tjstarA\txiA\tlamA\tjlimA\tjstarB\txiB\tlamB\tjlimB\tRohm\n');

% Loop through all results and write to file
for g = 1:nG
    res = RES(g);
    
    % Check P structure, assuming P is a struct containing fields
    if isstruct(res.P)
        % If it's a struct, access fields directly
        fprintf(fid, '%.2f\t%s\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f\t%.2e\n', ...
            res.pH, res.modelTag, ...
            res.P.jstarA, res.P.xiA, res.P.lamA, res.P.jlimA, ...
            res.P.jstarB, res.P.xiB, res.P.lamB, res.P.jlimB, ...
            res.P.Rohm);
    else
        % If P is an array or matrix, need to access based on actual storage structure
        fprintf('Warning: P for pH=%.2f is not a struct, using backup output method\n', res.pH);
        
        % Backup plan: adjust indices according to your actual data structure
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

% Overall fitting comparison plot
plotOpt2 = struct();
plotOpt2.titleTotal = 'Overall fitting vs data (all pH)';
plotOpt2.figFolder = 'figures';
plotOpt2.saveFigures = true;
plot_polarization_comparison_multi_pH(RES, plotOpt2, maxExp, ...
    st2.newtonMaxIter, st2.newtonTol);

%% ---------------- Cai-Smith plot analysis ----------------
caiOpt = struct();
caiOpt.useEtaEff = true;
caiOpt.phaseWeight = 'cos';
caiOpt.showGrid = true;
caiOpt.capGammaTo1 = false;
caiOpt.title = 'Cai-Smith plot (multi-pH overlay)';
caiOpt.plotStorageCurve = true;
caiOpt.storageUseEtaEff = true;
caiOpt.storageYScale = 'log';
caiOpt.nEta = 900;
caiOpt.edgeFrac = 0.01;
caiOpt.etaPctRange = [1 99];

% Assume the original program has been run, RES obtained
maxExp = 70;
newtonMaxIter = 10;
newtonTol = 5e-11;

% Call plotting function
plot_waveguide_electrochem_summary_figure(RES, maxExp, newtonMaxIter, newtonTol);

opt = struct();
opt.colorMode = 'Pi_dens';      % Or 'I_store' / 'rho' / 'theta'
opt.perTileColorbar = true;
opt.zoomToData = true;
opt.showUnitCircle = false;
opt.savePng = 'CS_electrochem.png';
plot_caismith_electrochem_waveguide(RES, opt);

%% ---------------- Waveguide theory advanced analysis ----------------
fprintf('\n=== Performing waveguide theory advanced analysis ===\n');

% 1. Verify four basic laws
% verify_four_laws(RES, maxExp, st2.newtonMaxIter, st2.newtonTol);

% 2. Feedback factor and critical coupling analysis
% analyze_feedback_coupling(RES, maxExp, st2.newtonMaxIter, st2.newtonTol);

% 3. Resonance and attenuation dynamics analysis
analyze_resonance_dynamics(RES, maxExp, st2.newtonMaxIter, st2.newtonTol);

% 4. Asymmetry parameter ξ analysis
% plot_asymmetry_parameter_analysis(RES);

% 5. 3D Cai-Smith visualization
% plot_3d_caismith(RES, caiOpt, maxExp, st2.newtonMaxIter, st2.newtonTol);

% Create custom plotting options
plotOpt = struct();
plotOpt.useEtaEff = true;           % Use effective overpotential
plotOpt.nPoints = 1000;            % Number of points per curve
plotOpt.etaRangePercentile = [1, 99]; % Use η 1%-99% range
plotOpt.saveFigures = true;        % Save figures
plotOpt.savePath = 'my_gamma_results'; % Save path
plotOpt.dpi = 600;                 % High resolution output

%% ---------------- Calculate waveguide relativistic parameters ----------------
% if ~isempty(RES)
%     % Select one pH condition for calculation
%     g = 3;  % Select first pH
%     P = RES(g).P;
% 
%     % Create overpotential grid
%     eta_range = linspace(min(RES(g).eta), max(RES(g).eta), 1000)';
% 
%     % Calculate new parameters
%     [Gmag, Gph, I_store, eta_eff, U, S_rel, beta_w, phi, ...
%      D_plus, D_minus, eta_trans] = ...
%         calc_caismith_from_params_improved(P, eta_range, true, 100, 50, 1e-6, 'none');
% 
%     %% ---------------- Visualization ----------------
%     visualize_waveguide_relativity_comprehensive(Gmag, Gph, I_store, eta_eff, ...
%                                                  U, S_rel, beta_w, phi, ...
%                                                  D_plus, D_minus, eta_trans);
% end
% optL = struct();
% optL.savePng = 'FigureS_laws_echem.png';
% optL.exampleIdx = 1;     % Which pH to use for Law4 scatter example
% LAW = verify_four_laws_electrochem_polarization(RES, optL);
% disp(LAW.table);

optC = struct();
optC.fixedEtaWindow = [-1.5 0.5];
optC.nEta = 700;
optC.eta0 = 0.10;
optC.cathodicOnly = true;

% Make Law4 more "significant" nonreciprocal adjustment (adjustable)
optC.nonrecipRotGain   = 0.35*pi;
optC.nonrecipPhaseGain = 0.25*pi;

% Tlaw4 = verify_law4_correlations_echem(RES, optC);

fprintf('\n=== Analysis complete ===\n');

%% =====================================================================
% Helper functions
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
        error('Some pH values could not be parsed, please check pH column.');
    end
end

function [jstar2, xi2, lam2, jlim2] = seed_second_mode_from_residual(eta, dj, f, maxExp, jmax)
    % Estimate initial parameters for second mode from residuals
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
% Basic visualization functions
%% =====================================================================

function plot_polarization_decomp(RES, plotOpt, maxExp)
% ===== Global defaults (put at the top of your script/function) =====
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
        
        % Format RMSE, keep 2 significant figures
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
        
        % Display pH and RMSE in top-left corner of subplot
        % text(0.05, 0.95, sprintf('pH=%.0f (RMSE=%s mA/cm²)', RES(g).pH, RMSE_str), ...
        %     'Units','normalized', 'FontSize', 10, 'FontWeight','normal', ...
        %     'VerticalAlignment','top', 'HorizontalAlignment','left', ...
        %     'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'Margin', 2);
        
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
    
    % If number of subplots <8, display RMSE statistics in last subplot
    if nG < 8
        % Determine position of last subplot
        last_tile = nG;
        
        % Get axis of last subplot
        ax_last = nexttile(last_tile);
        
        % Add RMSE statistics to current axis
        hold(ax_last, 'on');
        
        % Calculate statistics
        mean_RMSE = mean(RMSE_values);
        min_RMSE = min(RMSE_values);
        max_RMSE = max(RMSE_values);
        std_RMSE = std(RMSE_values);
        
        % Format statistics, keep 2 significant figures
        formatValue = @(x) sprintf('%.2f', x);
        
        % Create statistics text
        stats_text = {};
        
        % Add RMSE for each pH
        for g = 1:nG
            stats_text{end+1} = sprintf('  pH=%.1f: %s mA/cm²', RES(g).pH, formatValue(RMSE_values(g)));
        end
        
        % Add text to right side of axis
        text(ax_last, 1.2, 0.8, stats_text, ...
            'Units','normalized', 'FontSize', 12, 'FontWeight','normal', ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', ...
            'BackgroundColor', [1 1 1 0.9], 'EdgeColor', 'k', 'Margin', 3);
        
    else
        % If 8 or more subplots, display statistics in subplot 8
        stats_tile = min(8, nG);  % Use subplot 8, if nG<8 use last
        
        % Get axis of subplot
        ax_stats = nexttile(stats_tile);
        
        % Clear current axis
        cla(ax_stats);
        set(ax_stats, 'Visible', 'off');
        
        % Calculate statistics
        mean_RMSE = mean(RMSE_values);
        min_RMSE = min(RMSE_values);
        max_RMSE = max(RMSE_values);
        std_RMSE = std(RMSE_values);
        
        % Format statistics, keep 2 significant figures
        formatValue = @(x) sprintf('%.2g', x);
        
        % Create statistics text
        stats_text = {
            'RMSE statistics summary';
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
        
        % Add text at axis center
        % text(ax_stats, 0.5, 0.5, stats_text, ...
        %     'Units','normalized', 'FontSize', 10, 'FontWeight','normal', ...
        %     'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
        %     'BackgroundColor', [1 1 1 0.9], 'EdgeColor', 'k', 'Margin', 3);
        % 
        % Add label for statistics subplot
        text(ax_stats, -0.20, 0.98, char('a' + stats_tile - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
    end
    % Adjust aspect ratio of all subplots so statistics subplot not too small
    set(tlo, 'TileSpacing', 'compact', 'Padding', 'compact');
    filename = 'FigureS14';
print(gcf, '-dpng', '-r600', [filename, '.png']);
disp(['Figure saved as ', filename, '.png']);
end

function plot_polarization_comparison_multi_pH(RES, plotOpt, maxExp, newtonMaxIter, newtonTol)
    if nargin < 2 || isempty(plotOpt), plotOpt = struct(); end
    plotOpt = set_default(plotOpt, 'nFine', 700);
    plotOpt = set_default(plotOpt, 'saveFigures', true);
    plotOpt = set_default(plotOpt, 'figFolder', 'figures');
    plotOpt = set_default(plotOpt, 'titleTotal', 'Overall fitting vs data (all pH)');
    
    nG = numel(RES);
    C  = lines(nG);
    
    % Overall overlay plot
    fig1 = figure('Color','w', 'Name','Overall fitting overlay', 'Position',[90 70 980 700]);
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
         text(-0.20, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
        % ======== Display Mode A / Mode B Tafel slope in each subplot ========
% Show content: beta (mV/dec) + optional R2/N
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

% Place at top-left, two lines
x0 = 0.02; y0 = 0.20; dy = 0.08;  % Normalized coordinates
text(x0, y0,      txtA, 'Units','normalized', 'HorizontalAlignment','left', ...
    'VerticalAlignment','top', 'Color',[0.85 0 0], 'FontSize',10, 'FontWeight','bold');

text(x0, y0-dy,   txtB, 'Units','normalized', 'HorizontalAlignment','left', ...
    'VerticalAlignment','top', 'Color',[0 0.25 0.9], 'FontSize',10, 'FontWeight','bold');

        
        % title(sprintf('pH=%.3g (%s)', RES(g).pH, RES(g).modelTag), 'Interpreter','none');
        xlabel(xlab, 'Interpreter','latex');
        ylabel('$\log_{10}|j|$', 'Interpreter','latex');
        ylim([-12 2]);
        
        if g == 1
            legend1=legend({'$j_\mathrm{tot}$','$j_\mathrm{A}$','$j_\mathrm{B}$',...
                'Tafel $j_\mathrm{A}','Tafel $j_\mathrm{B}'}, 'Location','southwest',...
                'Interpreter','latex','FontSize',12);
            legend boxoff
            set(legend1,...
         'Position',[0.830256472856571 0.203370245546103 0.079556626539964 0.15634765625]);
        end
        
        pH_vals(g) = RES(g).pH;
    end
    
    TAF = table(pH_vals, betaA, R2A, NA, betaB, R2B, NB, ...
        'VariableNames', {'pH','betaA_mVdec','R2A','NA','betaB_mVdec','R2B','NB'});
    filename = 'FigureS15';
print(gcf, '-dpng', '-r600', [filename, '.png']);
disp(['Figure saved as ', filename, '.png']);
end

%% =====================================================================
% Cai-Smith plot analysis (improved version)
%% =====================================================================

function CaiTAB = plot_caismith_multi_pH_overlay_improved(RES, caiOpt, maxExp, newtonMaxIter, newtonTol)
    % Default parameter settings
    caiOpt = set_default(caiOpt, 'useEtaEff', true);
    caiOpt = set_default(caiOpt, 'phaseWeight', 'cos');
    caiOpt = set_default(caiOpt, 'showGrid', true);
    caiOpt = set_default(caiOpt, 'capGammaTo1', false);
    caiOpt = set_default(caiOpt, 'title', 'Cai-Smith plot (multi-pH overlay)');
    caiOpt = set_default(caiOpt, 'plotStorageCurve', true);
    caiOpt = set_default(caiOpt, 'storageUseEtaEff', true);
    caiOpt = set_default(caiOpt, 'storageYScale', 'log');
    caiOpt = set_default(caiOpt, 'nEta', 900);
    caiOpt = set_default(caiOpt, 'edgeFrac', 0.01);
    caiOpt = set_default(caiOpt, 'etaPctRange', [1 99]);
    
    nG = numel(RES);
    C  = lines(nG);
    LS = {'-','--',':','-.'};
    
    % Update CUR struct definition to include all new parameters
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
        
        % Store all parameters in CUR struct
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
        title(ax1, 'Normalized integrated storage intensity vs overpotential (normalized per pH)', 'Interpreter', 'none');
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
            ylabel(ax2, '$I_{\mathrm{store}}$ (absolute, log scale)', 'Interpreter', 'latex');
            title(ax2, 'Absolute storage intensity vs overpotential (semi-log)', 'Interpreter', 'none');
        else
            ylabel(ax2, '$I_{\mathrm{store}}$ (absolute)', 'Interpreter', 'latex');
            title(ax2, 'Absolute storage intensity vs overpotential', 'Interpreter', 'none');
        end
        
        sgtitle(tl, 'Storage intensity comparison across pH', 'Interpreter', 'none');
    end
    
    % Output table - updated to include complete new parameters
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
        % Keep old variable names for compatibility with calling code
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

    if nargin < 7 || isempty(phaseWeight), phaseWeight = 'none'; end %#ok<NASGU>

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
    chi = (Jox-Jred ) ./ (Jred + Jox + eps0);     % in [-1,1]
    chi = max(min(chi, 1), -1);

    varphi = acos(chi);                 % [0, pi]
    theta  = varphi - pi/2;             % [-pi/2, +pi/2] for Cai-Smith disk display

    rho = sqrt( min(Jred, Jox) ./ (max(Jred, Jox) + eps0) );  % [0,1]
    rho = sqrt((1-abs(beta_w))./(1+abs(beta_w)));
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
% Waveguide theory advanced analysis functions
%% =====================================================================

function verify_four_laws(RES, maxExp, newtonMaxIter, newtonTol)
% verify_four_laws (IMPROVED, sign-gated directional tests)
% ------------------------------------------------------------
% Key fixes vs previous version:
%   (1) Directional reflections are ONLY evaluated on their valid branches:
%       cathodic incidence: eta_eff < -eta_db
%       anodic   incidence: eta_eff > +eta_db
%       This avoids artificial |Gamma|>1 "warnings".
%   (2) Reciprocity pairing uses NEGATIVE branch (cathodic) vs POSITIVE branch (anodic)
%       at matched |eta_eff|, not (pos vs neg) which collapses to trivial zeros.
%   (3) Law 2 uses adaptive quantile bins in |beta_w|, and reports N/A if insufficient
%       bidirectional coverage exists.
%
% Outputs:
%   L1_A/L1_B : modal reciprocity mismatch (paired by |eta_eff|)
%   L2_beta   : equal-power-profile mismatch (binned by |beta_w|)
%   L3_total  : global modal balance mismatch (paired by |eta_eff|)
%   passivBad : true passivity issues within the VALID branch only
%
% Requires: safe_two_mode_model, clip, prctile_fallback

    fprintf('\n=== Verify four waveguide laws (directional, sign-gated) ===\n');

    if nargin < 2 || isempty(maxExp), maxExp = 70; end
    if nargin < 3 || isempty(newtonMaxIter), newtonMaxIter = 10; end
    if nargin < 4 || isempty(newtonTol), newtonTol = 5e-11; end

    eps0   = 1e-30;
    eta_db = 0.03;   % deadband to avoid near-zero ambiguous region
    nBins  = 6;      % adaptive bins in |beta_w|
    minPerBin = 4;   % relaxed threshold (was 5)
    minPair   = 10;  % min paired samples for |eta| tests

    nG = numel(RES);
    pH_vals   = nan(nG,1);

    L1_A     = nan(nG,1);
    L1_B     = nan(nG,1);
    L2_beta  = nan(nG,1);
    L3_total = nan(nG,1);
    L4_total = nan(nG,1);   % same as L3 here (directional reciprocity on totals)
    passivBad = zeros(nG,1);
    coverage  = strings(nG,1);

    EX = struct('ok',false); % for example plots

    for g = 1:nG
        P = RES(g).P;
        pH_vals(g) = RES(g).pH;

        eta = RES(g).eta(:);
        jTot = safe_two_mode_model(P, eta, maxExp, newtonMaxIter, newtonTol, []);
        eta_eff = eta - P.Rohm.*(jTot./1000);

        % --- unsaturated directional fluxes (mode-wise) ---
        JoxA  = P.jstarA .* exp( clip(P.aA.*eta_eff, maxExp) );
        JredA = P.jstarA .* exp( clip(P.bA.*eta_eff, maxExp) );
        JoxB  = P.jstarB .* exp( clip(P.aB.*eta_eff, maxExp) );
        JredB = P.jstarB .* exp( clip(P.bB.*eta_eff, maxExp) );

        % totals and beta_w
        Jox = JoxA + JoxB;
        Jred = JredA + JredB;
        U = Jox + Jred;
        S = Jox - Jred;
        beta = S ./ (U + eps0);
        beta = max(min(beta, 1-1e-12), -1+1e-12);
        betaAbs = abs(beta);

        % valid branch masks
        mC = (eta_eff < -eta_db); % cathodic branch
        mA = (eta_eff >  eta_db); % anodic branch

        if ~any(mC) && ~any(mA)
            coverage(g) = "none";
            continue;
        elseif any(mC) && ~any(mA)
            coverage(g) = "cathodic-only";
        elseif ~any(mC) && any(mA)
            coverage(g) = "anodic-only";
        else
            coverage(g) = "bidirectional";
        end

        % -------- directional amplitude reflections (ONLY on their branch) --------
        % cathodic incidence on mC: Gamma_c = sqrt(Jox/Jred) should be <=1 typically
        GamC_A = sqrt( max(JoxA,0) ./ (max(JredA,0)+eps0) );
        GamC_B = sqrt( max(JoxB,0) ./ (max(JredB,0)+eps0) );
        % anodic incidence on mA: Gamma_a = sqrt(Jred/Jox) should be <=1 typically
        GamA_A = sqrt( max(JredA,0) ./ (max(JoxA,0)+eps0) );
        GamA_B = sqrt( max(JredB,0) ./ (max(JoxB,0)+eps0) );

        % absorptance per mode (clip), but evaluate on the correct branch only
        aC_A = nan(size(eta_eff)); aC_B = aC_A;
        aA_A = nan(size(eta_eff)); aA_B = aC_A;

        aC_A(mC) = clip01(1 - min(1, GamC_A(mC).^2));
        aC_B(mC) = clip01(1 - min(1, GamC_B(mC).^2));
        aA_A(mA) = clip01(1 - min(1, GamA_A(mA).^2));
        aA_B(mA) = clip01(1 - min(1, GamA_B(mA).^2));

        aC_tot = clip01(aC_A + aC_B);
        aA_tot = clip01(aA_A + aA_B);

        % true passivity warnings: only count |Gamma|>1 INSIDE the valid branch
% ---- robust passivity warnings: force column vectors before logical ops ----
GcA = GamC_A(:);  GcB = GamC_B(:);
GaA = GamA_A(:);  GaB = GamA_B(:);
mC  = mC(:);      mA  = mA(:);

badC = false(size(mC));
badA = false(size(mA));

if any(mC)
    badC = (GcA(mC) > 1+1e-6) | (GcB(mC) > 1+1e-6);
end
if any(mA)
    badA = (GaA(mA) > 1+1e-6) | (GaB(mA) > 1+1e-6);
end

passivBad(g) = nnz(badC) + nnz(badA);


        % -------- Laws 1/3/4: pair NEG branch (cathodic) vs POS branch (anodic) by |eta_eff| --------
        if any(mC) && any(mA)
            [eta_abs, cA, aA] = pair_neg_pos_by_absEta(eta_eff, aC_A, aA_A);
            [~,      cB, aB] = pair_neg_pos_by_absEta(eta_eff, aC_B, aA_B);
            [~,     cT, aT]  = pair_neg_pos_by_absEta(eta_eff, aC_tot, aA_tot);

            if numel(eta_abs) >= minPair
                L1_A(g)     = rel_mismatch(cA, aA);
                L1_B(g)     = rel_mismatch(cB, aB);
                L3_total(g) = rel_mismatch(cT, aT);
                L4_total(g) = L3_total(g);
            end
        end

        % -------- Law 2: equal-power-profile via adaptive quantile bins in |beta_w| (branch-wise) --------
        if any(mC) && any(mA)
            L2_beta(g) = beta_binned_mismatch_quantile(betaAbs, mC, mA, aC_tot, aA_tot, nBins, minPerBin);
        end

        % store one example (first bidirectional)
        if ~EX.ok && any(mC) && any(mA)
            [eta_abs, cT, aT] = pair_neg_pos_by_absEta(eta_eff, aC_tot, aA_tot);
            [bb, mc, ma, acb, aab] = beta_binned_curves_quantile(betaAbs, mC, mA, aC_tot, aA_tot, nBins, minPerBin);
            EX.ok = true;
            EX.pH = RES(g).pH;
            EX.eta_abs = eta_abs; EX.cT = cT; EX.aT = aT;
            EX.bb = bb; EX.acb = acb; EX.aab = aab;
        end

        fprintf('pH=%.3g: cov=%s | L1(A)=%.3e L1(B)=%.3e L2=%.3e L3=%.3e passivBad=%d\n', ...
            RES(g).pH, coverage(g), L1_A(g), L1_B(g), L2_beta(g), L3_total(g), passivBad(g));
    end

    % ---------------- plotting ----------------
    figure('Color','w','Name','Electrochemical waveguide-law verification (improved)','Position',[80 60 1500 900]);
    x = 1:nG;
    xt = arrayfun(@(v) sprintf('%.2g',v), pH_vals, 'UniformOutput', false);

    subplot(2,3,1); hold on; grid on; box on;
    bh1 = bar(x-0.18, L1_A, 0.35); %#ok<NASGU>
    bh2 = bar(x+0.18, L1_B, 0.35); %#ok<NASGU>
    set(gca,'XTick',x,'XTickLabel',xt); xtickangle(45);
    ylabel('relative mismatch'); title('Law 1: modal reciprocity (neg vs pos, matched |η_{eff}|)','Interpreter','tex');
    legend({'mode A','mode B'},'Location','best');

    subplot(2,3,2); hold on; grid on; box on;
    bar(x, L2_beta, 0.65);
    set(gca,'XTick',x,'XTickLabel',xt); xtickangle(45);
    ylabel('relative mismatch'); title('Law 2: equal-power-profile (adaptive bins in |β_w|)','Interpreter','tex');

    subplot(2,3,3); hold on; grid on; box on;
    bar(x, L3_total, 0.65);
    set(gca,'XTick',x,'XTickLabel',xt); xtickangle(45);
    ylabel('relative mismatch'); title('Law 3: global balance (Σ_p α_c,p vs Σ_p α_a,p)','Interpreter','tex');

    subplot(2,3,4); hold on; grid on; box on;
    bar(x, passivBad, 0.65);
    set(gca,'XTick',x,'XTickLabel',xt); xtickangle(45);
    ylabel('count'); title('passivity warnings within valid branches','Interpreter','none');

    subplot(2,3,5); hold on; grid on; box on;
    if EX.ok && ~isempty(EX.eta_abs)
        scatter(EX.cT, EX.aT, 30, EX.eta_abs, 'filled');
        plot([0,1],[0,1],'k--','LineWidth',1.5);
        xlabel('\alpha_c(|\eta_{eff}|)  (cathodic branch)','Interpreter','tex');
        ylabel('\alpha_a(|\eta_{eff}|)  (anodic branch)','Interpreter','tex');
        title(sprintf('Example pH=%.2g: Law 4 reciprocity scatter', EX.pH),'Interpreter','none');
        axis equal; xlim([0,1]); ylim([0,1]);
        colorbar;
    else
        axis off; text(0.5,0.5,'No bidirectional coverage (HER-only dataset): reciprocity not testable','HorizontalAlignment','center');
    end

    subplot(2,3,6); hold on; grid on; box on;
    if EX.ok && ~isempty(EX.bb)
        plot(EX.bb, EX.acb, 'o-', 'LineWidth', 2);
        plot(EX.bb, EX.aab, 's-', 'LineWidth', 2);
        xlabel('|β_w| bin center','Interpreter','tex');
        ylabel('<α>'); ylim([0,1]);
        title(sprintf('Example pH=%.2g: Law 2 (profile equivalence)', EX.pH),'Interpreter','none');
        legend({'cathodic branch','anodic branch'},'Location','best');
    else
        axis off; text(0.5,0.5,'No valid beta-binned curves (insufficient bidirectional points)','HorizontalAlignment','center');
    end

    sgtitle('Electrochemical waveguide-law verification via directional reflections (sign-gated)', ...
        'FontSize', 14, 'FontWeight','bold');

    % summary
    fprintf('\n=== Summary (omit NaN) ===\n');
    fprintf('Law1(A): %.3e ± %.3e\n', mean(L1_A,'omitnan'), std(L1_A,'omitnan'));
    fprintf('Law1(B): %.3e ± %.3e\n', mean(L1_B,'omitnan'), std(L1_B,'omitnan'));
    fprintf('Law2(beta): %.3e ± %.3e\n', mean(L2_beta,'omitnan'), std(L2_beta,'omitnan'));
    fprintf('Law3(total): %.3e ± %.3e\n', mean(L3_total,'omitnan'), std(L3_total,'omitnan'));
    fprintf('Passivity warnings (median): %.1f\n', median(passivBad));

    % ------------ helpers (nested) ------------
    function y = clip01(x), y = min(max(x,0),1); end

    function m = rel_mismatch(a, b)
        a = a(:); b = b(:);
        good = isfinite(a) & isfinite(b);
        a = a(good); b = b(good);
        if numel(a) < minPair, m = NaN; return; end
        denom = max(max(a,b), 1e-12);
        m = mean(abs(a-b)./denom, 'omitnan');
    end

    function [eta_abs, c_neg, a_pos] = pair_neg_pos_by_absEta(eta_eff, alpha_c, alpha_a)
        eta_eff = eta_eff(:); alpha_c = alpha_c(:); alpha_a = alpha_a(:);
        neg = eta_eff < 0; pos = eta_eff > 0;
        if ~any(neg) || ~any(pos), eta_abs=[]; c_neg=[]; a_pos=[]; return; end

        % cathodic uses NEG side (abs value)
        en = -eta_eff(neg);  cn = alpha_c(neg);
        ep =  eta_eff(pos);  ap = alpha_a(pos);

        % unique for interp1
        [en, in] = unique(en,'stable'); cn = cn(in);
        [ep, ip] = unique(ep,'stable'); ap = ap(ip);

        lo = max(min(en), min(ep));
        hi = min(max(en), max(ep));
        if ~(isfinite(lo)&&isfinite(hi)&&hi>lo), eta_abs=[]; c_neg=[]; a_pos=[]; return; end

        eta_abs = linspace(lo, hi, 80).';
        c_neg = interp1(en, cn, eta_abs, 'linear', 'extrap');
        a_pos = interp1(ep, ap, eta_abs, 'linear', 'extrap');

        c_neg = clip01(c_neg); a_pos = clip01(a_pos);
    end

    function mismatch = beta_binned_mismatch_quantile(betaAbs, mC, mA, aC_tot, aA_tot, nBins, minPerBin)
        [bb, ~, ~, acb, aab] = beta_binned_curves_quantile(betaAbs, mC, mA, aC_tot, aA_tot, nBins, minPerBin);
        if isempty(bb), mismatch = NaN; return; end
        mismatch = rel_mismatch(acb, aab);
    end

    function [bins, mC_use, mA_use, aC_bin, aA_bin] = beta_binned_curves_quantile(betaAbs, mC, mA, aC_tot, aA_tot, nBins, minPerBin)
        betaAbs = betaAbs(:); aC_tot = aC_tot(:); aA_tot = aA_tot(:);

        mC_use = mC & isfinite(betaAbs) & isfinite(aC_tot);
        mA_use = mA & isfinite(betaAbs) & isfinite(aA_tot);

        if nnz(mC_use) < (minPerBin*2) || nnz(mA_use) < (minPerBin*2)
            bins=[]; aC_bin=[]; aA_bin=[]; return;
        end

        % quantile edges based on pooled |beta|
        bb_pool = betaAbs(mC_use | mA_use);
        q = linspace(0,1,nBins+1);
        edges = quantile(bb_pool, q);
        edges(1)=0; edges(end)=1;
        edges = unique(edges);
        if numel(edges) < 3
            bins=[]; aC_bin=[]; aA_bin=[]; return;
        end

        centers = 0.5*(edges(1:end-1)+edges(2:end));
        aC_bin = nan(numel(centers),1);
        aA_bin = nan(numel(centers),1);

        for k = 1:numel(centers)
            inC = mC_use & (betaAbs >= edges(k)) & (betaAbs < edges(k+1));
            inA = mA_use & (betaAbs >= edges(k)) & (betaAbs < edges(k+1));
            if nnz(inC) >= minPerBin, aC_bin(k) = mean(aC_tot(inC),'omitnan'); end
            if nnz(inA) >= minPerBin, aA_bin(k) = mean(aA_tot(inA),'omitnan'); end
        end

        good = isfinite(aC_bin) & isfinite(aA_bin);
        bins = centers(good).';
        aC_bin = aC_bin(good);
        aA_bin = aA_bin(good);
    end
end

% ======================================================================
% Power absorption rate based on Gamma_g: alpha = 1 - |Gamma_g|^2, epsilon = alpha
% ======================================================================
% ======================================================================
% Power absorption rate based on Gamma_g (j as field amplitude):
%   Af = Jf, Ab = Jb
%   eta_eff < 0:  Gamma_g = Af/Ab
%   eta_eff > 0:  Gamma_g = Ab/Af
%   alpha = 1 - |Gamma_g|^2,  epsilon = alpha
% ======================================================================
function [alpha, epsilon] = calculate_absorption_emission_gammaPower(P, eta_eff, mode, maxExp)
    % mode='A' or 'B': do single-channel Gamma on polarization-only definition
    eta_eff = eta_eff(:);
    eps0 = 1e-30;

    if mode=='A'
        Jox  = P.jstarA .* exp(clip(P.aA.*eta_eff, maxExp));
        Jred = P.jstarA .* exp(clip(P.bA.*eta_eff, maxExp));
    else
        Jox  = P.jstarB .* exp(clip(P.aB.*eta_eff, maxExp));
        Jred = P.jstarB .* exp(clip(P.bB.*eta_eff, maxExp));
    end

    rho = sqrt( min(Jred,Jox) ./ (max(Jred,Jox)+eps0) );
    rho(~isfinite(rho)) = 0;

    alpha = 1 - rho.^2;
    alpha = min(max(alpha,0),1);
    epsilon = alpha;
end

% ======================================================================
% Law 4: Compare alpha_+(|η|) with epsilon_-(|η|) after aligning |η|
% Return average relative difference
% ======================================================================
function mismatch = reciprocity_mismatch_absEta(eta_eff, alpha_tot, eps_tot)
    eta_eff   = eta_eff(:);
    alpha_tot = alpha_tot(:);
    eps_tot   = eps_tot(:);

    pos = eta_eff > 0;
    neg = eta_eff < 0;

    if ~any(pos) || ~any(neg)
        mismatch = NaN;
        return;
    end

    [eta_abs_grid, alpha_pos_i, eps_neg_i] = pair_pos_neg_by_absEta(eta_eff, alpha_tot, eps_tot);

    denom = max(alpha_pos_i, eps_neg_i) + eps;
    mismatch = mean(abs(alpha_pos_i - eps_neg_i) ./ denom, 'omitnan');
end

function [eta_abs_grid, alpha_pos_i, eps_neg_i] = pair_pos_neg_by_absEta(eta_eff, alpha_tot, eps_tot)
    pos = eta_eff > 0;
    neg = eta_eff < 0;

    eta_pos = eta_eff(pos);
    a_pos   = alpha_tot(pos);

    eta_neg_abs = -eta_eff(neg);      % Take |η| on negative side
    e_neg  = eps_tot(neg);

    % Deduplicate and sort (avoid interp1 error)
    [eta_pos_u, ia] = unique(eta_pos, 'stable'); a_pos_u = a_pos(ia);
    [eta_neg_u, in] = unique(eta_neg_abs, 'stable'); e_neg_u = e_neg(in);

    eta_abs_grid = linspace(max(min(eta_pos_u), min(eta_neg_u)), ...
                            min(max(eta_pos_u), max(eta_neg_u)), 80).';

    if numel(eta_abs_grid) < 5 || any(~isfinite(eta_abs_grid))
        eta_abs_grid = [];
        alpha_pos_i = [];
        eps_neg_i = [];
        return;
    end

    alpha_pos_i = interp1(eta_pos_u, a_pos_u, eta_abs_grid, 'linear', 'extrap');
    eps_neg_i   = interp1(eta_neg_u, e_neg_u, eta_abs_grid, 'linear', 'extrap');

    % clip to [0,1] to avoid numerical extrapolation out of bounds
    alpha_pos_i = min(max(alpha_pos_i, 0), 1);
    eps_neg_i   = min(max(eps_neg_i, 0), 1);
end

function analyze_feedback_coupling(RES, maxExp, newtonMaxIter, newtonTol)
    % Analyze feedback factor κ and critical coupling
    
    fprintf('\n=== Feedback factor and critical coupling analysis ===\n');
    
    figure('Color','w', 'Name','Feedback and critical coupling analysis', 'Position',[100 100 1200 500]);
    
    for g = 1:numel(RES)
        P = RES(g).P;
        eta = RES(g).eta;
        
        % Calculate effective overpotential and current
        jTot = safe_two_mode_model(P, eta, maxExp, newtonMaxIter, newtonTol, []);
        eta_eff = eta - P.Rohm.*(jTot./1000);
        
        % Calculate generalized reflection coefficient
        [Gmag, ~, ~] = calc_caismith_from_params_improved(P, eta, true, maxExp, newtonMaxIter, newtonTol, 'none');
        
        % Calculate feedback factor κ = |Γ| * exp(-2R)
        % R is attenuation, approximated by absolute value of η_eff
        R = mean(abs(eta_eff)) * 0.1; % Simplified estimate
        kappa = Gmag .* exp(-2 * R);
        
        % Find critical coupling point (κ closest to 1)
        [~, critical_idx] = min(abs(kappa - 1));
        
        subplot(2, ceil(numel(RES)/2), g);
        hold on; grid on; box on;
        
        plot(eta_eff, kappa, 'b-', 'LineWidth', 2);
        if ~isempty(critical_idx)
            plot(eta_eff(critical_idx), kappa(critical_idx), 'ro', ...
                'MarkerSize', 10, 'MarkerFaceColor', 'r');
            text(eta_eff(critical_idx), kappa(critical_idx), 'Critical coupling', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        end
        
        xlabel('$\eta_{\mathrm{eff}}$ (V)', 'Interpreter', 'latex');
        ylabel('$\kappa$ (feedback factor)', 'Interpreter', 'latex');
        title(sprintf('pH=%.3g - Feedback factor evolution', RES(g).pH));
        ylim([0 1.5]);
        plot(xlim, [1 1], 'k--', 'LineWidth', 1);
    end
    
    sgtitle('Feedback factor κ evolution and critical coupling point identification');
end

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

    eta0 = 0.10;
    fixedEtaWindow = [-1.5, 0.5];
    nEta = 700;

    nG = numel(RES);
    nShow = min(7, nG);

    OUT = repmat(struct('pH',NaN,'eta_eff_s',[],'rho_s',[],'tau_ec',[], ...
                        'Istore_s',[],'dchi_abs',[],'eta_eff2s',[],'Pi_dens',[], ...
                        'maskC',[],'etaStar',NaN), nShow, 1);

    for g = 1:nShow
        P = RES(g).P;

        etaGrid = linspace(fixedEtaWindow(1), fixedEtaWindow(2), nEta).';

        % NOTE: we keep this call to obtain rho, theta, I_store, eta_eff quickly
        [rho, theta, I_store, eta_eff] = ...
            calc_caismith_from_params_improved(P, etaGrid, true, maxExp, newtonMaxIter, newtonTol, 'none');

        rho     = max(min(rho(:),1),0);
        theta   = unwrap(theta(:)); %#ok<NASGU>
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
        chi = compute_chi_from_params_unsat(P, eta_eff, maxExp); % same length as eta_eff
        chi_s = chi(ord);
        dchi = gradient(chi_s, eta_eff_s);
        dchi_abs = abs(dchi);
        % optional smoothing:
        % dchi_abs = movmean(dchi_abs, 9);

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
    nrow = 2;  % fixed 2x4
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

        % if isfinite(etaStar)
        %     xline(ax, etaStar, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off');
        % end

        text(ax, -0.25, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k'); %#ok<*NBRAK>
        % (If MATLAB complains about "Color mocking", change to 'Color','k')
    end

    % ---- tile 8 legend (Fig 1) ----
    axL1 = nexttile(tlo1, 8);
    axis(axL1,'off'); hold(axL1,'on');

    h1 = plot(axL1, NaN, NaN, 'b-', 'LineWidth', 2);
    h2 = plot(axL1, NaN, NaN, 'r-', 'LineWidth', 2);
    % h3 = plot(axL1, NaN, NaN, 'k--', 'LineWidth', 2.0);

    lg1 = legend(axL1, [h1 h2], ...
        {'$\rho=|\Gamma_g|$', '$\tilde{\tau}_{\rm ec}=1/(-\ln\rho)$'}, ...
        'Interpreter','latex', 'Location','best');
    lg1.Box = 'off';
    lg1.FontSize = 12;

    filename = 'FigureS17';
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

    % ---- tile 8 legend (Fig 2) ----
    axL2 = nexttile(tlo2, 8);
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

    filename = 'FigureS18';
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, [filename,'.png']);

end

% =====================================================================
% helper: chi from UNSATURATED directional fluxes (consistent with Cai–Smith defs)
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

function plot_asymmetry_parameter_analysis(RES)
    % Asymmetry parameter ξ analysis
    
    figure('Color','w', 'Name','Asymmetry parameter analysis', 'Position',[100 100 1000 400]);
    
    subplot(1, 2, 1);
    xiA_values = arrayfun(@(x) x.P.xiA, RES);
    plot([RES.pH], xiA_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('pH');
    ylabel('$\xi_A$', 'Interpreter', 'latex');
    title('Mode A asymmetry parameter $\xi$', 'Interpreter', 'latex');
    grid on; box on;
    
    subplot(1, 2, 2);
    xiB_values = arrayfun(@(x) x.P.xiB, RES);
    plot([RES.pH], xiB_values, 's-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('pH');
    ylabel('$\xi_B$', 'Interpreter', 'latex');
    title('Mode B asymmetry parameter $\xi$', 'Interpreter', 'latex');
    grid on; box on;
    
    sgtitle('Asymmetry parameter $\xi$ vs pH', 'Interpreter', 'latex');
end

function plot_3d_caismith(RES, caiOpt, maxExp, newtonMaxIter, newtonTol)
    % Plot 3D Cai-Smith plot
    
    figure('Color','w', 'Name','3D Cai-Smith visualization', 'Position',[100 100 1200 800]);
    
    for g = 1:min(4, numel(RES))
        P = RES(g).P;
        eta = linspace(min(RES(g).eta), max(RES(g).eta), 200)';
        
        [Gmag, Gph, S] = calc_caismith_from_params_improved(P, eta, true, maxExp, newtonMaxIter, newtonTol, 'cos');
        
        % Convert to Cartesian coordinates
        x = Gmag .* cos(Gph);
        y = Gmag .* sin(Gph);
        z = S / max(S); % Normalized storage intensity
        
        subplot(2, 2, g);
        scatter3(x, y, z, 40, eta, 'filled');
        hold on;
        
        % Plot unit circle
        theta = linspace(0, 2*pi, 100);
        plot3(cos(theta), sin(theta), zeros(size(theta)), 'k--', 'LineWidth', 1);
        
        xlabel('Re($\Gamma_g$)', 'Interpreter', 'latex');
        ylabel('Im($\Gamma_g$)', 'Interpreter', 'latex');
        zlabel('Normalized storage intensity');
        title(sprintf('pH=%.3g - 3D Cai-Smith plot', RES(g).pH));
        grid on; box on;
        view(45, 30);
        colormap(jet);
        colorbar;
    end
    
    sgtitle('3D Cai-Smith plot: Generalized reflection coefficient vs storage intensity');
end

%% =====================================================================
% Fitting functions (unchanged)
%% =====================================================================

function P1 = fit_stage1_noohm(eta, jObs, Iref, f, maxExp, st1)
    % Stage 1: fitting without ohmic drop
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
    % Stage 2: full implicit fitting
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
% Model evaluation functions
%% =====================================================================

function jpred = safe_two_mode_model(P, eta, maxExp, newtonMaxIter, newtonTol, j0)
    % Enhanced error handling model evaluation function
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
    % Vectorized Newton method to solve implicit equation
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
    
    % Check if there is a mode swap flag
    if isfield(P, 'modesSwapped') && P.modesSwapped
        % Modes have been swapped, mode A is now the original mode B
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
% Single mode fitting functions
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
% AIC and diagnostic functions
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
% Tafel fitting helper functions
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
% Mathematical and utility functions
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
        warning('Cannot create folder "%s": %s. Using temp directory.', folderAbs, msg);
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
% Waveguide mass-energy mapping function (optimized)
%% =====================================================================

function [j, djdx] = wg_mass_energy_map_optimized(x, jlim)
    % Optimized mass-energy mapping function
    jlim = max(jlim, 1e-12);
    z = x ./ jlim;
    
    % Use hybrid method for better numerical stability
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
% Simple getter functions
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

function visualize_gamma_results(RES, maxExp, newtonMaxIter, newtonTol, plotOpt)
% Visualize Gmag and Gph results
% Input:
%   RES: struct array containing fitting results
%   maxExp, newtonMaxIter, newtonTol: model calculation parameters
%   plotOpt: plotting options struct

    if nargin < 5
        plotOpt = struct();
    end
    
    % Set default plotting options
    plotOpt = set_default(plotOpt, 'useEtaEff', true);
    plotOpt = set_default(plotOpt, 'nPoints', 500);
    plotOpt = set_default(plotOpt, 'etaRangePercentile', [5, 95]);
    plotOpt = set_default(plotOpt, 'saveFigures', false);
    plotOpt = set_default(plotOpt, 'savePath', 'gamma_visualization');
    plotOpt = set_default(plotOpt, 'dpi', 300);
    
    nG = numel(RES);
    colors = lines(nG);  % Assign color to each pH
    
    % Create output directory
    if plotOpt.saveFigures && ~exist(plotOpt.savePath, 'dir')
        mkdir(plotOpt.savePath);
    end
    
    %% Figure 1: Gmag and Gph vs η_eff (separate subplot for each pH)
    fig1 = figure('Name', 'Gmag and Gph vs Eta_eff (by pH)', ...
                  'Color', 'w', 'Position', [100, 100, 1400, 800]);
    
    % Calculate subplot layout
    nRows = ceil(sqrt(nG));
    nCols = ceil(nG / nRows);
    
    % Store all data min/max for uniform axes
    allGmag = [];
    allGphDeg = [];
    allEtaEff = [];
    
    for g = 1:nG
        % Get raw η data for this pH
        eta_raw = RES(g).eta;
        
        % Determine plotting range (ignore edge outliers)
        lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
        hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        
        % Create refined η grid
        eta_grid = linspace(lo, hi, plotOpt.nPoints)';
        
        % Calculate Gmag and Gph
        [Gmag, Gph, ~, eta_eff] = calc_caismith_from_params_improved(...
            RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
        
        % Convert phase to degrees (-180° to 180°)
        Gph_deg = Gph * 180/pi;
        Gph_deg = mod(Gph_deg + 180, 360) - 180;  % Wrap to [-180, 180]
        
        % Store data for later uniform axes
        allGmag = [allGmag; Gmag];
        allGphDeg = [allGphDeg; Gph_deg];
        allEtaEff = [allEtaEff; eta_eff];
        
        % Create subplot
        subplot(nRows, nCols, g);
        
        % Plot Gmag (left y-axis)
        yyaxis left;
        plot(eta_eff, Gmag, 'b-', 'LineWidth', 2);
        ylabel('|\Gamma_g|', 'FontSize', 12);
        ylim([0, 1.1]);  % Gmag theoretical range [0,1]
        grid on;
        
        % Plot Gph (right y-axis)
        yyaxis right;
        plot(eta_eff, Gph_deg, 'r-', 'LineWidth', 2);
        ylabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
        
        % Set title and labels
        title(sprintf('pH = %.2f', RES(g).pH), 'FontSize', 14);
        xlabel('\eta_{eff} (V)', 'FontSize', 12);
        
        % Add grid and legend
        grid on;
        legend({'|\Gamma_g|', '\angle\Gamma_g'}, 'Location', 'best');
        
        % Add horizontal reference lines
        yyaxis left;
        hold on;
        plot(xlim, [1, 1], 'b--', 'LineWidth', 0.5);  % Gmag=1 reference
        yyaxis right;
        plot(xlim, [0, 0], 'r--', 'LineWidth', 0.5);  % Gph=0 reference
    end
    
    sgtitle('Generalized reflection coefficient \Gamma_g vs effective overpotential', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure 1
    if plotOpt.saveFigures
        saveas(fig1, fullfile(plotOpt.savePath, 'gamma_by_pH.png'));
        savefig(fig1, fullfile(plotOpt.savePath, 'gamma_by_pH.fig'));
    end
    
    %% Figure 2: All pH Gmag and Gph overlay
    fig2 = figure('Name', 'Gmag and Gph Overlay (all pH)', ...
                  'Color', 'w', 'Position', [100, 100, 1200, 500]);
    
    % Subplot 1: All Gmag overlay
    subplot(1, 2, 1);
    hold on; grid on; box on;
    
    for g = 1:nG
        eta_raw = RES(g).eta;
        lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
        hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        eta_grid = linspace(lo, hi, plotOpt.nPoints)';
        
        [Gmag, ~, ~, eta_eff] = calc_caismith_from_params_improved(...
            RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
        
        plot(eta_eff, Gmag, '-', 'LineWidth', 2, ...
             'Color', colors(g, :), 'DisplayName', sprintf('pH=%.2f', RES(g).pH));
    end
    
    xlabel('\eta_{eff} (V)', 'FontSize', 12);
    ylabel('|\Gamma_g|', 'FontSize', 12);
    title('Generalized reflection coefficient magnitude |\Gamma_g|', 'FontSize', 14);
    legend('Location', 'best');
    ylim([0, 1.1]);
    plot(xlim, [1, 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Subplot 2: All Gph overlay (converted to degrees)
    subplot(1, 2, 2);
    hold on; grid on; box on;
    
    for g = 1:nG
        eta_raw = RES(g).eta;
        lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
        hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        eta_grid = linspace(lo, hi, plotOpt.nPoints)';
        
        [~, Gph, ~, eta_eff] = calc_caismith_from_params_improved(...
            RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
        
        Gph_deg = Gph * 180/pi;
        Gph_deg = mod(Gph_deg + 180, 360) - 180;
        
        plot(eta_eff, Gph_deg, '-', 'LineWidth', 2, ...
             'Color', colors(g, :), 'DisplayName', sprintf('pH=%.2f', RES(g).pH));
    end
    
    xlabel('\eta_{eff} (V)', 'FontSize', 12);
    ylabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
    title('Generalized reflection coefficient phase \angle\Gamma_g', 'FontSize', 14);
    legend('Location', 'best');
    plot(xlim, [0, 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    sgtitle('Generalized reflection coefficient for all pH conditions', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure 2
    if plotOpt.saveFigures
        saveas(fig2, fullfile(plotOpt.savePath, 'gamma_overlay.png'));
        savefig(fig2, fullfile(plotOpt.savePath, 'gamma_overlay.fig'));
    end
    
    %% Figure 3: Gmag-Gph parameter space trajectory (each pH)
    fig3 = figure('Name', 'Gamma Parametric Trajectories', ...
                  'Color', 'w', 'Position', [100, 100, 1200, 800]);
    
    for g = 1:nG
        subplot(nRows, nCols, g);
        hold on; grid on; box on;
        
        eta_raw = RES(g).eta;
        lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
        hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        eta_grid = linspace(lo, hi, plotOpt.nPoints)';
        
        [Gmag, Gph, ~, eta_eff] = calc_caismith_from_params_improved(...
            RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
        
        % Convert to degrees
        Gph_deg = Gph * 180/pi;
        Gph_deg = mod(Gph_deg + 180, 360) - 180;
        
        % Parametric plot: color represents η_eff
        scatter(Gph_deg, Gmag, 30, eta_eff, 'filled');
        
        % Add unit circle
        theta = linspace(0, 2*pi, 100);
        plot(cos(theta)*180/pi, ones(size(theta)), 'k--', 'LineWidth', 1);
        
        xlabel('\angle\Gamma_g (degrees)', 'FontSize', 10);
        ylabel('|\Gamma_g|', 'FontSize', 10);
        title(sprintf('pH = %.2f', RES(g).pH), 'FontSize', 12);
        
        % Set axis ranges
        xlim([-180, 180]);
        ylim([0, 1.1]);
        
        % Add grid
        grid on;
        
        % Add axes
        plot(xlim, [0, 0], 'k-', 'LineWidth', 0.5);
        plot([0, 0], ylim, 'k-', 'LineWidth', 0.5);
        
        % Add colorbar
        colormap(jet);
        colorbar;
        caxis([min(eta_eff), max(eta_eff)]);
    end
    
    sgtitle('Generalized reflection coefficient parameter space trajectory (|\Gamma_g| vs \angle\Gamma_g)', ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure 3
    if plotOpt.saveFigures
        saveas(fig3, fullfile(plotOpt.savePath, 'gamma_trajectories.png'));
        savefig(fig3, fullfile(plotOpt.savePath, 'gamma_trajectories.fig'));
    end
    
    %% Figure 4: Gmag and Gph statistical distributions
fig4 = figure('Name', 'Gamma Statistics', ...
              'Color', 'w', 'Position', [100, 100, 1400, 600]);

% Collect Gmag and Gph data for all pH
all_Gmag_data = cell(nG, 1);
all_Gph_data = cell(nG, 1);
pH_values = zeros(nG, 1);

for g = 1:nG
    eta_raw = RES(g).eta;
    lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
    hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
    if ~(isfinite(lo) && isfinite(hi) && hi > lo)
        lo = min(eta_raw);
        hi = max(eta_raw);
    end
    eta_grid = linspace(lo, hi, plotOpt.nPoints)';
    
    [Gmag, Gph, ~, ~] = calc_caismith_from_params_improved(...
        RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
    
    % Convert phase to degrees and wrap
    Gph_deg = Gph * 180/pi;
    Gph_deg = mod(Gph_deg + 180, 360) - 180;
    
    all_Gmag_data{g} = Gmag;
    all_Gph_data{g} = Gph_deg;
    pH_values(g) = RES(g).pH;
end

% Subplot 1: Gmag boxplot
subplot(2, 3, 1);
hold on; grid on; box on;

% Create grouped data for boxplot
box_data = [];
group_labels = [];
for g = 1:nG
    box_data = [box_data; all_Gmag_data{g}];
    group_labels = [group_labels; g * ones(length(all_Gmag_data{g}), 1)];
end

% Use grouping method to plot boxplot
boxplot(box_data, group_labels, 'Labels', cellstr(num2str(pH_values, '%.2f')));
ylabel('|\Gamma_g|', 'FontSize', 12);
title('|\Gamma_g| distribution (boxplot)', 'FontSize', 14);
ylim([0, 1.1]);

% Subplot 2: Gph boxplot
subplot(2, 3, 2);
hold on; grid on; box on;

box_data = [];
group_labels = [];
for g = 1:nG
    box_data = [box_data; all_Gph_data{g}];
    group_labels = [group_labels; g * ones(length(all_Gph_data{g}), 1)];
end

boxplot(box_data, group_labels, 'Labels', cellstr(num2str(pH_values, '%.2f')));
ylabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
title('\angle\Gamma_g distribution (boxplot)', 'FontSize', 14);

% Subplot 3: Gmag histogram (all data)
subplot(2, 3, 3);
all_Gmag = [];
for g = 1:nG
    all_Gmag = [all_Gmag; all_Gmag_data{g}];
end
histogram(all_Gmag, 30, 'FaceColor', 'blue', 'EdgeColor', 'black', 'Normalization', 'probability');
xlabel('|\Gamma_g|', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
title('|\Gamma_g| histogram (all pH)', 'FontSize', 14);
grid on;

% Subplot 4: Gph histogram (all data)
subplot(2, 3, 4);
all_Gph = [];
for g = 1:nG
    all_Gph = [all_Gph; all_Gph_data{g}];
end
histogram(all_Gph, 30, 'FaceColor', 'red', 'EdgeColor', 'black', 'Normalization', 'probability');
xlabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
ylabel('Probability density', 'FontSize', 12);
title('\angle\Gamma_g histogram (all pH)', 'FontSize', 14);
grid on;

% Subplot 5: Gmag vs pH (mean and standard deviation)
subplot(2, 3, 5);
hold on; grid on; box on;

Gmag_mean = zeros(nG, 1);
Gmag_std = zeros(nG, 1);

for g = 1:nG
    Gmag_mean(g) = mean(all_Gmag_data{g});
    Gmag_std(g) = std(all_Gmag_data{g});
end

errorbar(1:nG, Gmag_mean, Gmag_std, 'bo-', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'blue');

xlabel('pH', 'FontSize', 12);
ylabel('|\Gamma_g| mean ± std', 'FontSize', 12);
title('|\Gamma_g| statistics vs pH', 'FontSize', 14);
set(gca, 'XTick', 1:nG, 'XTickLabel', cellstr(num2str(pH_values, '%.2f')));
xtickangle(45);
ylim([0, 1.1]);

% Subplot 6: Gph vs pH (mean and standard deviation)
subplot(2, 3, 6);
hold on; grid on; box on;

Gph_mean = zeros(nG, 1);
Gph_std = zeros(nG, 1);

for g = 1:nG
    Gph_mean(g) = mean(all_Gph_data{g});
    Gph_std(g) = std(all_Gph_data{g});
end

errorbar(1:nG, Gph_mean, Gph_std, 'ro-', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'red');

xlabel('pH', 'FontSize', 12);
ylabel('\angle\Gamma_g mean ± std', 'FontSize', 12);
title('\angle\Gamma_g statistics vs pH', 'FontSize', 14);
set(gca, 'XTick', 1:nG, 'XTickLabel', cellstr(num2str(pH_values, '%.2f')));
xtickangle(45);

sgtitle('Generalized reflection coefficient statistical distribution analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    %% Figure 5: 3D visualization - Gmag, Gph, η_eff
    fig5 = figure('Name', '3D Gamma Visualization', ...
                  'Color', 'w', 'Position', [100, 100, 1200, 500]);
    
    % Subplot 1: 3D scatter plot
    subplot(1, 2, 1);
    hold on; grid on; box on;
    
    for g = 1:nG
        eta_raw = RES(g).eta;
        lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
        hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
        if ~(isfinite(lo) && isfinite(hi) && hi > lo)
            lo = min(eta_raw);
            hi = max(eta_raw);
        end
        eta_grid = linspace(lo, hi, 100)';  % Reduce points for better 3D rendering
        
        [Gmag, Gph, ~, eta_eff] = calc_caismith_from_params_improved(...
            RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
        
        % Convert phase to degrees
        Gph_deg = Gph * 180/pi;
        Gph_deg = mod(Gph_deg + 180, 360) - 180;
        
        % 3D scatter plot
        scatter3(Gph_deg, Gmag, eta_eff, 30, colors(g, :), 'filled', ...
                'DisplayName', sprintf('pH=%.2f', RES(g).pH));
    end
    
    xlabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
    ylabel('|\Gamma_g|', 'FontSize', 12);
    zlabel('\eta_{eff} (V)', 'FontSize', 12);
    title('3D parameter space: \Gamma_g vs \eta_{eff}', 'FontSize', 14);
    legend('Location', 'best');
    view(45, 30);  % Set viewpoint
    
    % Subplot 2: 3D surface plot (single pH example)
    subplot(1, 2, 2);
    hold on; grid on; box on;
    
    % Select one pH as example (e.g., first)
    g = min(1, nG);
    
    eta_raw = RES(g).eta;
    lo = prctile(eta_raw, plotOpt.etaRangePercentile(1));
    hi = prctile(eta_raw, plotOpt.etaRangePercentile(2));
    if ~(isfinite(lo) && isfinite(hi) && hi > lo)
        lo = min(eta_raw);
        hi = max(eta_raw);
    end
    
    % Create grid
    eta_grid = linspace(lo, hi, 50)';
    
    % Calculate Gmag and Gph
    [Gmag, Gph, ~, eta_eff] = calc_caismith_from_params_improved(...
        RES(g).P, eta_grid, plotOpt.useEtaEff, maxExp, newtonMaxIter, newtonTol, 'none');
    
    % Convert phase to degrees
    Gph_deg = Gph * 180/pi;
    Gph_deg = mod(Gph_deg + 180, 360) - 180;
    
    % 3D surface plot (using Gmag as height)
    [X, Z] = meshgrid(linspace(min(Gph_deg), max(Gph_deg), 20), ...
                      linspace(min(eta_eff), max(eta_eff), 20));
    
    % Interpolate Y values (Gmag)
    Y = griddata(Gph_deg, eta_eff, Gmag, X, Z, 'cubic');
    
    surf(X, Y, Z, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    
    xlabel('\angle\Gamma_g (degrees)', 'FontSize', 12);
    ylabel('|\Gamma_g|', 'FontSize', 12);
    zlabel('\eta_{eff} (V)', 'FontSize', 12);
    title(sprintf('pH=%.2f: |\Gamma_g| surface', RES(g).pH), 'FontSize', 14);
    view(45, 30);
    colormap(jet);
    colorbar;
    
    sgtitle('Generalized reflection coefficient 3D visualization', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure 5
    if plotOpt.saveFigures
        saveas(fig5, fullfile(plotOpt.savePath, 'gamma_3d.png'));
        savefig(fig5, fullfile(plotOpt.savePath, 'gamma_3d.fig'));
    end
    
    %% Output statistical summary
    fprintf('\n=== Generalized reflection coefficient statistical summary ===\n');
    fprintf('%-10s %-15s %-15s %-15s %-15s\n', ...
            'pH', 'Gmag mean', 'Gmag std', 'Gph mean(°)', 'Gph std(°)');
    fprintf('%s\n', repmat('-', 70, 1));
    
    for g = 1:nG
        fprintf('%-10.2f %-15.4f %-15.4f %-15.4f %-15.4f\n', ...
                RES(g).pH, Gmag_mean(g), Gmag_std(g), Gph_mean(g), Gph_std(g));
    end
    
    fprintf('\nAll pH combined statistics:\n');
    fprintf('  Gmag: mean = %.4f, std = %.4f, range = [%.4f, %.4f]\n', ...
            mean(all_Gmag), std(all_Gmag), min(all_Gmag), max(all_Gmag));
    fprintf('  Gph:  mean = %.2f°, std = %.2f°, range = [%.2f°, %.2f°]\n', ...
            mean(all_Gph), std(all_Gph), min(all_Gph), max(all_Gph));
    
    fprintf('\nVisualization complete! Figures saved to: %s\n', plotOpt.savePath);
end

function P = ensure_modeA_primary(P)
    % Ensure mode A jstar >= mode B jstar
    % If jstarA < jstarB, swap all parameters of the two modes
    
    if ~isfield(P, 'jstarA') || ~isfield(P, 'jstarB')
        return;  % If not dual mode, return directly
    end
    
    % Check if swap needed
    if P.jstarA >= P.jstarB
        return;  % Already mode A primary, no swap
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

function check_mode_consistency(RES)
    % Check mode parameter consistency across all pH conditions
    fprintf('\n=== Mode consistency check ===\n');
    
    % Get nG from RES
    nG = numel(RES);
    
    all_swapped = false(nG, 1);
    primary_mode_ratio = zeros(nG, 1);
    
    for g = 1:nG
        P = RES(g).P;
        
        if isfield(P, 'modesSwapped')
            all_swapped(g) = P.modesSwapped;
        end
        
        % Calculate primary/secondary mode current ratio
        if isfield(P, 'jstarA') && isfield(P, 'jstarB')
            primary_mode_ratio(g) = P.jstarA / max(P.jstarB, eps);
            fprintf('pH=%.3g: jstarA=%.2e, jstarB=%.2e, ratio=%.2f', ...
                RES(g).pH, P.jstarA, P.jstarB, primary_mode_ratio(g));
            
            if all_swapped(g)
                fprintf(' (modes swapped)\n');
            else
                fprintf(' (modes not swapped)\n');
            end
        end
    end
    
    % Statistics
    swapped_count = sum(all_swapped);
    fprintf('\nStatistics: %d/%d pH conditions had modes swapped\n', swapped_count, nG);
    fprintf('Primary/secondary mode average ratio: %.2f\n', mean(primary_mode_ratio));
    
    % Visualize mode swap situation
    figure('Color', 'w', 'Name', 'Mode swap check', 'Position', [100, 100, 800, 400]);
    
    subplot(1, 2, 1);
    bar([sum(~all_swapped), swapped_count]);
    set(gca, 'XTickLabel', {'Not swapped', 'Swapped'});
    ylabel('Count');
    title('Mode swap statistics');
    grid on;
    
    subplot(1, 2, 2);
    pH_vals = [RES.pH];
    plot(pH_vals, primary_mode_ratio, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(pH_vals(all_swapped), primary_mode_ratio(all_swapped), 'ro', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('pH');
    ylabel('jstarA/jstarB ratio');
    title('Primary/secondary mode ratio');
    grid on;
    legend({'All points', 'Swapped points'}, 'Location', 'best');
    
    sgtitle('Mode parameter consistency analysis');
end

function analyze_waveguide_relativity(phi, beta_w, I_store, eta_trans, eta_eff)
    % Analyze waveguide relativistic parameters
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. Rapidity vs overpotential
    subplot(2, 3, 1);
    plot(eta_eff, phi, 'b-', 'LineWidth', 2);
    xlabel('Effective overpotential \eta_{eff} (V)');
    ylabel('Rapidity \phi');
    title('Rapidity distribution');
    grid on;
    
    % Mark characteristic regions
    hold on;
    plot([min(eta_eff), max(eta_eff)], [0, 0], 'k--');
    text(mean(eta_eff), 0.5, '\phi>0: forward dominant', 'Color', 'r');
    text(mean(eta_eff), -0.5, '\phi<0: backward dominant', 'Color', 'b');
    
    % 2. Waveguide velocity vs overpotential
    subplot(2, 3, 2);
    plot(eta_eff, beta_w, 'r-', 'LineWidth', 2);
    xlabel('Effective overpotential \eta_{eff} (V)');
    ylabel('Waveguide velocity \beta_w');
    title('Waveguide velocity distribution');
    ylim([-1.1, 1.1]);
    grid on;
    
    % 3. Storage intensity vs overpotential
    subplot(2, 3, 3);
    plot(eta_eff, I_store, 'g-', 'LineWidth', 2);
    xlabel('Effective overpotential \eta_{eff} (V)');
    ylabel('Normalized storage intensity I_{store}');
    title('Storage intensity distribution');
    ylim([0, 1.1]);
    grid on;
    
    % 4. Transmission efficiency vs overpotential
    subplot(2, 3, 4);
    plot(eta_eff, eta_trans, 'm-', 'LineWidth', 2);
    xlabel('Effective overpotential \eta_{eff} (V)');
    ylabel('Transmission efficiency \eta_{trans}');
    title('Transmission efficiency distribution');
    ylim([0, 1.1]);
    grid on;
    
    % 5. Phase plot: storage intensity vs transmission efficiency
    subplot(2, 3, 5);
    scatter(I_store, eta_trans, 30, eta_eff, 'filled');
    xlabel('Storage intensity I_{store}');
    ylabel('Transmission efficiency \eta_{trans}');
    title('Storage-transmission phase plot');
    colorbar;
    grid on;
    
    % 6. Rapidity-storage intensity relationship
    subplot(2, 3, 6);
    scatter(phi, I_store, 30, eta_eff, 'filled');
    xlabel('Rapidity \phi');
    ylabel('Storage intensity I_{store}');
    title('Rapidity-storage intensity relationship');
    colorbar;
    grid on;
    
    % Calculate statistics
    fprintf('===== Waveguide relativistic parameter statistics =====\n');
    fprintf('Average storage intensity: %.4f\n', mean(I_store));
    fprintf('Average transmission efficiency: %.4f\n', mean(eta_trans));
    fprintf('Maximum storage intensity: %.4f at η_eff = %.4f V\n', ...
        max(I_store), eta_eff(I_store == max(I_store)));
    fprintf('Maximum transmission efficiency: %.4f at η_eff = %.4f V\n', ...
        max(eta_trans), eta_eff(eta_trans == max(eta_trans)));
    
    % Identify operating modes
    high_storage = I_store > 0.8;
    high_transport = eta_trans > 0.8;
    if any(high_storage)
        fprintf('Strong storage mode region: η_eff ∈ [%.4f, %.4f] V\n', ...
            min(eta_eff(high_storage)), max(eta_eff(high_storage)));
    end
    if any(high_transport)
        fprintf('Strong transport mode region: η_eff ∈ [%.4f, %.4f] V\n', ...
            min(eta_eff(high_transport)), max(eta_eff(high_transport)));
    end
end

function visualize_waveguide_relativity_comprehensive(Gmag, Gph, I_store, eta_eff, ...
                                                      U, S_rel, beta_w, phi, ...
                                                      D_plus, D_minus, eta_trans)
    % Comprehensive visualization of waveguide relativistic parameters
    % Input parameters from calc_caismith_from_params_improved function
    
    % Set up figure window
    figure('Position', [50, 50, 1600, 1000], 'Name', 'Waveguide relativistic parameter visualization', ...
           'NumberTitle', 'off');
    
    %% Subplot 1: Trajectory on Cai-Smith unit disk
    subplot(3, 4, [1, 2, 5, 6]);
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k-', 'LineWidth', 1.5); % Unit circle
    hold on;
    
    % Color map represents overpotential
    cmap = jet(256);
    color_indices = floor(interp1([min(eta_eff), max(eta_eff)], [1, 256], eta_eff));
    color_indices = max(1, min(256, color_indices));
    
    % Plot trajectory
    for i = 1:length(Gmag)-1
        x1 = Gmag(i) * cos(Gph(i));
        y1 = Gmag(i) * sin(Gph(i));
        x2 = Gmag(i+1) * cos(Gph(i+1));
        y2 = Gmag(i+1) * sin(Gph(i+1));
        
        % Use overpotential color
        color1 = cmap(color_indices(i), :);
        color2 = cmap(color_indices(i+1), :);
        
        % Plot line segment
        plot([x1, x2], [y1, y2], 'Color', color1, 'LineWidth', 1.5);
        
        % Mark start and end points
        if i == 1
            plot(x1, y1, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
            text(x1, y1, 'Start', 'VerticalAlignment', 'bottom', ...
                 'HorizontalAlignment', 'right', 'FontSize', 10);
        elseif i == length(Gmag)-1
            plot(x2, y2, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            text(x2, y2, 'End', 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left', 'FontSize', 10);
        end
    end
    
    % Mark important points
    [~, idx_max_store] = max(I_store);
    x_store = Gmag(idx_max_store) * cos(Gph(idx_max_store));
    y_store = Gmag(idx_max_store) * sin(Gph(idx_max_store));
    plot(x_store, y_store, 'ms', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    text(x_store, y_store, '  Max storage', 'FontSize', 10);
    
    [~, idx_max_trans] = max(eta_trans);
    x_trans = Gmag(idx_max_trans) * cos(Gph(idx_max_trans));
    y_trans = Gmag(idx_max_trans) * sin(Gph(idx_max_trans));
    plot(x_trans, y_trans, 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    text(x_trans, y_trans, '  Max transmission', 'FontSize', 10);
    
    % Mark disk center
    plot(0, 0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
    text(0, 0, ' Γ_g=0', 'VerticalAlignment', 'top', 'FontSize', 10);
    
    % Set axes
    axis equal;
    axis([-1.2, 1.2, -1.2, 1.2]);
    xlabel('Re(Γ_g)');
    ylabel('Im(Γ_g)');
    title('Cai-Smith unit disk trajectory (color represents overpotential)');
    grid on;
    
    % Add colorbar
    colormap(cmap);
    c = colorbar;
    c.Label.String = 'η_{eff} (V)';
    
    %% Subplot 2: Rapidity and waveguide velocity
    subplot(3, 4, 3);
    yyaxis left;
    plot(eta_eff, phi, 'b-', 'LineWidth', 2);
    ylabel('Rapidity φ');
    ylim([min(phi)-0.5, max(phi)+0.5]);
    
    yyaxis right;
    plot(eta_eff, beta_w, 'r-', 'LineWidth', 2);
    ylabel('Waveguide velocity β_w');
    ylim([-1.1, 1.1]);
    
    xlabel('Effective overpotential η_{eff} (V)');
    title('Rapidity and waveguide velocity');
    grid on;
    legend('φ (rapidity)', 'β_w (waveguide velocity)', 'Location', 'best');
    
    %% Subplot 3: Doppler factors
    subplot(3, 4, 4);
    semilogy(eta_eff, D_plus, 'g-', 'LineWidth', 2, 'DisplayName', 'D_+');
    hold on;
    semilogy(eta_eff, D_minus, 'm-', 'LineWidth', 2, 'DisplayName', 'D_-');
    semilogy(eta_eff, D_plus .* D_minus, 'k--', 'LineWidth', 1, 'DisplayName', 'D_+·D_-');
    
    xlabel('Effective overpotential η_{eff} (V)');
    ylabel('Doppler factor (log scale)');
    title('Doppler factors');
    grid on;
    legend('Location', 'best');
    
    %% Subplot 4: Storage intensity and transmission efficiency
    subplot(3, 4, 7);
    yyaxis left;
    plot(eta_eff, I_store, 'c-', 'LineWidth', 2);
    ylabel('Storage intensity I_{store}');
    ylim([0, 1.1]);
    
    yyaxis right;
    plot(eta_eff, eta_trans, 'm-', 'LineWidth', 2);
    ylabel('Transmission efficiency η_{trans}');
    ylim([0, 1.1]);
    
    xlabel('Effective overpotential η_{eff} (V)');
    title('Storage intensity vs transmission efficiency');
    grid on;
    legend('I_{store}', 'η_{trans}', 'Location', 'best');
    
    %% Subplot 5: Energy-type U and power-flow-type S
    subplot(3, 4, 8);
    yyaxis left;
    semilogy(eta_eff, U, 'b-', 'LineWidth', 2);
    ylabel('Energy-type U (log)');
    
    yyaxis right;
  semilogy(eta_eff, S_rel, 'r-', 'LineWidth', 2);
    ylabel('Power-flow-type S');
    
    xlabel('Effective overpotential η_{eff} (V)');
    title('U and S parameters');
    grid on;
    legend('U', 'S', 'Location', 'best');
    
    %% Subplot 6: Storage-transmission phase plot
    subplot(3, 4, 9);
    scatter(I_store, eta_trans, 40, eta_eff, 'filled');
    xlabel('Storage intensity I_{store}');
    ylabel('Transmission efficiency η_{trans}');
    title('Storage-transmission phase plot');
    xlim([0, 1.1]);
    ylim([0, 1.1]);
    grid on;
    colormap(cmap);
    c2 = colorbar;
    c2.Label.String = 'η_{eff} (V)';
    
    % Add operating region markers
    hold on;
    % Draw storage mode rectangle
x1 = [0.7, 1.0, 1.0, 0.7, 0.7];
y1 = [0, 0, 0.3, 0.3, 0];
plot(x1, y1, 'g--', 'LineWidth', 2, 'DisplayName', 'Storage mode');

% Draw transmission mode rectangle
x2 = [0, 0.3, 0.3, 0, 0];
y2 = [0.7, 0.7, 1.0, 1.0, 0.7];
plot(x2, y2, 'b--', 'LineWidth', 2, 'DisplayName', 'Transmission mode');

% Draw mixed mode rectangle
x3 = [0.5, 0.7, 0.7, 0.5, 0.5];
y3 = [0.5, 0.5, 0.7, 0.7, 0.5];
plot(x3, y3, 'r--', 'LineWidth', 2, 'DisplayName', 'Mixed mode');

% Adjust legend to ensure we don't add duplicate entries for other curves
% Since we only plot these three rectangles in this subplot, we can use legend directly
legend('Location', 'southeast');
    
    %% Subplot 7: Rapidity vs storage intensity relationship
    subplot(3, 4, 10);
    scatter(phi, I_store, 40, eta_eff, 'filled');
    xlabel('Rapidity φ');
    ylabel('Storage intensity I_{store}');
    title('Rapidity vs storage intensity');
    grid on;
    colormap(cmap);
    c3 = colorbar;
    c3.Label.String = 'η_{eff} (V)';
    
    %% Subplot 8: Waveguide velocity vs transmission efficiency relationship
    subplot(3, 4, 11);
    plot(beta_w, eta_trans, 'b-', 'LineWidth', 2);
    xlabel('Waveguide velocity β_w');
    ylabel('Transmission efficiency η_{trans}');
    title('β_w vs η_{trans}');
    xlim([-1.1, 1.1]);
    ylim([0, 1.1]);
    grid on;
    
    %% Subplot 9: Radial evolution (|Γ_g| vs η_eff)
    subplot(3, 4, 12);
    plot(eta_eff, Gmag, 'k-', 'LineWidth', 2);
    xlabel('Effective overpotential η_{eff} (V)');
    ylabel('|Γ_g| (generalized reflection magnitude)');
    title('Radial evolution');
    ylim([0, 1.1]);
    grid on;
    
    % Add important point markers
    hold on;
    plot(eta_eff(idx_max_store), Gmag(idx_max_store), 'mo', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'm');
    plot(eta_eff(idx_max_trans), Gmag(idx_max_trans), 'bo', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    %% Add statistics text
    annotation('textbox', [0.02, 0.02, 0.96, 0.05], ...
               'String', sprintf(['Statistics: Average storage intensity = %.3f, Average transmission efficiency = %.3f, ', ...
                                  'Max storage at η=%.3fV (I_store=%.3f), ', ...
                                  'Max transmission at η=%.3fV (η_trans=%.3f)'], ...
                                 mean(I_store), mean(eta_trans), ...
                                 eta_eff(idx_max_store), max(I_store), ...
                                 eta_eff(idx_max_trans), max(eta_trans)), ...
               'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'left', ...
               'BackgroundColor', [0.95, 0.95, 0.95]);
end

function plot_waveguide_electrochem_summary_figure(RES, maxExp, newtonMaxIter, newtonTol)
% Summary figure (2x2) with optimum defined by MAX useful-output density.
%   Pi_use(eta_eff)  = max(0,-jTot) * (1 - rho^2)
%   Pi_dens(eta_eff) = Pi_use / (|eta_eff| + eta0)
% Optimum in panel (c): argmax Pi_dens on cathodic side (eta_eff<0).
%
% Requires:
%   safe_two_mode_model, calc_caismith_from_params_improved, prctile_fallback.

    if nargin < 2 || isempty(maxExp), maxExp = 70; end
    if nargin < 3 || isempty(newtonMaxIter), newtonMaxIter = 10; end
    if nargin < 4 || isempty(newtonTol), newtonTol = 5e-11; end

    nG = numel(RES);
    C  = lines(nG);

    % ---------- options ----------
    opt = struct();
    opt.useEtaEff = true;
    opt.nEta = 500;
    opt.etaPctRange = [5 95];
    opt.useFixedWindow = true;
    opt.fixedEtaWindow = [-1.5 0.5];
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
        % etaGrid = linspace(-1.5, 0.5, opt.nEta).';

        [Gmag, Gph, ~, eta_eff_map, U, S_rel, beta_w] = ...
            calc_caismith_from_params_improved(RES(g).P, etaGrid, opt.useEtaEff, ...
                                               maxExp, newtonMaxIter, newtonTol, 'none'); %#ok<ASGLU>

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
    lg1 = legend(ax1, 'Location','east', 'FontSize', 12); %#ok<NASGU>
    legend(ax1,'boxoff');

    % ==================== (b) rho and chi vs eta_eff ====================
ax2 = nexttile(tlo,2);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
set(ax2,'FontSize',10,'Layer','top','Clipping','on');

% ---- Left axis: rho ----
yyaxis(ax2,'left');
for g = 1:nG
    d = storage_data{g};
    plot(ax2, d.eta_eff_map, d.rho, '-', ...
        'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
end
ylabel(ax2, '$\rho=|\Gamma_g|$', 'Interpreter','latex', 'FontSize',11, 'Color','k');
ylim(ax2, [0, 1.05]);
ax2.YColor = 'k';

% x-axis (shared by both sides)
    xlabel(ax2, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex','FontSize',12);
if opt.useFixedWindow, xlim(ax2, opt.fixedEtaWindow); end

% ---- Right axis: chi ----
yyaxis(ax2,'right');
for g = 1:nG
    d = storage_data{g};
    chi = max(min(d.chi,1),-1);
    plot(ax2, d.eta_eff_map, chi, '--', ...
        'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');
end
ylabel(ax2, '$\beta_{w}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.85 0.33 0.10]);
ylim(ax2, [-1.05, 1.05]);
ax2.YColor = [0.85 0.33 0.10];

% ---- Legend: Use dummy lines, only show line style meaning ----
yyaxis(ax2,'left');  hR = plot(ax2, NaN, NaN, 'k-',  'LineWidth',2.0, 'DisplayName', '$\rho=|\Gamma_g|$');
yyaxis(ax2,'right'); hC = plot(ax2, NaN, NaN, 'k--', 'LineWidth',2.0, 'DisplayName', '$\beta_{w}$');
legend1=legend(ax2, [hR hC], 'Location','best', 'Interpreter','latex', 'FontSize',12, 'Box','off');
set(legend1,...
    'Position',[0.621377893541425 0.717865915229878 0.114346409217984 0.0559640529108982],...
    'Interpreter','latex',...
    'FontSize',12);




    % ==================== (c) useful output + density optimum ====================
    ax3 = nexttile(tlo,3);
    hold(ax3,'on'); box(ax3,'on'); 
    grid(ax3,'on');

    eta0 = opt.eta0;
    eps0 = 1e-30;

    for g = 1:nG
        P = RES(g).P;
        d = storage_data{g};
        etaGrid = d.etaGrid_raw(:);

        jTot = safe_two_mode_model(P, etaGrid, maxExp, newtonMaxIter, newtonTol, []);  % mA/cm^2
        eta_eff = etaGrid - P.Rohm.*(jTot./1000);                                      % V
        rho = d.rho(:);

        n = min([numel(eta_eff), numel(rho), numel(jTot)]);
        eta_eff = eta_eff(1:n); rho = rho(1:n); jTot = jTot(1:n);

        Puse = max(0, -jTot) .* (1 - rho.^2);     % mA/cm^2
        Puse =abs(jTot) .* (1 - rho.^2);     % mA/cm^2
        Pdens = Puse ./ (abs(eta_eff) + eta0);    % mA/cm^2 / V  (density)

        if opt.cathodicOnly
            mask = (eta_eff < 0) & isfinite(eta_eff) & isfinite(Puse) & isfinite(Pdens);
        else
            mask = isfinite(eta_eff) & isfinite(Puse) & isfinite(Pdens);
        end
        if nnz(mask) < 8, continue; end

        x = eta_eff(mask);
        y = Puse(mask);
        z = Pdens(mask);

        [x, ord] = sort(x);
        y = y(ord);
        z = z(ord);

        % plot useful output curve
        plot(ax3, x, z, '-', 'Color', C(g,:), 'LineWidth', 2.0, 'HandleVisibility','off');

        % density optimum
        [zmax, im] = max(z);
        xstar = x(im);
        ystar = z(im);

        % plot(ax3, xstar, ystar, 'o', 'MarkerSize', 7, ...
        %     'MarkerFaceColor', C(g,:), 'MarkerEdgeColor','k', 'LineWidth', 1.0, ...
        %     'DisplayName', sprintf('pH=%.1f: $\\eta^*_{\\mathrm{eff}}=%.2f$ V', RES(g).pH, xstar));
                plot(ax3, xstar, ystar, 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', C(g,:), 'MarkerEdgeColor','k', 'LineWidth', 1.0, ...
            'DisplayName', sprintf('$\\eta^*_{\\mathrm{eff}}=%.2f$ V', xstar));
        fprintf('pH=%.2f: density optimum (eta0=%.2f): eta_eff*=%.3f V, Pi_use*=%.3g mA/cm^2, Pi_dens*=%.3g\n', ...
            RES(g).pH, eta0, xstar, ystar, zmax);
    end
    set(ax3,'YScale','log');
    ax3.YAxis.MinorTick = 'off';  % Add this line to hide minor y-ticks
    set(ax3, 'MinorGridLineStyle', 'none');
    % set(ax3, 'YMinorGrid', 'off');
    xlabel(ax3, '$\eta_{\mathrm{eff}}$ (V)', 'Interpreter','latex');
    ylabel(ax3, '$\Pi_{\mathrm{dens}}$ (mA cm$^{-2}$ V$^{-1}$)', ...
    'Interpreter','latex');

    % title(ax3, sprintf('Useful output vs $\\eta_{\\mathrm{eff}}$ (circles: max $\\Pi_{\\rm dens}$, $\\eta_0$=%.2f V)', eta0), 'Interpreter','latex');
    if opt.useFixedWindow, xlim(ax3, opt.fixedEtaWindow); end
    lg3 = legend(ax3, 'Location','southwest', 'FontSize', 12, 'Interpreter','latex'); %#ok<NASGU>
    set(lg3,'NumColumns',2);
    legend(ax3,'boxoff');
    set(lg3,...
    'Position',[0.0858296896907499 0.0859403890721938 0.299076816334444 0.158660140692019],...
    'NumColumns',2,...
    'Interpreter','latex');

    % ==================== (d) log10(U) and log10(|S|) ====================
    ax4 = nexttile(tlo,4);
    hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');

    for g = 1:nG
        d = storage_data{g};
        x = d.eta_eff_map;
        % yU = log10(max(d.U, eps0));
        % yS = log10(max(abs(d.S), eps0));
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
    lg4 = legend(ax4, [hU hS], 'Location','southwest', 'FontSize', 12, 'Interpreter','latex'); %#ok<NASGU>
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
% filename = 'Figure04';
% set(gcf,'Units','centimeters');                      % Ensure consistent units
% set(gcf,'PaperPositionMode','auto');                 % Key: paper = screen
% exportgraphics(gcf, [filename,'.png'], 'Resolution', 600);
filename = 'Figure04';
set(gcf,'Units','centimeters');                      % Ensure consistent units
set(gcf,'PaperPositionMode','auto');                 % Key: paper = screen
saveas(gcf, [filename,'.png']);

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
%
% Input:
%   RES(g).pH, RES(g).eta (raw eta), RES(g).P (fit params)
%   P must contain: jstarA, xiA, lamA, jlimA, jstarB, xiB, lamB, jlimB, Rohm
%   (aA,bA,aB,bB optional; will be auto-built if absent and f is provided)
%
% Options (opt):
%   opt.useEtaEff      = true
%   opt.nEta           = 700
%   opt.useFixedWindow = true
%   opt.fixedEtaWindow = [-1.5 0.5]
%   opt.etaPctRange    = [5 95]      (used only if useFixedWindow=false)
%   opt.cathodicOnly   = true
%   opt.eta0           = 0.10        (regularization in Pi_dens; MUST be >0 for stability)
%   opt.zoomToData     = true
%   opt.zoomMargin     = 0.12
%   opt.showUnitCircle = true
%   opt.showGrid       = false
%   opt.maxExp         = 70
%   opt.newtonMaxIter  = 10
%   opt.newtonTol      = 5e-11
%   opt.pointSize      = 14
%   opt.savePng        = ''          (e.g., 'CS_normdens.png')
%   opt.colormapName   = 'turbo'     ('turbo'|'parula'|...)
%
% Requires:
%   safe_two_mode_model, prctile_fallback, clip (you already have these in your file)
%
% ------------------------------------------------------------

if nargin < 2, opt = struct(); end
opt = setdef(opt,'useEtaEff',true);
opt = setdef(opt,'nEta',700);
opt = setdef(opt,'etaPctRange',[5 95]);
opt = setdef(opt,'useFixedWindow',true);
opt = setdef(opt,'fixedEtaWindow',[-1.5 0.5]);
opt = setdef(opt,'cathodicOnly',true);
opt = setdef(opt,'eta0',0.10);
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
% --- normalize by per-pH maximum ---
r = Pi_dens ./ (Pi_dens_max + 1e-300);   % in [0,1] on mask
r(~mask) = NaN;
r = min(max(r,0),1);

% --- dynamic-range compression (choose ONE) ---
% (A) log-compressed to 0..1 (best for long-tail)
epsr = 1e-4;   % between 1e-3~1e-5; smaller "brightens" low values
Cnorm = (log10(r + epsr) - log10(epsr)) / (log10(1 + epsr) - log10(epsr));

% (B) alternatively gamma (simpler)
% gamma = 0.25;   % recommended 0.2~0.5
% Cnorm = r.^gamma;

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
nShow = min(nG,7);

ncol = 4;
nrow = 2; % fixed 2x4 for publication
width_cm = 24;
height_cm = 18;

fig = figure('Color','w','Units','centimeters','Position',[2 2 width_cm height_cm]);
tlo = tiledlayout(fig, nrow, ncol, 'Padding','compact', 'TileSpacing','compact');

% colormap
apply_colormap(fig, opt.colormapName);

% title
% title(tlo, 'Cai--Smith disks colored by normalized output density $\Pi_{\mathrm{dens}}/\max(\Pi_{\mathrm{dens}})$', ...
%     'Interpreter','latex', 'FontSize', 14);

axisCol = [0.65 0.65 0.65];
trajCol = [0.80 0.80 0.80];

for g = 1:nShow
    ax = nexttile(tlo, g);
    hold(ax,'on'); box(ax,'on'); %axis(ax,'equal');

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

        txt = sprintf('$\\eta^*_{\\rm eff}=%.2f\\,\\mathrm{V}$\n$\\Pi_{\\rm dens}^{\\max}=%.2g$', etaStar, PmaxAbs);
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
    % title(ax, sprintf('pH=%.2g', OUT(g).pH), 'Interpreter','none', 'FontSize', 12);
            % Add label (a, b, c, ...) to each subplot
    text(-0.25, 0.98, char('a' + g - 1), 'Units','normalized', ...
            'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
            'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
            'Color','k');
end

% tile 8: EMPTY (keep layout)
ax8 = nexttile(tlo, 8);
axis(ax8,'off');

% -------- single global colorbar on the far right --------
% Use an invisible axes in figure coordinates to host the colorbar
axCB = axes(fig, 'Position',[0.80 0.10 0.02 0.38], 'Visible','off'); %#ok<LAXES>
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
print(fig,'-dpng','-r600','FigureS16.png');
fprintf('Saved: FigureS16.png\n');
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

function LAW = verify_four_laws_electrochem_polarization(RES, opt)
% verify_four_laws_electrochem_polarization_general2port
% ------------------------------------------------------------
% General 2-port (U != V) waveguide-law verification from electrochemical polarization mapping.
%
% Build a general passive 2-port at each pH and evaluate four laws at an operating point eta_eff*:
%   S(eta) = U(eta) * diag(GammaA(eta), GammaB(eta)) * V(eta)^H
% where GammaA/B are complex "directional-reflection" factors extracted from UNSATURATED fluxes.
%
% Laws (at eta_eff*):
%   Law1: eigenchannel equality  alpha_Mp = epsilon_Mp = 1 - sigma_p^2
%   Law2: equal-power-profile (matched eigenchannel weights) => alpha(i)=epsilon(o)
%   Law3: trace invariance under basis rotation
%   Law4: reciprocity pairing  alpha(i) = epsilon(i*)  (generally NOT zero if U != V / nonreciprocal)
%
% Output:
%   LAW.table + LAW.metrics + Figure saved (if opt.savePng not empty)
%
% Requires (already in your project):
%   safe_two_mode_model, calc_caismith_from_params_improved
%
% ------------------------------------------------------------

fprintf('\n=== Verify four waveguide laws (electrochem general 2-port) ===\n');

if nargin < 2, opt = struct(); end
opt = setdef(opt,'maxExp',70);
opt = setdef(opt,'newtonMaxIter',10);
opt = setdef(opt,'newtonTol',5e-11);

opt = setdef(opt,'useEtaEff',true);
opt = setdef(opt,'fixedEtaWindow',[-1.5 0.5]);
opt = setdef(opt,'nEta',700);
opt = setdef(opt,'eta0',0.10);                 % for Pi_dens
opt = setdef(opt,'cathodicOnly',true);
opt = setdef(opt,'evalAt','Pi_dens_max');      % 'Pi_dens_max' | 'I_store_max' | 'eta0' (closest to 0-)
opt = setdef(opt,'tolPass',1e-9);

% Monte-Carlo sizes (keep modest for speed; increase if needed)
opt = setdef(opt,'nPairs',500);   % law2
opt = setdef(opt,'nBases',150);   % law3
opt = setdef(opt,'nTest',600);    % law4

% mixing design strength (tune if you want more/less nonreciprocity)
opt = setdef(opt,'mixGain_w',1.0);        % how strongly mode-weight imbalance controls mixing
opt = setdef(opt,'mixGain_beta',0.35);    % how strongly beta_w offsets U vs V (break reciprocity)
opt = setdef(opt,'phaseGain',0.8);        % how strongly phase-diff controls mixing phase

opt = setdef(opt,'savePng','FigureS_laws_echem_general2port.png');

nG = numel(RES);
pHs = arrayfun(@(s) s.pH, RES(:)).';

L1 = nan(nG,1); L2 = nan(nG,1); L3 = nan(nG,1); L4 = nan(nG,1);
passViol = nan(nG,1);   % max(0, max(sigma)-1) at eta*
passBadCount = nan(nG,1);

% store example scatter for a non-empty panel (worst mismatch)
alpha_i_all = cell(nG,1);
eps_o_all   = cell(nG,1);
alpha_t_all = cell(nG,1);
eps_t_all   = cell(nG,1);
ex_etaStar  = nan(nG,1);
ex_tag      = strings(nG,1);

for g = 1:nG
    P = ensure_ab_fields(RES(g).P);

    % ---- eta grid ----
    etaGrid = linspace(opt.fixedEtaWindow(1), opt.fixedEtaWindow(2), opt.nEta).';
    % ---- mapping diagnostics (global; for eta_eff and beta_w etc.) ----
    [rhoG, thetaG, I_storeG, eta_effG, Uproxy, Sproxy, beta_w] = ...
        calc_caismith_from_params_improved(P, etaGrid, opt.useEtaEff, opt.maxExp, ...
                                           opt.newtonMaxIter, opt.newtonTol, 'none');

    eta_effG = eta_effG(:);
    beta_w = beta_w(:);

    % ---- compute UNSATURATED directional fluxes per mode (needed for GammaA/B + weights) ----
    JoxA  = P.jstarA .* exp( clip_local(P.aA.*eta_effG, opt.maxExp) );
    JredA = P.jstarA .* exp( clip_local(P.bA.*eta_effG, opt.maxExp) );
    JoxB  = P.jstarB .* exp( clip_local(P.aB.*eta_effG, opt.maxExp) );
    JredB = P.jstarB .* exp( clip_local(P.bB.*eta_effG, opt.maxExp) );

    eps0 = 1e-30;
    UA = JoxA + JredA;
    UB = JoxB + JredB;
    w = (UA - UB) ./ (UA + UB + eps0);          % mode-energy imbalance in [-1,1]
    w = max(min(w,1),-1);

    % ---- complex Gamma per mode (passive by construction) ----
    [GamA, rhoA, thA] = gamma_from_flux(JoxA, JredA);
    [GamB, rhoB, thB] = gamma_from_flux(JoxB, JredB);

    % ---- choose operating point eta_eff* ----
    % compute Pi_dens on this same grid using fitted jTot (implicit)
    jTot = safe_two_mode_model(P, etaGrid, opt.maxExp, opt.newtonMaxIter, opt.newtonTol, []);
    if opt.useEtaEff
        eta_eff2 = etaGrid - P.Rohm.*(jTot./1000);
    else
        eta_eff2 = etaGrid;
    end
    % interpolate (GamA/B are on eta_effG; map to eta_eff2 for consistent argmax)
    [eta_effG_s, ordG] = sort(eta_effG);
    GamA_s = GamA(ordG); GamB_s = GamB(ordG);
    rhoA_s = rhoA(ordG); rhoB_s = rhoB(ordG);
    w_s    = w(ordG);
    beta_s = beta_w(ordG);

    [eta_eff2s, ord2] = sort(eta_eff2(:));
    jTot_s = jTot(ord2);

    % use global rho (from sum flux) for Pi_use absorption weighting; but ok to use mode-combined rho
    rhoG_s = interp1(eta_effG_s, clamp01(rhoG(ordG)), eta_eff2s, 'linear', 'extrap');
    rhoG_s = clamp01(rhoG_s);

    Pi_use  = max(0, -jTot_s) .* max(1 - rhoG_s.^2, 0);
    Pi_dens = Pi_use ./ (abs(eta_eff2s) + opt.eta0);

    maskC = isfinite(Pi_dens) & isfinite(eta_eff2s);
    if opt.cathodicOnly
        maskC = maskC & (eta_eff2s < 0);
    end

    if ~any(maskC)
        fprintf('pH=%.3g: no valid cathodic points for eta*.\n', RES(g).pH);
        continue;
    end

    switch lower(opt.evalAt)
        case 'pi_dens_max'
            tmp = Pi_dens; tmp(~maskC) = -Inf;
            [~, im] = max(tmp);
            etaStar = eta_eff2s(im);
        case 'i_store_max'
            % pick max I_store on cathodic side (from mapping grid)
            tmp = I_storeG(:);
            tmp(~(eta_effG<0)) = -Inf;
            [~, im0] = max(tmp);
            etaStar = eta_effG(im0);
        otherwise
            % closest to 0- (cathodic)
            neg = eta_eff2s(maskC);
            [~,ii] = min(abs(neg)); etaStar = neg(ii);
    end
    ex_etaStar(g) = etaStar;

    % find nearest index on eta_effG grid
    [~, kStar] = min(abs(eta_effG_s - etaStar));
    ga = GamA_s(kStar);
    gb = GamB_s(kStar);

    % guard passivity: enforce |Gamma|<=1 (numerical)
    ga = clamp_mag(ga, 1-1e-12);
    gb = clamp_mag(gb, 1-1e-12);

    % ---- build GENERAL 2-port unitary U,V from electrochem-driven mixing ----
    ww   = w_s(kStar);
    bet  = beta_s(kStar);

    % mixing angles in [0, pi/2]
    aU = 0.25*pi * (1 - opt.mixGain_w*ww);   % ww=+1 => small mix, ww=-1 => strong mix
    aU = min(max(aU, 0), 0.5*pi);
    % break reciprocity by offsetting V mixing via beta_w
    aV = aU + opt.mixGain_beta * 0.25*pi * tanh(2*bet);
    aV = min(max(aV, 0), 0.5*pi);

    % mixing phases from mode phase difference (adds structure)
    dphi = wrapToPi_local(angle(ga) - angle(gb));
    bU = opt.phaseGain * dphi;
    bV = -opt.phaseGain * dphi + 0.6*wrapToPi_local(angle(ga)+angle(gb));

    % build unitary
    U = unitary2x2(aU, bU);
    V = unitary2x2(aV, bV);

    % diagonal "eigenchannel" reflection factors (complex)
    D = diag([ga, gb]);

    % general 2-port scattering
    S = U * D * (V');  % V' = V^H in MATLAB for complex

    % ---- passivity check (at eta*) ----
    sig = svd(S);
    viol = max(0, max(sig) - 1);
    passViol(g) = viol;

    % also count passivity warnings across the grid (optional, diagnostic)
    % Here, because U,V are unitary and |Gamma|<=1, should be 0 unless numerical issues.
    passBadCount(g) = 0;

    % ---- build Ab, Em ----
    Ab = eye(2) - (S')*S;
    Em = eye(2) - S*(S');
    Ab = (Ab + Ab')/2;
    Em = (Em + Em')/2;

    % ---- Law1: eigenvalues vs 1-sigma^2 (best pairing) ----
    alpha_M = max(min(1 - sig.^2, 1), 0);
    eAb = sort(real(eig(Ab)), 'descend');
    eEm = sort(real(eig(Em)), 'descend');
    eAb = max(min(eAb,1),0);
    eEm = max(min(eEm,1),0);

    L1(g) = max([bestperm2_maxabs(alpha_M, eAb), bestperm2_maxabs(alpha_M, eEm)]);

    % ---- Use SVD of S for Laws 2/4 sampling ----
    [Us, Ss, Vs] = svd(S,'econ'); %#ok<ASGLU>
    sig2 = diag(Ss);
    alpha_M2 = max(min(1 - sig2.^2, 1), 0); %#ok<NASGU>

    % ---- Law2: equal-power-profile (matched eigenchannel weights) ----
    alpha_i = zeros(opt.nPairs,1);
    eps_o   = zeros(opt.nPairs,1);
    for k = 1:opt.nPairs
        c = randn(2,1)+1i*randn(2,1); c = c/norm(c);          % weights
        mags = abs(c); phs = 2*pi*rand(2,1);
        d = mags .* exp(1i*phs); d = d/norm(d);              % same mags, random phases
        i_vec = Vs*c;
        o_vec = Us*d;
        alpha_i(k) = real(i_vec' * Ab * i_vec);
        eps_o(k)   = real(o_vec' * Em * o_vec);
    end
    alpha_i = max(min(alpha_i,1),0);
    eps_o   = max(min(eps_o,1),0);

    denom2 = max(mean(alpha_i), 1e-12);
    L2(g) = sqrt(mean((alpha_i - eps_o).^2)) / denom2;

    alpha_i_all{g} = alpha_i;
    eps_o_all{g}   = eps_o;

    % ---- Law3: trace invariance under basis rotation ----
    TrAb = real(trace(Ab));
    TrEm = real(trace(Em));
    errTr = abs(TrAb - TrEm);
    % also test random basis invariance
    sumA = zeros(opt.nBases,1); sumE = zeros(opt.nBases,1);
    for k = 1:opt.nBases
        Qin  = rand_unitary2();
        Qout = rand_unitary2();
        sumA(k) = real(trace(Qin'  * Ab * Qin));
        sumE(k) = real(trace(Qout' * Em * Qout));
    end
    denom3 = max(abs(TrAb), 1e-12);
    L3(g) = (mean(abs(sumA - TrAb)) + mean(abs(sumE - TrEm)) + errTr) / denom3;

    % ---- Law4: reciprocity pairing alpha(i)=epsilon(i*) (NOT guaranteed in general 2-port) ----
    alpha_t = zeros(opt.nTest,1);
    eps_t   = zeros(opt.nTest,1);
    for k = 1:opt.nTest
        i = randn(2,1)+1i*randn(2,1); i = i/norm(i);
        istar = conj(i);
        alpha_t(k) = real(i'     * Ab * i);
        eps_t(k)   = real(istar' * Em * istar);
    end
    alpha_t = max(min(alpha_t,1),0);
    eps_t   = max(min(eps_t,1),0);

    denom4 = max(mean(alpha_t), 1e-12);
    L4(g) = sqrt(mean((alpha_t - eps_t).^2)) / denom4;

    alpha_t_all{g} = alpha_t;
    eps_t_all{g}   = eps_t;

    fprintf('pH=%.3g: L1=%.2e, L2=%.2e, L3=%.2e, L4=%.2e, passViol=%.2e\n', ...
        RES(g).pH, L1(g), L2(g), L3(g), L4(g), passViol(g));
end

% ===================== plotting =====================
fig = figure('Color','w','Units','centimeters','Position',[2 2 22.5 18]);
tlo = tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');

% helper for log-safe
logsafe = @(x) log10(max(x, 1e-16));

% (a) Law1
ax1 = nexttile(tlo,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, pHs, logsafe(L1), 'o-', 'LineWidth',2);
xlabel(ax1,'pH'); ylabel(ax1,'$\log_{10}(\mathrm{mismatch})$','Interpreter','latex');
title(ax1,'Law 1: $\alpha_{Mp}=\epsilon_{Mp}=1-\sigma_p^2$','Interpreter','latex');

% (b) Law2
ax2 = nexttile(tlo,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, pHs, logsafe(L2), 'o-', 'LineWidth',2);
xlabel(ax2,'pH'); ylabel(ax2,'$\log_{10}(\mathrm{RMS}/\mathrm{mean})$','Interpreter','latex');
title(ax2,'Law 2: equal-power-profile','Interpreter','latex');

% (c) Law3
ax3 = nexttile(tlo,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
plot(ax3, pHs, logsafe(L3), 'o-', 'LineWidth',2);
xlabel(ax3,'pH'); ylabel(ax3,'$\log_{10}(\mathrm{trace\ mismatch})$','Interpreter','latex');
title(ax3,'Law 3: trace invariance','Interpreter','latex');

% (d) Law4
ax4 = nexttile(tlo,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
plot(ax4, pHs, logsafe(L4), 'o-', 'LineWidth',2);
xlabel(ax4,'pH'); ylabel(ax4,'$\log_{10}(\mathrm{RMS}/\mathrm{mean})$','Interpreter','latex');
title(ax4,'Law 4: reciprocity pairing','Interpreter','latex');

% (e) passivity violation only (non-empty even when zero)
ax5 = nexttile(tlo,5); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');
plot(ax5, pHs, passViol, 's-', 'LineWidth',2);
yline(ax5,0,'k-','LineWidth',1,'HandleVisibility','off');
xlabel(ax5,'pH'); ylabel(ax5,'$\max(0,\max\sigma(S)-1)$','Interpreter','latex');
title(ax5,'passivity violation (at $\eta_{\rm eff}^*$)','Interpreter','latex');

% (f) example scatter: pick worst (max of L2 or L4)
[~, iw2] = max(L2);
[~, iw4] = max(L4);
useLaw2 = L2(iw2) >= L4(iw4);
if useLaw2
    gex = iw2;
    x = alpha_i_all{gex}; y = eps_o_all{gex};
    ex_tag(gex) = "worst Law2";
    ttl = sprintf('Example (%s): pH=%.3g', "worst Law2", pHs(gex));
    xlab = '$\alpha_i$'; ylab = '$\epsilon_o$';
else
    gex = iw4;
    x = alpha_t_all{gex}; y = eps_t_all{gex};
    ex_tag(gex) = "worst Law4";
    ttl = sprintf('Example (%s): pH=%.3g', "worst Law4", pHs(gex));
    xlab = '$\alpha(i)$'; ylab = '$\epsilon(i^*)$';
end

ax6 = nexttile(tlo,6); hold(ax6,'on'); box(ax6,'on'); grid(ax6,'on');
if ~isempty(x) && ~isempty(y)
    scatter(ax6, x, y, 14, 'filled', 'MarkerFaceAlpha',0.55);
    mm = max([x(:); y(:); 1e-12]);
    plot(ax6, [0 mm],[0 mm],'k--','LineWidth',1.4);
    axis(ax6,[0 mm 0 mm]);
else
    text(ax6,0.5,0.5,'No valid samples','Units','normalized','HorizontalAlignment','center');
end
xlabel(ax6, xlab, 'Interpreter','latex');
ylabel(ax6, ylab, 'Interpreter','latex');
title(ax6, ttl, 'Interpreter','none');

% panel labels
put_panel_label(ax1,'a'); put_panel_label(ax2,'b'); put_panel_label(ax3,'c');
put_panel_label(ax4,'d'); put_panel_label(ax5,'e'); put_panel_label(ax6,'f');

if ~isempty(opt.savePng)
    set(fig,'PaperPositionMode','auto');
    exportgraphics(fig, opt.savePng, 'Resolution', 600);
    fprintf('Saved: %s\n', opt.savePng);
end

% output table
T = table(pHs, L1, L2, L3, L4, passViol, 'VariableNames', ...
    {'pH','Law1','Law2','Law3','Law4','PassivityViolation'});
LAW = struct();
LAW.table = T;
LAW.metrics = struct('L1',L1,'L2',L2,'L3',L3,'L4',L4,'passViol',passViol,'etaStar',ex_etaStar);
disp(T);

end

% ===================== helper: build Gamma from (Jox,Jred) =====================
function [Gamma, rho, theta] = gamma_from_flux(Jox, Jred)
eps0 = 1e-30;
Jox = max(Jox,0); Jred = max(Jred,0);
chi = (Jred - Jox) ./ (Jred + Jox + eps0);
chi = max(min(chi,1),-1);
theta = acos(chi) - pi/2;                               % [-pi/2, +pi/2]
rho = sqrt( min(Jred,Jox) ./ (max(Jred,Jox) + eps0) );   % [0,1]
rho(~isfinite(rho)) = 0;
rho = max(min(rho,1),0);
Gamma = rho .* exp(1i*theta);
end

% ===================== helper: SU(2)-like 2x2 unitary =====================
function U = unitary2x2(a, b)
% U = [cos a, e^{i b} sin a; -e^{-i b} sin a, cos a]
ca = cos(a); sa = sin(a);
U = [ca, exp(1i*b)*sa; -exp(-1i*b)*sa, ca];
end

% ===================== helper: best pairing mismatch for 2 entries =====================
function m = bestperm2_maxabs(x, y)
x = x(:); y = y(:);
if numel(x)~=2 || numel(y)~=2, m = NaN; return; end
m1 = max(abs(x(1)-y(1)), abs(x(2)-y(2)));
m2 = max(abs(x(1)-y(2)), abs(x(2)-y(1)));
m = min(m1,m2);
end

% ===================== helper: random 2x2 unitary =====================
function Q = rand_unitary2()
Z = randn(2) + 1i*randn(2);
[Q,R] = qr(Z);
d = diag(R);
ph = d ./ max(abs(d), 1e-12);
Q = Q * diag(conj(ph));
end

% ===================== helper: passivity-safe clamp magnitude =====================
function z = clamp_mag(z, rmax)
r = abs(z);
bad = r > rmax;
if any(bad)
    z(bad) = z(bad) .* (rmax ./ (r(bad) + 1e-30));
end
end

function y = clamp01(y)
y = max(min(y,1),0);
end

function put_panel_label(ax, ch)
text(ax, -0.15, 0.98, ch, 'Units','normalized', ...
    'FontSize',16, 'FontWeight','bold', 'FontName','Helvetica', ...
    'Interpreter','none', 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');
end

function z = wrapToPi_local(x)
z = mod(x + pi, 2*pi) - pi;
end

function z = clip_local(x, maxExp)
z = min(max(x, -maxExp), maxExp);
end

function T = verify_law4_correlations_echem(RES, opt)
% verify_law4_correlations_echem
% ------------------------------------------------------------
% Correlate Law-4 reciprocity-pairing mismatch with operating-point diagnostics:
%   beta_w(eta*),  w(eta*)=atanh(beta_w),  rho(eta*),  I_store(eta*)
% where eta* is the cathodic optimum defined by maximizing Pi_dens.
%
% Requires in your path/workspace:
%   safe_two_mode_model
%   calc_caismith_from_params_improved
%   prctile_fallback
%   clip
%
% Output:
%   T: table with pH, etaStar, Law4 mismatch, beta*, w*, rho*, Istore*, Pi_dens_max, passivity margin

    if nargin < 2, opt = struct(); end
    opt = setdef(opt,'maxExp',70);
    opt = setdef(opt,'newtonMaxIter',10);
    opt = setdef(opt,'newtonTol',5e-11);

    opt = setdef(opt,'eta0',0.10);                 % Pi_dens regularization (must >0)
    opt = setdef(opt,'fixedEtaWindow',[-1.5 0.5]); % same as your plots
    opt = setdef(opt,'nEta',700);
    opt = setdef(opt,'cathodicOnly',true);

    % ----- Law4 Monte Carlo -----
    opt = setdef(opt,'nTest',800);     % random vectors for Law4 mismatch
    opt = setdef(opt,'rngSeed',11);

    % ----- "general 2-port" mixing controls (U != V) -----
    % These are the knobs that make Law4 nontrivial.
    opt = setdef(opt,'rot0', 0.20*pi);                 % base mixing angle
    opt = setdef(opt,'rotGain', 0.25*pi);              % mixing depends on beta*
    opt = setdef(opt,'nonrecipRotGain', 0.35*pi);      % makes U and V differ
    opt = setdef(opt,'phase0', 0.10*pi);               % diagonal phase in U/V
    opt = setdef(opt,'nonrecipPhaseGain', 0.25*pi);    % phase asymmetry between U and V

    nG = numel(RES);

    % storage
    pH = nan(nG,1);
    etaStar = nan(nG,1);
    Pi_dens_max = nan(nG,1);

    betaStar = nan(nG,1);
    wStar    = nan(nG,1);
    rhoStar  = nan(nG,1);
    IStar    = nan(nG,1);

    law4_mis = nan(nG,1);
    passMargin = nan(nG,1);  % max(0, sigma_max(S)-1) at eta*

    rng(opt.rngSeed);

    for g = 1:nG
        P = RES(g).P;
        pH(g) = RES(g).pH;

        % --- eta grid ---
        etaGrid = linspace(opt.fixedEtaWindow(1), opt.fixedEtaWindow(2), opt.nEta).';

        % --- waveguide diagnostics on this grid ---
        % outputs: [rho, theta, I_store, eta_eff, U, S_rel, beta_w, ...]
        [rho, theta, I_store, eta_eff, Uex, Srel, beta_w] = ...
            calc_caismith_from_params_improved(P, etaGrid, true, opt.maxExp, opt.newtonMaxIter, opt.newtonTol, 'none'); %#ok<ASGLU>

        rho     = max(min(rho(:),1),0);
        theta   = theta(:);
        I_store = max(min(I_store(:),1),0);
        eta_eff = eta_eff(:);
        beta_w  = max(min(beta_w(:),1-1e-12),-1+1e-12);

        % --- compute Pi_dens on the SAME etaGrid ---
        jTot = safe_two_mode_model(P, etaGrid, opt.maxExp, opt.newtonMaxIter, opt.newtonTol, []);
        eta_eff2 = etaGrid - P.Rohm.*(jTot(:)./1000);

        % Align rho to eta_eff2 (because eta_eff from mapping may be slightly different ordering)
        [eta_eff_s, ord] = sort(eta_eff);
        rho_s    = rho(ord);
        I_s      = I_store(ord);
        beta_s   = beta_w(ord);

        [eta_eff2s, ord2] = sort(eta_eff2);
        jTot_s = jTot(ord2);

        rho_i  = interp1(eta_eff_s, rho_s,  eta_eff2s, 'linear', 'extrap');
        I_i    = interp1(eta_eff_s, I_s,    eta_eff2s, 'linear', 'extrap');
        beta_i = interp1(eta_eff_s, beta_s, eta_eff2s, 'linear', 'extrap');

        rho_i  = max(min(rho_i,1),0);
        I_i    = max(min(I_i,1),0);
        beta_i = max(min(beta_i,1-1e-12),-1+1e-12);

        Pi_use  = max(0, -jTot_s) .* max(1 - rho_i.^2, 0);
        Pi_dens = Pi_use ./ (abs(eta_eff2s) + opt.eta0);

        if opt.cathodicOnly
            mask = (eta_eff2s < 0) & isfinite(Pi_dens) & (Pi_dens > 0);
        else
            mask = isfinite(Pi_dens) & (Pi_dens > 0);
        end

        if nnz(mask) < 10
            continue;
        end

        tmp = Pi_dens; tmp(~mask) = -Inf;
        [Pi_dens_max(g), kStar] = max(tmp);
        etaStar(g) = eta_eff2s(kStar);

        % values at eta* (from interpolated arrays on eta_eff2s)
        betaStar(g) = beta_i(kStar);
        wStar(g)    = atanh(betaStar(g));
        rhoStar(g)  = rho_i(kStar);
        IStar(g)    = I_i(kStar);

        % --- Build a general passive 2-port S at eta* ---
        % Use two modal "reflection" amplitudes GammaA/GammaB.
        % Here we take a simple physically consistent choice:
        %   GammaA = rhoA * exp(i*thetaA), GammaB = rhoB * exp(i*thetaB)
        % but we only have total rho/theta from combined fluxes.
        % So we approximate by splitting amplitude using "mode fractions" from UNSAT fluxes:
        %   fracA = U_A/(U_A+U_B), fracB = 1-fracA
        % and set rhoA=rho*sqrt(fracA), rhoB=rho*sqrt(fracB) (keeps passivity).
        %
        % If you already have per-mode GammaA/B in your code, replace this block with that.
        %
        % ---- compute UNSAT directional fluxes at etaStar (use eta_eff=etaStar) ----
        etaS = etaStar(g);
        % Make sure P has aA,bA,aB,bB (your fit struct does)
        JoxA  = P.jstarA .* exp( clip(P.aA.*etaS, opt.maxExp) );
        JredA = P.jstarA .* exp( clip(P.bA.*etaS, opt.maxExp) );
        JoxB  = P.jstarB .* exp( clip(P.aB.*etaS, opt.maxExp) );
        JredB = P.jstarB .* exp( clip(P.bB.*etaS, opt.maxExp) );

        UA = (JoxA + JredA);
        UB = (JoxB + JredB);
        fracA = UA / max(UA + UB, 1e-30);
        fracA = max(min(fracA,1),0);
        fracB = 1 - fracA;

        thetaS = interp1(eta_eff_s, theta(ord), etaS, 'linear', 'extrap'); % phase at eta*
        if ~isfinite(thetaS), thetaS = 0; end

        rhoS = rhoStar(g);
        rhoA = rhoS * sqrt(fracA);
        rhoB = rhoS * sqrt(fracB);

        GammaA = rhoA * exp(1i*thetaS);
        GammaB = rhoB * exp(1i*thetaS);

        D = diag([GammaA, GammaB]);

        % Unitary mixing U and V driven by betaStar (nonreciprocal if U != V)
        b = betaStar(g);
        aU = opt.rot0 + opt.rotGain*b;
        aV = aU + opt.nonrecipRotGain*b;

        phU = opt.phase0;
        phV = opt.phase0 + opt.nonrecipPhaseGain*b;

        Uu = unitary2x2(aU, phU);
        Vv = unitary2x2(aV, phV);

        S = Uu * D * (Vv');   % passive by construction if |GammaA|,|GammaB|<=1

        % passivity margin: max(svd(S))-1
        sig = svd(S);
        passMargin(g) = max(0, max(real(sig)) - 1);

        % --- Law4 mismatch at eta* ---
        law4_mis(g) = law4_mismatch( S, opt.nTest );
    end

    % ---- assemble table ----
    T = table(pH, etaStar, Pi_dens_max, betaStar, wStar, rhoStar, IStar, law4_mis, passMargin, ...
        'VariableNames', {'pH','etaEffStar','PiDensMax','betaStar','wStar','rhoStar','IstoreStar','Law4Mismatch','PassivityMargin'});

    disp(T);

    % ---- correlations & plots ----
    make_corr_figure(T);

end

% ===================== helpers =====================

function m = law4_mismatch(S, nTest)
% Law4 mismatch metric (RMS/mean) using Ab and Em:
%   Ab = I - S^H S,   Em = I - S S^H
%   alpha(i) = i^H Ab i,   epsilon(i*) = (i*)^H Em (i*)
% For reciprocal pairing, alpha(i) ≈ epsilon(i*). Deviation quantifies breaking.

    Ab = eye(2) - (S')*S;
    Em = eye(2) - S*(S');

    Ab = (Ab + Ab')/2;
    Em = (Em + Em')/2;

    a = zeros(nTest,1);
    e = zeros(nTest,1);

    for k = 1:nTest
        i = randn(2,1) + 1i*randn(2,1);
        i = i / norm(i);

        istar = conj(i);

        a(k) = real(i' * Ab * i);
        e(k) = real(istar' * Em * istar);
    end

    % clip to [0,1]
    a = min(max(a,0),1);
    e = min(max(e,0),1);

    denom = mean(0.5*(a+e)) + 1e-12;
    m = sqrt(mean((a - e).^2)) / denom;
end

function make_corr_figure(T)
    % filter finite rows
    good = isfinite(T.Law4Mismatch) & isfinite(T.betaStar) & isfinite(T.wStar) & ...
           isfinite(T.rhoStar) & isfinite(T.IstoreStar);
    TT = T(good,:);

    if height(TT) < 3
        warning('Not enough valid pH points for correlation plots.');
        return;
    end

    y = TT.Law4Mismatch;

    X = {TT.betaStar, TT.wStar, TT.rhoStar, TT.IstoreStar};
    xlab = {'$\beta_w(\eta^*)$', '$w(\eta^*)=\mathrm{atanh}(\beta_w)$', '$\rho(\eta^*)$', '$I_{\rm store}(\eta^*)$'};
    ttl  = {'Law4 vs $\beta_w^*$', 'Law4 vs $w^*$', 'Law4 vs $\rho^*$', 'Law4 vs $I_{\rm store}^*$'};

    fig = figure('Color','w','Units','centimeters','Position',[2 2 22.5 18]);
    tlo = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

    for k = 1:4
        ax = nexttile(tlo,k);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        x = X{k};

        % scatter
        scatter(ax, x, y, 60, 'filled');

        % linear fit
        p = polyfit(x, y, 1);
        xx = linspace(min(x), max(x), 50);
        yy = polyval(p, xx);
        plot(ax, xx, yy, 'k--', 'LineWidth', 1.8);

        % correlations
        R = corrcoef(x, y);
        rP = R(1,2);
        rS = corr(x, y, 'Type', 'Spearman');

        xlabel(ax, xlab{k}, 'Interpreter','latex');
        ylabel(ax, 'Law4 mismatch (RMS/mean)', 'Interpreter','none');
        title(ax, ttl{k}, 'Interpreter','latex');

        txt = sprintf('Pearson r=%.2f\nSpearman r=%.2f', rP, rS);
        text(ax, 0.02, 0.98, txt, 'Units','normalized', ...
            'HorizontalAlignment','left','VerticalAlignment','top', ...
            'FontSize',10,'BackgroundColor',[1 1 1 0.85], 'EdgeColor',[0.85 0.85 0.85]);

    end

    sgtitle(tlo, 'Correlating Law4 reciprocity mismatch with operating-point diagnostics at $\eta^*_{\rm eff}$', ...
        'Interpreter','latex', 'FontSize', 13);

end