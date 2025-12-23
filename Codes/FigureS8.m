% =========================================================================
% SCRIPT: FigureS8.m
%
% PURPOSE: 
%   Generates Supplementary Figure S8: Pre vs Final-Day Synergy Profiles.
%
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------

% This ensures it works on any computer without hardcoding 'C:\Users...'
baseDir = 'C:\Users\mypre\Documents\Manuscripts\Revision\post acceptance revision\Philipp_eLife_2025';

% Input Directory (Data/synergy)
matDir = fullfile(baseDir, 'Data', 'synergy'); 

% Output Directory (outputFigures_FigS8)
outFigDir = fullfile(baseDir, 'outputFigures_FigS8');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify Data Path Exists
if ~exist(matDir, 'dir')
    error('Data folder not found at: %s\nPlease check directory structure.', matDir);
end

% Settings
monkeyNameList = {'Yachimun', 'Seseki'};
synergyLabels  = {'A', 'B', 'C', 'D'};
interpN        = 200;   
alphaVal       = 0.05;  
doBonferroni   = true;  

% Colors
colPre  = [0.00, 0.45, 0.74]; 
colPost = [0.85, 0.33, 0.10]; 

%% 2. MAIN PROCESSING LOOP
% -------------------------------------------------------------------------
all_data = cell(2, 4); 

for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);
    
    % Load Data from the specific synergy directory
    dataPre  = loadSynergyData(monkeyName, 'PreAll', matDir);
    dataPost = loadSynergyData(monkeyName, 'PostAll', matDir);
    
    for synIdx = 1:4
        % --- A. Gather Pre-Surgery Data ---
        preTrials = [];
        if ~isempty(dataPre) && synIdx <= length(dataPre) && iscell(dataPre(synIdx).dataAll)
            for d = 1:length(dataPre(synIdx).dataAll)
                dailyData = dataPre(synIdx).dataAll{d};
                if ~isempty(dailyData)
                    preTrials = [preTrials, normalizeTime(dailyData, interpN)]; %#ok<AGROW>
                end
            end
        end
        
        % --- B. Gather Post-Surgery Data (Last Day) ---
        postTrials = [];
        if ~isempty(dataPost) && synIdx <= length(dataPost) && iscell(dataPost(synIdx).dataAll)
            % Check valid days
            validDays = find(~cellfun(@isempty, dataPost(synIdx).dataAll));
            if ~isempty(validDays)
                lastDayIdx = validDays(end);
                dailyData = dataPost(synIdx).dataAll{lastDayIdx};
                postTrials = normalizeTime(dailyData, interpN);
            end
        end
        
        % --- C. Statistical Testing ---
        p_values = nan(1, interpN);
        sigMask = zeros(1, interpN);
        
        if ~isempty(preTrials) && ~isempty(postTrials)
            sigThresh = alphaVal;
            if doBonferroni, sigThresh = alphaVal / interpN; end
            for t = 1:interpN
                p_values(t) = ranksum(preTrials(t,:), postTrials(t,:));
            end
            sigMask = p_values < sigThresh;
        end
        
        % Store Results
        S.preMean  = mean(preTrials, 2, 'omitnan');
        S.preSD    = std(preTrials, 0, 2, 'omitnan');
        S.postMean = mean(postTrials, 2, 'omitnan');
        S.postSD   = std(postTrials, 0, 2, 'omitnan');
        S.sigMask  = sigMask;
        S.timeAxis = linspace(-15, 15, interpN);
        
        all_data{mIndx, synIdx} = S;
    end
end

%% 3. PLOTTING
% -------------------------------------------------------------------------
fprintf('Generating Figure S8...\n');

figWidth_cm  = 18;  
figHeight_cm = 12; 
fig = figure('Name', 'Figure S8: Synergy Activation', 'Units', 'centimeters', ...
             'Position', [5, 5, figWidth_cm, figHeight_cm], 'Color', 'w');

t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for mIndx = 1:2
    if mIndx == 1, lbl = 'Monkey A'; else, lbl = 'Monkey B'; end
    
    for synIdx = 1:4
        tileIdx = (mIndx-1)*4 + synIdx;
        ax = nexttile(tileIdx);
        
        S = all_data{mIndx, synIdx};
        if isempty(S) || isempty(S.preMean)
            axis(ax, 'off'); 
            continue; 
        end
        
        hold(ax, 'on');
        plotShaded(ax, S.timeAxis, S.preMean, S.preSD, colPre);
        plotShaded(ax, S.timeAxis, S.postMean, S.postSD, colPost);
        
        yMax = max([S.preMean+S.preSD; S.postMean+S.postSD], [], 'all', 'omitnan');
        if isempty(yMax), yMax = 1; end
        sigLevel = yMax * 1.1; 
        
        if any(S.sigMask)
            plot(ax, S.timeAxis(S.sigMask), repmat(sigLevel, 1, sum(S.sigMask)), ...
                 '.', 'Color', 'k', 'MarkerSize', 5);
        end
        
        title(ax, ['Synergy ' synergyLabels{synIdx}], 'FontWeight', 'bold', 'FontSize', 9, 'FontName', 'Arial');
        if synIdx == 1
            ylabel(ax, sprintf('%s\nActivation (a.u.)', lbl), 'FontWeight', 'bold', 'FontSize', 9, 'FontName', 'Arial');
        else
            set(ax, 'YTickLabel', []);
        end
        
        if mIndx == 2
            xlabel(ax, 'Task Range (%)', 'FontSize', 9, 'FontName', 'Arial');
        else
            set(ax, 'XTickLabel', []);
        end
        
        xlim(ax, [-15 15]);
        xline(ax, 0, 'k:', 'LineWidth', 0.8); 
        ylim(ax, [0, sigLevel * 1.15]);
        
        set(ax, 'TickDir', 'out', 'Box', 'off', 'FontSize', 7, 'LineWidth', 0.75, 'FontName', 'Arial');
        hold(ax, 'off');
    end
end

% --- Safe Legend ---
hDum = axes(fig, 'Visible', 'off'); 
hold(hDum, 'on');
l1 = plot(hDum, nan, nan, 'Color', colPre, 'LineWidth', 2);
l2 = plot(hDum, nan, nan, 'Color', colPost, 'LineWidth', 2);
l3 = plot(hDum, nan, nan, '.', 'Color', 'k', 'MarkerSize', 10);

hL = legend(hDum, [l1, l2, l3], {'Pre-Surgery', 'Final Post-Day', 'Significant'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', ...
       'FontSize', 8, 'Box', 'off');
hL.Position = [0.35, 0.02, 0.3, 0.05]; 

% --- Export ---
set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto', 'PaperSize', [figWidth_cm figHeight_cm]);
print(fig, fullfile(outFigDir, 'FigureS8_SynergyActivation.svg'), '-dsvg', '-painters');
exportgraphics(fig, fullfile(outFigDir, 'FigureS8_SynergyActivation.png'), 'Resolution', 300);

fprintf('Figure S8 saved to: %s\n', outFigDir);

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotShaded(ax, x, y, sd, col)
    x = x(:)'; y = y(:)'; sd = sd(:)';
    y_upper = y + sd; y_lower = y - sd; y_lower(y_lower < 0) = 0; 
    x_fill = [x, fliplr(x)]; y_fill = [y_upper, fliplr(y_lower)];
    nanIdx = isnan(y_fill); x_fill(nanIdx) = []; y_fill(nanIdx) = [];
    if ~isempty(x_fill)
        fill(ax, x_fill, y_fill, col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(ax, x, y, 'Color', col, 'LineWidth', 1.5);
end

function dataOut = normalizeTime(dataIn, nPoints)
    if size(dataIn, 2) > size(dataIn, 1) && size(dataIn, 2) > 50, dataIn = dataIn'; end
    [nT, nTrials] = size(dataIn);
    if nT == nPoints
        dataOut = dataIn;
    else
        dataOut = zeros(nPoints, nTrials);
        t_orig = linspace(0, 1, nT);
        t_new  = linspace(0, 1, nPoints);
        for tr = 1:nTrials
            dataOut(:,tr) = interp1(t_orig, dataIn(:,tr), t_new, 'linear');
        end
    end
end

function dataStruct = loadSynergyData(monkeyName, condition, targetDir)
    fNameShort = ['synergyData_' monkeyName '_' condition '.mat'];
    fullPath = fullfile(targetDir, fNameShort);
    
    if ~exist(fullPath, 'file')
        warning('Data file NOT found: %s', fullPath);
        dataStruct = [];
        return;
    end
    
    % Load Data
    try
        D = load(fullPath, 'synergyData');
        if isfield(D, 'synergyData')
            rawCell = D.synergyData;
            for i = 1:4
                dataStruct(i).dataAll = rawCell(i, :);
            end
        else
            warning('Variable "synergyData" missing in %s', fullPath);
            dataStruct = [];
        end
    catch ME
        warning('Error loading %s: %s', fullPath, ME.message);
        dataStruct = [];
    end
end