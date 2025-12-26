% =========================================================================
% SCRIPT: Figure7.m
%
% PURPOSE: 
%   Generates Figure 7: Primary Synergies (A & B).
%   - Spatial Weights, Temporal Activation, Cosine Similarity.
%   - UPDATES: Robust path detection added.
%
% LAYOUT: 5 Rows x 4 Columns
%
% AUTHOR: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
% --- DYNAMIC PATH SETUP ---
scriptPath = fileparts(mfilename('fullpath'));
if isempty(scriptPath), scriptPath = pwd; end % Fallback for running sections
baseDir = fileparts(scriptPath); 

fprintf('Detected Base Directory: %s\n', baseDir);

matDir = fullfile(baseDir, 'Data', 'synergy'); 
outFigDir = fullfile(baseDir, 'outputFigures_Fig7');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

if ~exist(matDir, 'dir')
    error('Data folder not found at: %s\n(Did you download the ''Data'' folder from GitHub?)', matDir); 
end

monkeyNameList = {'Yachimun', 'Seseki'};
synergyList    = {'A', 'B'}; 
muscleLabels   = {'BRD', 'ECR', 'ECU', 'ED23', 'EDC', 'FCR', 'FCU', 'FDP', 'FDS', 'PL'};

% Landmark Days
LMD_A = [29, 64, 69, 79, 99];
LMD_B = [22, 36, 44, 48, 64];
colorsLMD = lines(5); 

colPreCos  = [0.00, 0.45, 0.74]; 
colPostCos = [0.85, 0.33, 0.10]; 

figWidth_cm  = 24; 
figHeight_cm = 30; 
fig = figure('Name', 'Figure 7: Primary Synergies', 'Units', 'centimeters', ...
             'Position', [2, 2, figWidth_cm, figHeight_cm], 'Color', 'w');
t = tiledlayout(5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 2. MAIN LOOPS
% -------------------------------------------------------------------------
for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);
    
    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        colOffset = 0; lbl = 'Monkey A'; currLMD = LMD_A; gapEnd = 29;
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        colOffset = 2; lbl = 'Monkey B'; currLMD = LMD_B; gapEnd = 22;
    end
    
    % Load Data
    dataPreSel  = loadMatData('Synergy', monkeyName, 'PreSelect',  matDir, synergyList, ttDate);
    dataPostSel = loadMatData('Synergy', monkeyName, 'PostSelect', matDir, synergyList, ttDate);
    dataPreAll  = loadSynergyWeights(monkeyName, 'PreAll', matDir, ttDate, synergyList);
    dataPostAll = loadSynergyWeights(monkeyName, 'PostAll', matDir, ttDate, synergyList);
    
    W_ref = cell(1,2);
    for i = 1:2
        % Calculate Ref from PreAll (Robust Mean)
        if ~isempty(dataPreAll(i).wAll)
            if iscell(dataPreAll(i).wAll)
                tempW = [];
                for d=1:length(dataPreAll(i).wAll)
                    if ~isempty(dataPreAll(i).wAll{d}), tempW = [tempW, mean(dataPreAll(i).wAll{d},2)]; end
                end
                W_ref{i} = mean(tempW, 2);
            else
                W_ref{i} = mean(dataPreAll(i).wAll, 2);
            end
        end
    end

    for synIdx = 1:2
        synLabel = synergyList{synIdx};
        
        % Rows 1-2: Spatial Weights
        row = synIdx;
        ax = nexttile((row-1)*4 + colOffset + 1);
        plotWeights(ax, dataPreSel(synIdx).dataWM, muscleLabels, 'Pre', synLabel, mIndx==1 && synIdx==1, synIdx==2);
        ax = nexttile((row-1)*4 + colOffset + 2);
        plotWeights(ax, dataPostSel(synIdx).dataWM, muscleLabels, 'Post', synLabel, false, synIdx==2);

        % Rows 3-4: Temporal Activation
        row = synIdx + 2;
        ax = nexttile((row-1)*4 + colOffset + 1);
        plotActivation(ax, dataPreSel(synIdx), 'Pre', synLabel, mIndx==1 && synIdx==1, false, [], []);
        ax = nexttile((row-1)*4 + colOffset + 2);
        plotActivation(ax, dataPostSel(synIdx), 'Post', synLabel, false, synIdx==2, currLMD, colorsLMD);

        % Row 5: Cosine Similarity
        row = 5;
        ax = nexttile((row-1)*4 + colOffset + synIdx);
        plotCosine(ax, dataPreAll(synIdx), dataPostAll(synIdx), W_ref{synIdx}, ...
                   synLabel, mIndx==1 && synIdx==1, gapEnd, colPreCos, colPostCos);
    end
end

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [figWidth_cm figHeight_cm], 'PaperPositionMode', 'auto');
try
    print(fig, fullfile(outFigDir, 'Figure7_PrimarySynergies.svg'), '-dsvg', '-painters');
    fprintf('Figure 7 saved to: %s\n', outFigDir);
catch
    exportgraphics(fig, fullfile(outFigDir, 'Figure7_PrimarySynergies.png'), 'Resolution', 300);
    fprintf('SVG failed, saved as PNG to: %s\n', outFigDir);
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotWeights(ax, wData, labels, prePost, synLbl, showYLab, showXLab)
    bar(ax, wData, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
    if showYLab, ylabel(ax, 'weight contribution [a.u.]', 'FontSize', 8); end
    if showXLab
        set(ax, 'XTick', 1:length(labels), 'XTickLabel', labels, 'XTickLabelRotation', 90, 'FontSize', 7);
    else
        set(ax, 'XTick', []);
    end
    title(ax, sprintf('Synergy %s (%s)', synLbl, prePost), 'FontSize', 9, 'FontWeight', 'bold');
    ylim(ax, [0 1.8]); box(ax, 'off'); set(ax, 'TickDir', 'out');
end

function plotActivation(ax, dataStruct, prePost, synLbl, showYLab, showLegend, LMD, colors)
    hold(ax, 'on');
    t = dataStruct.t;
    
    % REMOVED: Gray shadings for task phases
    xline(ax, 0, 'k:');
    
    if strcmp(prePost, 'Pre')
        plot(ax, t, dataStruct.dataM, 'k-', 'LineWidth', 1.5);
    else
        L = []; legTxt = {};
        for k = 1:size(dataStruct.dataM, 2)
            L(k) = plot(ax, t, dataStruct.dataM(:,k), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
            legTxt{k} = sprintf('Day %d', LMD(k));
        end
        if showLegend
            legend(L, legTxt, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');
        end
    end
    xlim(ax, [-15 15]);
    if showYLab, ylabel(ax, 'Amplitude [a.u.]', 'FontSize', 8); end
    if strcmp(prePost, 'Pre'), title(ax, sprintf('Synergy %s', synLbl), 'FontSize', 9, 'FontWeight', 'bold'); end
    xlabel(ax, 'Task range [%]', 'FontSize', 8);
    box(ax, 'off'); set(ax, 'TickDir', 'out');
end

function plotCosine(ax, dataPre, dataPost, wRef, synLbl, showYLab, gapEnd, colPre, colPost)
    hold(ax, 'on');
    fill(ax, [-100 200 200 -100], [0.95 0.95 1.05 1.05], [0.95 0.95 0.95], 'EdgeColor', 'none');
    fill(ax, [0 gapEnd gapEnd 0], [0 0 1.2 1.2], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    if ~isempty(dataPre.wAll) && ~isempty(wRef)
        sims = calculateCosineSeries(dataPre.wAll, wRef);
        plot(ax, dataPre.daysDiff, sims, '-o', 'Color', colPre, 'MarkerFaceColor', colPre, 'MarkerSize', 4, 'LineWidth', 1);
    end
    if ~isempty(dataPost.wAll) && ~isempty(wRef)
        sims = calculateCosineSeries(dataPost.wAll, wRef);
        plot(ax, dataPost.daysDiff, sims, '-o', 'Color', colPost, 'MarkerFaceColor', colPost, 'MarkerSize', 4, 'LineWidth', 1);
    end
    
    yline(ax, 1.0, 'k:', 'LineWidth', 0.8); xline(ax, 0, 'k--', 'LineWidth', 0.75);
    ylim(ax, [0.45 1.02]); xlim(ax, [-20 130]);
    title(ax, ['Synergy ' synLbl], 'FontSize', 9, 'FontWeight', 'bold');
    if showYLab
        ylabel(ax, 'Cosine Similarity', 'FontSize', 8);
        legend(ax, {'Pre-Stability', 'Post-Adaptation'}, 'Location', 'southeast', 'FontSize', 7, 'Box', 'off');
    else
        set(ax, 'YTickLabel', []);
    end
    xlabel(ax, 'Post surgery days', 'FontSize', 8);
    box(ax, 'off'); set(ax, 'TickDir', 'out');
end

function dataAll = loadMatData(~, monkeyName, condition, matDir, focusList, ttDate)
    matPath  = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    dataFile = load(matPath, 'synergyData', 'percentTime', 'expDates','nameList', 'synergyWData');
    
    nList = length(focusList);
    focusListNum = zeros(1,nList);
    for i=1:nList, focusListNum(i) = find(strcmp(focusList{i}, dataFile.nameList)); end
    
    dataAll  = struct([]);
    for i = 1:nList
        expDates = datetime(cellstr(dataFile.expDates), 'InputFormat','yyMMdd');
        meanData = zeros(length(dataFile.percentTime), length(expDates));
        
        synBlockFirst = dataFile.synergyWData{focusListNum(i)};
        if iscell(synBlockFirst), sampleW = synBlockFirst{1}; else, sampleW = synBlockFirst; end
        meanW = zeros(size(sampleW, 1), length(expDates));
        
        for dayI = 1:length(expDates)
            meanData(:,dayI) = mean(dataFile.synergyData{focusListNum(i),dayI},1)';
            
            synBlock = dataFile.synergyWData{focusListNum(i)};
            if iscell(synBlock)
                if dayI <= length(synBlock), rawW = synBlock{dayI}; else, rawW = []; end
            else
                if dayI <= size(synBlock, 2), rawW = synBlock(:, dayI); else, rawW = []; end
            end
            
            if ~isempty(rawW)
                if size(rawW, 2) > 1, meanW(:,dayI) = mean(rawW,2); else, meanW(:,dayI) = rawW; end
            end
        end
        dataAll(i).t         = dataFile.percentTime;
        dataAll(i).dataM     = meanData;
        dataAll(i).dataWM    = meanW;
        dataAll(i).daysDiff  = days(expDates - ttDate);
    end
end

function dataAll = loadSynergyWeights(monkeyName, condition, matDir, ttDate, focusList)
    fname = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    dataAll = repmat(struct('wAll',[], 'daysDiff',[]), 1, length(focusList));
    try
        D = load(fname, 'synergyWData', 'expDates', 'nameList');
        expDates = datetime(cellstr(D.expDates), 'InputFormat','yyMMdd');
        daysDiff = days(expDates - ttDate);
        for i = 1:length(focusList)
            idx = find(strcmp(focusList{i}, D.nameList));
            dataAll(i).daysDiff = daysDiff;
            dataAll(i).wAll = D.synergyWData{idx}; 
        end
    catch
        warning('Error loading %s', fname);
    end
end

function sims = calculateCosineSeries(wData, wRef)
    if iscell(wData)
        nDays = length(wData);
    else
        nDays = size(wData, 2);
    end
    
    sims = nan(1, nDays);
    
    for d = 1:nDays
        if iscell(wData)
            vec = wData{d};
        else
            vec = wData(:, d);
        end
        
        if ~isempty(vec)
            if size(vec, 2) > 1, vec = mean(vec, 2); end
            if norm(vec) > 0 && norm(wRef) > 0
                sims(d) = (vec(:)' * wRef(:)) / (norm(vec) * norm(wRef));
            end
        end
    end
end