% =========================================================================
% SCRIPT: Plot_FigureS9_SynergyLag.m
%
% PURPOSE: 
%   Generates Figure S9: Synergy Lag Analysis.
%   Calculates the time lag of peak correlation between Pre-surgery 
%   synergy templates and Post-surgery synergy activations.
%
% LAYOUT:
%   - 4x2 Grid (Monkey A Top, Monkey B Bottom).
%   - Left Axis: Lag in % of movement cycle (Solid lines).
%   - Right Axis: Behavioral Recovery (Dashed lines).
%
% AUTHOR: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
baseDir = 'C:\Users\mypre\Documents\Manuscripts\Revision\post acceptance revision\Philipp_eLife_2025';

% Input Directories
matDir  = fullfile(baseDir, 'Data', 'synergy'); 
behDirA = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
behDirB = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

% Output Directory
outFigDir =  fullfile(baseDir, 'outputFigures_FigS9');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Settings
calcType = 'Synergy'; 
monkeyNameList = {'Yachimun', 'Seseki'};
focusList = {'A', 'B'};

% Landmark Days
LMD_A = [29, 64, 69, 79, 99];
LMD_B = [22, 36, 44, 48, 64];

%% 2. DATA ANALYSIS LOOP
% -------------------------------------------------------------------------
monkey_data = cell(length(monkeyNameList), 1);

for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);
    
    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd'); 
        behPath = behDirA;
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd'); 
        behPath = behDirB; 
    end
   
    % Load Synergy Data
    [dataPreMean,  dataPreAll]  = loadMatData(calcType, monkeyName, 'PreAll',  matDir, focusList, ttDate);
    [dataPostMean, dataPostAll] = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);
    
    % Calculate Lags
    monkey_data{mIndx}.lags = calcLagAll(dataPreMean, dataPostMean, dataPostAll, dataPreAll, focusList);
    
    % Load Behavior
    monkey_data{mIndx}.beh  = loadBehavioralData(monkeyName, behPath);
end

%% 3. PLOTTING
% -------------------------------------------------------------------------
fprintf('Generating Figure S9...\n');

% Figure Setup (A4 Vertical approx. 13.5cm x 27cm for compact columns)
fig = figure('Name', 'Figure S9: Synergy Lag', 'Units', 'centimeters', ...
             'Position', [2 1 13.5 27], 'Color', 'w');

t = tiledlayout(4, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

% Plot Order: [2, 4, 1, 3] corresponds to specific Synergy pairs (B-B, B-A, A-B, A-A)
plotOrder = [2, 4, 1, 3];
titles = {'Post-Syn A vs Pre-Syn B', 'Post-Syn B vs Pre-Syn B', ...
          'Post-Syn A vs Pre-Syn A', 'Post-Syn B vs Pre-Syn A'};

for mIndx = 1:2
    D = monkey_data{mIndx};
    
    if mIndx == 1
        marks    = LMD_A;
        rowStart = 0; 
        xLimVal  = [-45 120];
        behLbl   = 'Off-Target (s)';
    else
        marks    = LMD_B;
        rowStart = 4; 
        xLimVal  = [-10 70];
        behLbl   = 'Contact Time (s)';
    end
    
    for i = 1:4
        idx = plotOrder(i);
        ax = nexttile(rowStart + i);
        
        plotDualAxis(ax, D.lags.daysDiffList, ...
                     D.lags.emgPeakLagMList{idx}, D.lags.emgPeakLagSDList{idx}, ...
                     D.beh.days, D.beh.mean, D.beh.std, ...
                     titles{i}, marks, xLimVal, behLbl);
                 
        % Formatting Adjustments
        if mIndx == 1 && i <= 2
            set(ax, 'XTickLabel', []); % Hide X labels for top monkey
        end 
        
        if mod(i,2) == 1 
            % Left Column: Show Y-Label Left
            yyaxis(ax, 'left'); ylabel(ax, 'Lag (%)', 'FontWeight', 'bold', 'Color', 'k');
            yyaxis(ax, 'right'); set(ax, 'YTickLabel', []); ylabel(ax, '');
        else 
            % Right Column: Show Y-Label Right
            yyaxis(ax, 'left'); set(ax, 'YTickLabel', []); ylabel(ax, '');
            yyaxis(ax, 'right'); ylabel(ax, behLbl, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4]);
        end
    end
end

% Save Results
exportgraphics(fig, fullfile(outFigDir, 'Figure_S9_LagAnalysis.png'), 'Resolution', 300);
try
    print(fig, fullfile(outFigDir, 'Figure_S9_LagAnalysis.svg'), '-dsvg', '-painters');
    fprintf('Figure S9 saved to: %s\n', outFigDir);
catch
    fprintf('SVG save failed. Figure available in PNG format.\n');
end


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotDualAxis(ax, xLag, yLag, sdLag, xBeh, yBeh, sdBeh, titleTxt, marks, xLims, behLbl)
    % --- RIGHT AXIS: BEHAVIOR (Dashed Gray) ---
    yyaxis(ax, 'right'); hold(ax, 'on');
    hasBeh = ~isempty(xBeh) && any(~isnan(yBeh));
    
    if hasBeh
        % Pre-Surgery Behavior
        idxPreB = xBeh < 0; 
        if any(idxPreB)
            plotShadedSingle(ax, xBeh(idxPreB), yBeh(idxPreB), sdBeh(idxPreB), [0.6 0.6 0.6], 0.2); 
        end
        
        % Post-Surgery Behavior (Smoothed)
        idxPostB = xBeh >= 0; 
        if any(idxPostB)
            yB_post = movmean(yBeh(idxPostB), 5); 
            sdB_post = movmean(sdBeh(idxPostB), 5);
            plotShadedSingle(ax, xBeh(idxPostB), yB_post, sdB_post, [0.6 0.6 0.6], 0.2);
        end
        
        yMax = max(yBeh + sdBeh, [], 'all', 'omitnan') * 1.3; 
        if yMax <= 0, yMax = 1; end
        ylim(ax, [0 yMax]);
    else
        ylim(ax, [0 1]);
    end
    set(ax, 'YColor', [0.4 0.4 0.4]);
    
    % --- LEFT AXIS: LAG (Solid Color) ---
    yyaxis(ax, 'left'); hold(ax, 'on'); 
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    idxPre = xLag <= 0; idxPost = xLag > 0;
    
    % Pre-Surgery (Blue Tube)
    if any(idxPre)
        plotTube(ax, xLag(idxPre), yLag(idxPre), sdLag(idxPre), [0.2 0.2 0.8], 0.3); 
    end
    
    % Post-Surgery (Red Tube)
    if any(idxPost)
        plotTube(ax, xLag(idxPost), yLag(idxPost), sdLag(idxPost), [0.8 0.2 0.2], 0.3); 
    end
    
    % Landmarks
    plotable = marks(marks >= xLims(1) & marks <= xLims(2));
    if ~isempty(plotable)
        yPos = repmat(-14, size(plotable)); % Place markers at bottom
        plot(ax, plotable, yPos, 'v', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'Clipping', 'off');
    end
    
    title(ax, titleTxt, 'FontWeight', 'normal', 'FontSize', 10);
    xlim(ax, xLims); ylim(ax, [-15 15]); 
    set(ax, 'YColor', 'k');
    
    xline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); 
    box(ax, 'off'); hold(ax, 'off');
end

function plotTube(ax, x, y_line, y_std, col, alphaVal)
    % Plots Solid Line with Error Tube
    x = x(:)'; y_line = y_line(:)'; y_std = y_std(:)';
    y_upper = y_line + y_std; y_lower = y_line - y_std;
    
    x_fill = [x, fliplr(x)]; y_fill = [y_upper, fliplr(y_lower)];
    nan_idx = isnan(y_fill); x_fill(nan_idx) = []; y_fill(nan_idx) = [];
    
    if ~isempty(x_fill)
        fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    end
    plot(ax, x, y_line, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '-'); 
end

function plotShadedSingle(ax, x, y, sd, col, alphaVal)
    % Plots Dashed Line with Error Shade
    x = x(:)'; y = y(:)'; sd = sd(:)'; 
    x_fill = [x, fliplr(x)]; y_fill = [y+sd, fliplr(y-sd)];
    nanIdx = isnan(y_fill); x_fill(nanIdx)=[]; y_fill(nanIdx)=[];
    
    if ~isempty(x_fill)
        fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    end
    plot(ax, x, y, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '--'); 
end

% --- DATA LOADING ---
function [dataMean, dataAll] = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    fName = sprintf('synergyData_%s_%s.mat', monkeyName, condition);
    matPath = fullfile(matDir, fName);
    
    if ~exist(matPath, 'file')
        warning('File not found: %s', matPath); 
        dataMean=[]; dataAll=[]; return; 
    end
    
    dataFile = load(matPath);
    dataCell = dataFile.synergyData;
    nList = length(focusList);
    
    % Find Synergy Indices
    focusListNum = zeros(1,nList);
    for i = 1:nList
        match = find(contains(dataFile.nameList, focusList{i}));
        if ~isempty(match), focusListNum(i) = match(1); end
    end
    
    expDates = datetime(cellstr(dataFile.expDates), 'InputFormat','yyMMdd');
    daysDiff = days(expDates - ttDate);
    t = dataFile.percentTime;
    
    for i = 1:nList
        meanData = zeros(length(t), length(daysDiff));
        allData  = cell(1, length(daysDiff));
        
        for dayI = 1:length(daysDiff)
            idx = focusListNum(i);
            if idx > 0 && ~isempty(dataCell{idx, dayI})
                meanData(:,dayI) = mean(dataCell{idx, dayI}, 1)';
                allData{dayI}    = (dataCell{idx, dayI})';
            else
                meanData(:,dayI) = nan;
            end
        end
        dataMean(i).data = meanData; dataAll(i).dataAll = allData;
        dataMean(i).daysDiff = daysDiff; dataAll(i).daysDiff = daysDiff;
    end
end

function behData = loadBehavioralData(monkeyName, behDir)
    behData.days = []; behData.mean = []; behData.std = [];
    
    if ~exist(behDir, 'dir'), return; end
    
    % Load Timeline
    f_days = dir(fullfile(behDir, 'days*.csv'));
    if isempty(f_days), return; end
    daysRaw = csvread(fullfile(behDir, f_days(1).name), 0, 0);
    
    % Load Files
    if strcmpi(monkeyName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(fullfile(behDir, pat));
    
    [~, idx] = sort({files.name}); files = files(idx);
    
    % Unique File Prefixes
    c = {}; 
    for k=1:length(files)
        name = files(k).name;
        if length(name) >= 14, c{1,k} = name(1:14); else, c{1,k} = name; end
    end
    uFiles = unique(c, 'stable');
    
    meanC = nan(1, length(uFiles)); 
    stdC  = nan(1, length(uFiles));
    
    for j = 1:length(uFiles)
        mIdx = find(startsWith({files.name}, uFiles{j}), 1);
        if isempty(mIdx), continue; end
        try
            dat = csvread(fullfile(behDir, files(mIdx).name), rOff, 1); 
            dat = dat(:);
            meanC(j) = nanmean(dat); 
            stdC(j) = nanstd(dat);
        catch, end
    end
    
    minLen = min(length(meanC), length(daysRaw));
    behData.days = daysRaw(1:minLen)';
    behData.mean = meanC(1:minLen);
    behData.std  = stdC(1:minLen);
end

% --- CALCULATIONS ---
function results = calcLagAll(emgPreMean, emgPostMean, emgPostAll, emgPreAll, focusList)
    nFocus = length(focusList); 
    nTotal = nFocus^2; 
    
    emgPeakLagMList = cell(1, nTotal); 
    emgPeakLagSDList = cell(1, nTotal);
    
    n = 1; 
    interpN = 200; 
    calcTimeRange = round((0.5-0.15)*interpN):round((0.5+0.15)*interpN);
    
    for i = 1:nFocus 
        for preIndx = 1:nFocus 
            % Reference Signal (Pre-Surgery Mean)
            emgPre = emgPreMean(preIndx).data; 
            emgPreInterpM = mean(interp1(1:size(emgPre,1), emgPre, linspace(1, size(emgPre,1), interpN), 'linear'), 2);
            refSig = emgPreInterpM(calcTimeRange) - mean(emgPreInterpM(calcTimeRange));
            
            % Helper to calculate lags for a dataset
            calcLags = @(targetStruct) compute_lags(targetStruct, refSig, interpN, calcTimeRange);
            
            [lagPostM, lagPostSD] = calcLags(emgPostAll(i));
            [lagPreM,  lagPreSD]  = calcLags(emgPreAll(i));
            
            emgPeakLagMList{n}  = [lagPreM, lagPostM]; 
            emgPeakLagSDList{n} = [lagPreSD, lagPostSD];
            n = n + 1;
        end
    end
    
    d1 = emgPreMean(1).daysDiff; d2 = emgPostMean(1).daysDiff;
    results.daysDiffList = [d1(:)', d2(:)'];
    results.emgPeakLagMList = emgPeakLagMList; 
    results.emgPeakLagSDList = emgPeakLagSDList;
end

function [lagM, lagSD] = compute_lags(targetStruct, refSig, interpN, calcTimeRange)
    dataAll = targetStruct.dataAll;
    
    if isempty(dataAll) || isempty(dataAll{1})
        lagM = []; lagSD = []; return;
    end
    
    nDays = length(dataAll);
    maxTrials = 0; 
    for d = 1:nDays, if ~isempty(dataAll{d}), maxTrials = max(maxTrials, size(dataAll{d}, 2)); end; end
    
    lag_trials = nan(nDays, maxTrials);
    
    for d = 1:nDays
        if isempty(dataAll{d}), continue; end
        trials = dataAll{d};
        nTrials = size(trials, 2);
        lenT = size(trials, 1);
        tInterp = linspace(1, lenT, interpN);
        
        for t = 1:nTrials
            trialSig = interp1(1:lenT, trials(:,t), tInterp, 'linear')';
            sigSeg = trialSig(calcTimeRange) - mean(trialSig(calcTimeRange));
            
            [c, lags] = xcorr(sigSeg, refSig, 0.15*interpN, 'coeff');
            [~, maxIdx] = max(c); 
            
            % Lag in % of cycle
            lag_trials(d, t) = lags(maxIdx) / interpN * 100;
        end
    end
    
    lagM  = nanmean(lag_trials, 2)';
    lagSD = nanstd(lag_trials, 0, 2)';
end