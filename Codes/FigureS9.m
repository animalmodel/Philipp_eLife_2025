% =========================================================================
% SCRIPT: FigureS9.m
%
% PURPOSE:
%   Generates Supplementary Figure S9: Synergy Lag Analysis.
%   - Layout: Fixed 4x2 Grid (Monkey A Top, Monkey B Bottom).
%   - Page Format: DIN A4 (Vertical), Square Panels.
%   - Styling:
%       1. Lag Data (Left Axis): Blue/Red lines are SOLID ('-').
%       2. Behavior (Right Axis): Gray Envelope is DASHED ('--').
%       3. Markers: Triangles reduced to Size 8.
%       4. Gap: No visual connection across Day 0.
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

% Input Directories
matDir  = fullfile(baseDir, 'Data', 'synergy'); 
behDirA = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
behDirB = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

% Output Directory
outFigDir = fullfile(baseDir, 'outputFigures_FigS9');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify Data Paths
if ~exist(matDir, 'dir')
    error('Data folder not found at: %s\n(Did you download the ''Data'' folder from GitHub?)', matDir); 
end

% Settings
calcType = 'Synergy'; 
monkeyNameList = {'Yachimun', 'Seseki'};
focusList = {'A', 'B'};

%% 2. ANALYSIS
% -------------------------------------------------------------------------
monkey_data = cell(length(monkeyNameList), 1);
for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);
    
    if mIndx==1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd'); 
        behPath = behDirA;
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd'); 
        behPath = behDirB; 
    end
   
    [dataPreMean,  dataPreAll]  = loadMatData(calcType, monkeyName, 'PreAll',  matDir, focusList, ttDate);
    [dataPostMean, dataPostAll] = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);
    
    [monkey_data{mIndx}.lags] = calcLagAll(dataPreMean, dataPostMean, dataPostAll, dataPreAll, focusList);
    [monkey_data{mIndx}.beh]  = loadBehavioralData(monkeyName, behPath);
end

%% 3. PLOTTING
% -------------------------------------------------------------------------
fprintf('Generating Figure S9...\n');

landmarkDaysA = [29, 64, 69, 79, 99];
landmarkDaysB = [22, 36, 44, 48, 64];

% A4 Dimensions (13.5x27cm) to ensure square tiles in 4x2
fig = figure('Name', 'Figure S9: Synergy Lag', 'Units', 'centimeters', 'Position', [2 1 13.5 27], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

plotOrder = [2, 4, 1, 3];
titles = {'Post-Syn A vs Pre-Syn B', 'Post-Syn B vs Pre-Syn B', 'Post-Syn A vs Pre-Syn A', 'Post-Syn B vs Pre-Syn A'};

for mIndx = 1:2
    D = monkey_data{mIndx};
    if mIndx == 1
        marks = landmarkDaysA;
        rowStart = 0; 
        xLimVal = [-45 120];
        behLbl = 'Off-Target (s)';
    else
        marks = landmarkDaysB;
        rowStart = 4; 
        xLimVal = [-10 70];
        behLbl = 'Contact Time (s)';
    end
    
    for i = 1:4
        idx = plotOrder(i);
        ax = nexttile(rowStart + i);
        
        plotDualAxis(ax, D.lags.daysDiffList, ...
                     D.lags.emgPeakLagMList{idx}, D.lags.emgPeakLagSDList{idx}, ...
                     D.beh.days, D.beh.mean, D.beh.std, ...
                     titles{i}, marks, xLimVal, behLbl);
                 
        if mIndx == 1 && i <= 2, set(ax, 'XTickLabel', []); end 
        if mod(i,2) == 1 
            yyaxis(ax, 'left'); ylabel(ax, 'Lag (%)', 'FontWeight', 'bold', 'Color', 'k');
            yyaxis(ax, 'right'); set(ax, 'YTickLabel', []); ylabel(ax, '');
        else 
            yyaxis(ax, 'left'); set(ax, 'YTickLabel', []); ylabel(ax, '');
            yyaxis(ax, 'right'); ylabel(ax, behLbl, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4]);
        end
    end
end

% Save
try
    print(fig, fullfile(outFigDir, 'FigureS9_LagAnalysis.svg'), '-dsvg', '-painters');
    fprintf('Figure S9 saved to: %s\n', outFigDir);
catch
    exportgraphics(fig, fullfile(outFigDir, 'FigureS9_LagAnalysis.png'), 'Resolution', 300);
    fprintf('SVG save failed. Figure available in PNG at: %s\n', outFigDir);
end


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotDualAxis(ax, xLag, yLag, sdLag, xBeh, yBeh, sdBeh, titleTxt, marks, xLims, behLbl)
    % RIGHT AXIS: BEHAVIOR (DASHED)
    yyaxis(ax, 'right'); hold(ax, 'on');
    hasBeh = ~isempty(xBeh) && any(~isnan(yBeh));
    if hasBeh
        idxPreB = xBeh < 0; if any(idxPreB), plotShadedSingle(ax, xBeh(idxPreB), yBeh(idxPreB), sdBeh(idxPreB), [0.6 0.6 0.6], 0.2); end
        idxPostB = xBeh >= 0; 
        if any(idxPostB)
            yB_post = movmean(yBeh(idxPostB), 5); sdB_post = movmean(sdBeh(idxPostB), 5);
            plotShadedSingle(ax, xBeh(idxPostB), yB_post, sdB_post, [0.6 0.6 0.6], 0.2);
        end
        yMax = max(yBeh + sdBeh, [], 'all', 'omitnan') * 1.3; if yMax<=0, yMax=1; end; ylim(ax, [0 yMax]);
    else
        ylim(ax, [0 1]);
    end
    set(ax, 'YColor', [0.4 0.4 0.4]);
    
    % LEFT AXIS: LAG (SOLID)
    yyaxis(ax, 'left'); hold(ax, 'on'); yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    idxPre = xLag <= 0; idxPost = xLag > 0;
    if any(idxPre), plotTube(ax, xLag(idxPre), yLag(idxPre), sdLag(idxPre), [0.2 0.2 0.8], 0.3); end
    if any(idxPost), plotTube(ax, xLag(idxPost), yLag(idxPost), sdLag(idxPost), [0.8 0.2 0.2], 0.3); end
    
    plotable = marks(marks >= xLims(1) & marks <= xLims(2));
    if ~isempty(plotable)
        yL = [-15 15]; yPos = repmat(yL(1)+1, size(plotable)); 
        % REDUCED MARKER SIZE TO 8
        plot(ax, plotable, yPos, 'v', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'Clipping', 'off');
    end
    title(ax, titleTxt, 'FontWeight', 'normal', 'FontSize', 10);
    xlim(ax, xLims); ylim(ax, [-15 15]); set(ax, 'YColor', 'k');
    xline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); box(ax, 'off'); hold(ax, 'off');
end

function plotTube(ax, x, y_line, y_std, col, alphaVal)
    % EXPLICITLY SOLID LINE
    x = x(:)'; y_line = y_line(:)'; y_std = y_std(:)';
    y_upper = y_line + y_std; y_lower = y_line - y_std;
    x_fill = [x, fliplr(x)]; y_fill = [y_upper, fliplr(y_lower)];
    nan_idx = isnan(y_fill); x_fill(nan_idx) = []; y_fill(nan_idx) = [];
    if ~isempty(x_fill), fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); end
    plot(ax, x, y_line, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '-'); 
end

function plotShadedSingle(ax, x, y, sd, col, alphaVal)
    % EXPLICITLY DASHED LINE
    x = x(:)'; y = y(:)'; sd = sd(:)'; 
    x_fill = [x, fliplr(x)]; y_fill = [y+sd, fliplr(y-sd)];
    nanIdx = isnan(y_fill); x_fill(nanIdx)=[]; y_fill(nanIdx)=[];
    if ~isempty(x_fill), fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); end
    plot(ax, x, y, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '--'); 
end

function [results] = calcLagAll(emgPreMean, emgPostMean, emgPostAll, emgPreAll, focusList)
    nFocus = length(focusList); nTotal = nFocus * nFocus; 
    emgPeakLagMList = cell(1, nTotal); emgPeakLagSDList = cell(1, nTotal);
    n=1; interpN = 200; calcTimeRange = round((0.5-0.15)*interpN):round((0.5+0.15)*interpN);
    for i = 1:nFocus 
        for preIndx = 1:nFocus 
            emgPre = emgPreMean(preIndx).data; 
            emgPreInterpM = mean(interp1(1:size(emgPre,1), emgPre, linspace(1, size(emgPre,1), interpN), 'linear'), 2);
            emgTargetPost = emgPostAll(i).dataAll;
            if ~isempty(emgTargetPost) && ~isempty(emgTargetPost{1})
                expDateN = length(emgTargetPost); maxTrials=0; for d=1:expDateN, maxTrials=max(maxTrials,size(emgTargetPost{d},2)); end
                lag_trials = nan(expDateN, maxTrials);
                for d=1:expDateN
                    if isempty(emgTargetPost{d}), continue; end
                    for t=1:size(emgTargetPost{d},2)
                         postS = interp1(1:size(emgTargetPost{d},1), emgTargetPost{d}(:,t), linspace(1,size(emgTargetPost{d},1),interpN), 'linear');
                         [c, lags] = xcorr(postS(calcTimeRange)'-mean(postS(calcTimeRange)), emgPreInterpM(calcTimeRange)-mean(emgPreInterpM(calcTimeRange)), 0.15*interpN, 'coeff');
                         [~, idx] = max(c); lag_trials(d, t) = lags(idx)/interpN*100;
                    end
                end
                lagPostM = nanmean(lag_trials, 2)'; lagPostSD = nanstd(lag_trials, 0, 2)';
            else
                lagPostM = NaN(1, length(emgPostMean(i).daysDiff)); lagPostSD = NaN(1, length(emgPostMean(i).daysDiff));
            end
            emgTargetPre = emgPreAll(i).dataAll;
            if ~isempty(emgTargetPre) && ~isempty(emgTargetPre{1})
                expDateN = length(emgTargetPre); maxTrials=0; for d=1:expDateN, maxTrials=max(maxTrials,size(emgTargetPre{d},2)); end
                lag_trials = nan(expDateN, maxTrials);
                for d=1:expDateN
                    if isempty(emgTargetPre{d}), continue; end
                    for t=1:size(emgTargetPre{d},2)
                         preS = interp1(1:size(emgTargetPre{d},1), emgTargetPre{d}(:,t), linspace(1,size(emgTargetPre{d},1),interpN), 'linear');
                         [c, lags] = xcorr(preS(calcTimeRange)'-mean(preS(calcTimeRange)), emgPreInterpM(calcTimeRange)-mean(emgPreInterpM(calcTimeRange)), 0.15*interpN, 'coeff');
                         [~, idx] = max(c); lag_trials(d, t) = lags(idx)/interpN*100;
                    end
                end
                lagPreM = nanmean(lag_trials, 2)'; lagPreSD = nanstd(lag_trials, 0, 2)';
            else
                lagPreM = zeros(1, length(emgPreMean(1).daysDiff)); lagPreSD = zeros(1, length(emgPreMean(1).daysDiff));
            end
            emgPeakLagMList{n}  = [lagPreM, lagPostM]; emgPeakLagSDList{n} = [lagPreSD, lagPostSD];
            n=n+1;
        end
    end
    results.daysDiffList = [emgPreMean(1).daysDiff(:)', emgPostMean(1).daysDiff(:)'];
    results.emgPeakLagMList = emgPeakLagMList; results.emgPeakLagSDList = emgPeakLagSDList;
end

function [dataMean, dataAll] = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    matPath  = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    if ~exist(matPath, 'file'), error('File not found: %s', matPath); end
    dataFile = load(matPath);
    dataCell = dataFile.synergyData;
    nList = length(focusList);
    focusListNum = zeros(1,nList);
    for i=1:nList
        match = find(contains(dataFile.nameList, focusList{i}));
        focusListNum(i) = match(1);
    end
    expDates = datetime(cellstr(dataFile.expDates), 'InputFormat','yyMMdd');
    daysDiff = days(expDates - ttDate);
    t = dataFile.percentTime;
    for i = 1:nList
        meanData = zeros(length(t), length(daysDiff));
        allData = cell(1, length(daysDiff));
        for dayI = 1:length(daysDiff)
            if ~isempty(dataCell{focusListNum(i),dayI})
                meanData(:,dayI) = mean(dataCell{focusListNum(i),dayI},1)';
                allData{dayI}    = (dataCell{focusListNum(i),dayI})';
            else
                meanData(:,dayI) = nan;
            end
        end
        dataMean(i).data = meanData; dataAll(i).dataAll = allData; dataMean(i).daysDiff = daysDiff; dataAll(i).daysDiff = daysDiff;
    end
end

function behData = loadBehavioralData(monkeyName, behDir)
    currentDir = pwd;
    try, cd(behDir); catch, behData.days=[]; behData.mean=[]; behData.std=[]; return; end
    if strcmpi(monkeyName, 'Yachimun'), f_days = dir('days.csv'); if isempty(f_days), f_days = dir('daysoriginal.csv'); end
    else, f_days = dir('days.csv'); end
    if isempty(f_days), cd(currentDir); behData.days=[]; behData.mean=[]; behData.std=[]; return; end
    try, daysRaw = csvread(f_days(1).name, 0, 0); catch, cd(currentDir); behData.days=[]; behData.mean=[]; behData.std=[]; return; end
    if strcmpi(monkeyName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(pat); [~, idx] = sort({files.name}); files = files(idx);
    c = {}; for k=1:length(files), if length(files(k).name)>=14, c{1,k}=files(k).name(1:14); else, c{1,k}=files(k).name; end, end
    uFiles = unique(c, 'stable');
    meanC = nan(1, length(uFiles)); stdC = nan(1, length(uFiles));
    for j=1:length(uFiles)
        mIdx = find(startsWith({files.name}, uFiles{j}), 1);
        if isempty(mIdx), continue; end
        try
            dat = csvread(files(mIdx).name, rOff, 1); dat=dat(:);
            meanC(j) = nanmean(dat); stdC(j) = nanstd(dat);
        catch, end
    end
    minLen = min(length(meanC), length(daysRaw));
    behData.days = daysRaw(1:minLen)'; behData.mean = meanC(1:minLen); behData.std = stdC(1:minLen);
    cd(currentDir);
end