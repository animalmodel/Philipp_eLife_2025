% =========================================================================
% SCRIPT: FigureS1.m
%
% PURPOSE: 
%   Generates Figure S1: EMG Lag Analysis (Peak Correlation Timing).
%   Layout: Fixed 4x2 Grid (Monkey A Top, Monkey B Bottom).
%
% STYLING:
%   - Left Axis (Lag): Solid lines (Blue=Pre, Red=Post).
%   - Right Axis (Behavior): Dashed Gray Envelope.
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
matDir  = fullfile(baseDir, 'Data', 'emg', 'emg_mat'); 
% Note: Updated to include 'behavior' folder to match other scripts
behDirA = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
behDirB = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

% Output Directory
outFigDir = fullfile(baseDir, 'outputFigures_FigS1');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify Data
if ~exist(matDir, 'dir')
    error('Data folder not found at: %s\n(Did you download the ''Data'' folder from GitHub?)', matDir); 
end

% Settings
calcType = 'EMG'; 
monkeyNameList = {'Yachimun', 'Seseki'};

%% 2. MAIN ANALYSIS LOOP
% -------------------------------------------------------------------------
monkey_data = cell(length(monkeyNameList), 1);

for mIndx = 1:length(monkeyNameList)
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);
    
    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        focusList = {'EDCdist', 'FDSdist', 'FCR'}; 
        behPath = behDirA;
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        focusList = {'EDC', 'FDS', 'ECR'}; 
        behPath = behDirB;
    end
   
    % Load Data
    [dataPreMean,  dataPreAll]  = loadMatData(calcType, monkeyName, 'PreAll',  matDir, focusList, ttDate);
    [dataPostMean, dataPostAll] = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);
    
    % Calculations
    [monkey_data{mIndx}.lags] = calcLagAll(dataPreMean, dataPostMean, dataPostAll, dataPreAll, focusList);
    [monkey_data{mIndx}.beh]  = loadBehavioralData(monkeyName, behPath);
end

%% 3. PLOTTING SECTION
% -------------------------------------------------------------------------
landmarkDaysA = [29, 64, 69, 79, 99];
landmarkDaysB = [22, 36, 44, 48, 64];

% A4 Dimensions (Vertical 13.5 x 27 cm)
fig = figure('Name', 'Figure S1: EMG Lag', 'Units', 'centimeters', ...
             'Position', [2 1 13.5 27], 'Color', 'w');

t = tiledlayout(4, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

% Comparison Titles
titlesA = {'Post-FDS vs Pre-EDC', 'Post-EDC vs Pre-EDC', 'Post-FDS vs Pre-FDS', 'Post-EDC vs Pre-FDS'};
plotOrderA = [4, 1, 5, 2];

titlesB = {'', 'Post-EDC vs Pre-EDC', '', 'Post-EDC vs Pre-FDS'};
plotOrderB = [NaN, 1, NaN, 2]; 

for mIndx = 1:2
    D = monkey_data{mIndx};
    
    if mIndx == 1
        marks    = landmarkDaysA;
        rowStart = 0; 
        xLimVal  = [-45 120];
        orders   = plotOrderA;
        txts     = titlesA;
        behLbl   = 'Off-Target (s)';
    else
        marks    = landmarkDaysB;
        rowStart = 4; 
        xLimVal  = [-10 70];
        orders   = plotOrderB;
        txts     = titlesB;
        behLbl   = 'Contact Time (s)';
    end
    
    for i = 1:4
        idx = orders(i);
        ax = nexttile(rowStart + i);
        
        if isnan(idx)
            axis(ax, 'off'); 
            continue;
        end
        
        plotDualAxis(ax, D.lags.daysDiffList, ...
                     D.lags.emgPeakLagMList{idx}, D.lags.emgPeakLagSDList{idx}, ...
                     D.beh.days, D.beh.mean, D.beh.std, ...
                     txts{i}, marks, xLimVal, behLbl);
                 
        % Formatting
        if mIndx == 1 && i <= 2, set(ax, 'XTickLabel', []); end 
        
        if mod(i,2) == 1 
            % Left Column
            yyaxis(ax, 'left'); ylabel(ax, 'Lag (%)', 'FontWeight', 'bold', 'Color', 'k');
            yyaxis(ax, 'right'); set(ax, 'YTickLabel', []); ylabel(ax, '');
        else 
            % Right Column
            yyaxis(ax, 'left'); set(ax, 'YTickLabel', []); ylabel(ax, '');
            yyaxis(ax, 'right'); ylabel(ax, behLbl, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4]);
        end
    end
end

% Save
try
    print(fig, fullfile(outFigDir, 'FigureS1_EMGLag.svg'), '-dsvg', '-painters');
    fprintf('Figure S1 saved to: %s\n', outFigDir);
catch
    exportgraphics(fig, fullfile(outFigDir, 'FigureS1_EMGLag.png'), 'Resolution', 300);
    fprintf('SVG save failed. Saved PNG to: %s\n', outFigDir);
end


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotDualAxis(ax, xLag, yLag, sdLag, xBeh, yBeh, sdBeh, titleTxt, marks, xLims, behLbl)
    % RIGHT AXIS: BEHAVIOR (Dashed Gray)
    yyaxis(ax, 'right'); hold(ax, 'on');
    hasBeh = ~isempty(xBeh) && any(~isnan(yBeh));
    
    if hasBeh
        idxPreB = xBeh < 0; 
        if any(idxPreB), plotShadedSingle(ax, xBeh(idxPreB), yBeh(idxPreB), sdBeh(idxPreB), [0.6 0.6 0.6], 0.2); end
        
        idxPostB = xBeh >= 0; 
        if any(idxPostB)
            yB_post = movmean(yBeh(idxPostB), 5); sdB_post = movmean(sdBeh(idxPostB), 5);
            plotShadedSingle(ax, xBeh(idxPostB), yB_post, sdB_post, [0.6 0.6 0.6], 0.2);
        end
        
        yMax = max(yBeh + sdBeh, [], 'all', 'omitnan') * 1.3; 
        if yMax <= 0, yMax = 1; end
        ylim(ax, [0 yMax]);
    else
        ylim(ax, [0 1]);
    end
    set(ax, 'YColor', [0.4 0.4 0.4]);
    
    % LEFT AXIS: LAG (Solid)
    yyaxis(ax, 'left'); hold(ax, 'on');
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    idxPre = xLag <= 0; idxPost = xLag > 0;
    if any(idxPre), plotTube(ax, xLag(idxPre), yLag(idxPre), sdLag(idxPre), [0.2 0.2 0.8], 0.3); end
    if any(idxPost), plotTube(ax, xLag(idxPost), yLag(idxPost), sdLag(idxPost), [0.8 0.2 0.2], 0.3); end
    
    % Landmarks
    plotable = marks(marks >= xLims(1) & marks <= xLims(2));
    if ~isempty(plotable)
        yPos = repmat(-14, size(plotable)); 
        plot(ax, plotable, yPos, 'v', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'Clipping', 'off');
    end
    
    title(ax, titleTxt, 'FontWeight', 'normal', 'FontSize', 10);
    xlim(ax, xLims); ylim(ax, [-15 15]); set(ax, 'YColor', 'k');
    xline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    box(ax, 'off'); hold(ax, 'off');
end

function plotTube(ax, x, y_line, y_std, col, alphaVal)
    % Plots Solid Line with Error Tube
    x = x(:)'; y_line = y_line(:)'; y_std = y_std(:)';
    y_upper = y_line + y_std; y_lower = y_line - y_std;
    x_fill = [x, fliplr(x)]; y_fill = [y_upper, fliplr(y_lower)];
    
    if ~isempty(x_fill)
        fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    end
    plot(ax, x, y_line, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '-'); 
end

function plotShadedSingle(ax, x, y, sd, col, alphaVal)
    % Plots Dashed Line with Error Shade
    x = x(:)'; y = y(:)'; sd = sd(:)'; 
    x_fill = [x, fliplr(x)]; y_fill = [y+sd, fliplr(y-sd)];
    
    if ~isempty(x_fill)
        fill(ax, x_fill, y_fill, col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    end
    plot(ax, x, y, 'Color', col, 'LineWidth', 1.5, 'LineStyle', '--'); 
end

% --- CALCULATIONS ---
function [results] = calcLagAll(emgPreMean, emgPostMean, emgPostAll, emgPreAll, focusList)
    nFocus = length(focusList); nTotal = nFocus^2; 
    emgPeakLagMList = cell(1, nTotal); emgPeakLagSDList = cell(1, nTotal);
    
    n = 1; interpN = 200; 
    calcTimeRange = round((0.5-0.15)*interpN):round((0.5+0.15)*interpN);
    
    for i = 1:nFocus 
        for preIndx = 1:nFocus 
            emgPre = emgPreMean(preIndx).data; 
            emgPreInterpM = mean(interp1(1:size(emgPre,1), emgPre, linspace(1, size(emgPre,1), interpN), 'linear'), 2);
            refSig = emgPreInterpM(calcTimeRange) - mean(emgPreInterpM(calcTimeRange));
            
            calcLags = @(targetStruct) compute_lags(targetStruct, refSig, interpN, calcTimeRange);
            
            [lagPostM, lagPostSD] = calcLags(emgPostAll(i));
            [lagPreM, lagPreSD]   = calcLags(emgPreAll(i));
            
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
    if isempty(dataAll) || isempty(dataAll{1}), lagM=[]; lagSD=[]; return; end
    
    nDays = length(dataAll);
    maxTrials = 0; for d=1:nDays, if ~isempty(dataAll{d}), maxTrials=max(maxTrials,size(dataAll{d},2)); end; end
    
    lag_trials = nan(nDays, maxTrials);
    for d = 1:nDays
        if isempty(dataAll{d}), continue; end
        trials = dataAll{d};
        nTrials = size(trials, 2);
        lenT = size(trials, 1); tInterp = linspace(1, lenT, interpN);
        
        for t = 1:nTrials
            trialSig = interp1(1:lenT, trials(:,t), tInterp, 'linear')';
            sigSeg = trialSig(calcTimeRange) - mean(trialSig(calcTimeRange));
            
            [c, lags] = xcorr(sigSeg, refSig, 0.15*interpN, 'coeff');
            [~, maxIdx] = max(c); 
            lag_trials(d, t) = lags(maxIdx) / interpN * 100;
        end
    end
    lagM = nanmean(lag_trials, 2)'; lagSD = nanstd(lag_trials, 0, 2)';
end

% --- DATA LOADING ---
function [dataMean, dataAll] = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    fName = sprintf('emgData_%s_%s.mat', monkeyName, condition);
    matPath = fullfile(matDir, fName);
    
    if ~exist(matPath, 'file'), warning('File not found: %s', matPath); dataMean=[]; dataAll=[]; return; end
    
    D = load(matPath);
    dataCell = D.emgData;
    nList = length(focusList);
    
    indices = zeros(1,nList);
    for i = 1:nList
        idx = find(contains(D.nameList, focusList{i}), 1);
        if ~isempty(idx), indices(i) = idx; end
    end
    
    expDates = datetime(cellstr(D.expDates), 'InputFormat','yyMMdd');
    daysDiff = days(expDates - ttDate);
    t = D.percentTime;
    
    for i = 1:nList
        idx = indices(i);
        meanData = zeros(length(t), length(daysDiff));
        allData = cell(1, length(daysDiff));
        
        for d = 1:length(daysDiff)
            if idx>0 && ~isempty(dataCell{idx, d})
                meanData(:,d) = mean(dataCell{idx, d}, 1)';
                allData{d}    = dataCell{idx, d}';
            else
                meanData(:,d) = nan;
            end
        end
        dataMean(i).data = meanData; dataAll(i).dataAll = allData;
        dataMean(i).daysDiff = daysDiff; dataAll(i).daysDiff = daysDiff;
    end
end

function behData = loadBehavioralData(monkeyName, behDir)
    behData.days=[]; behData.mean=[]; behData.std=[];
    if ~exist(behDir, 'dir'), return; end
    
    dFile = dir(fullfile(behDir, 'days*.csv'));
    if isempty(dFile), return; end
    daysRaw = csvread(fullfile(behDir, dFile(1).name), 0, 0);
    
    if strcmp(monkeyName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(fullfile(behDir, pat));
    [~, idx] = sort({files.name}); files = files(idx);
    
    c = {}; for k=1:length(files), c{1,k} = files(k).name(1:min(14,end)); end
    uFiles = unique(c, 'stable');
    
    meanC = nan(1, length(uFiles)); stdC = nan(1, length(uFiles));
    for j = 1:length(uFiles)
        mIdx = find(startsWith({files.name}, uFiles{j}), 1);
        if isempty(mIdx), continue; end
        try
            dat = csvread(fullfile(behDir, files(mIdx).name), rOff, 1); dat=dat(:);
            meanC(j) = nanmean(dat); stdC(j) = nanstd(dat);
        catch, end
    end
    minLen = min(length(meanC), length(daysRaw));
    behData.days = daysRaw(1:minLen)'; behData.mean = meanC(1:minLen); behData.std = stdC(1:minLen);
end