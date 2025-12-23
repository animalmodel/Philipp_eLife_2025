% =========================================================================
% SCRIPT: Plot_Figure10_SynergyCD_v44.m
%
% PURPOSE: 
%   Generates Figure 10: Synergy C & D Cross-Correlation Analysis.
%   Layout: 4 Rows x 2 Columns (Portrait Mode).
%
% COMPARISONS:
%   - Top Row (Blue): Post-surgery Syn C/D vs Pre-surgery Synergy B (Extensor).
%   - Bottom Row (Red): Post-surgery Syn C/D vs Pre-surgery Synergy A (Flexor).
%   - Overlay: Behavioral Recovery (Ataxia/Contact Time) on Right Axis.
%
% AUTHOR: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
baseDir = 'C:\Users\mypre\Documents\Manuscripts\Revision\post acceptance revision\Philipp_eLife_2025';

% Input Directories
matDir    = fullfile(baseDir, 'Data', 'synergy'); 
behDirA   = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
behDirB   = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

% Output Directory
fullfile(baseDir, 'outputFigures_Fig10');

% Settings
calcType = 'Synergy'; 
monkeyNameList = {'Yachimun', 'Seseki'};

% Landmark Days
LMD_A = [29, 64, 69, 79, 99];
LMD_B = [22, 36, 44, 48, 64];

% Axis Ticks
ticksA = [-20, 0, 40, 80, 120];
ticksB = [0, 20, 40, 60];

% Colors
colorTop = [0.2 0.2 0.8]; % Blue (vs Syn B - Extensor)
colorBot = [0.8 0.2 0.2]; % Red  (vs Syn A - Flexor)

%% 2. MAIN ANALYSIS LOOP
% -------------------------------------------------------------------------
monkey_data_All  = cell(length(monkeyNameList), 1);
monkey_data_Mean = cell(length(monkeyNameList), 1);
behavior_data    = cell(length(monkeyNameList), 1);

for mIndx = 1:length(monkeyNameList)
    monkeyName = monkeyNameList{mIndx};
    fprintf('\n--- Processing %s ---\n', monkeyName);
    
    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        behPath = behDirA;
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        behPath = behDirB;
    end
    focusList = {'A', 'B', 'C', 'D'};
   
    % Load Data
    [dataPreMean,  dataPreAll]  = loadMatData(calcType, monkeyName, 'PreAll',  matDir, focusList, ttDate);
    [dataPostMean, dataPostAll] = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);
    
    % Calculate Correlations
    [monkey_data_All{mIndx}] = calcCorrAll(dataPreMean, dataPostMean, dataPostAll, dataPreAll, focusList);
    [monkey_data_Mean{mIndx}] = calcCorrMean_modified(dataPreMean, dataPostMean, focusList);
    
    % Load Behavior
    [behavior_data{mIndx}] = loadBehavioralData(monkeyName, behPath);
end

%% 3. PLOTTING SECTION (Portrait Layout)
% -------------------------------------------------------------------------
fprintf('\n--- Generating Figure 10 ---\n');

scrsz = get(0,'ScreenSize');
figPos = [100, 100, scrsz(3)*0.33, scrsz(4)*0.85]; % Portrait Aspect Ratio
fig = figure('Name', 'Figure 10: Synergy C/D Correlation', 'Position', figPos, 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

for mIndx = 1:2
    data_All  = monkey_data_All{mIndx};
    data_Mean = monkey_data_Mean{mIndx};
    dataBeh   = behavior_data{mIndx};
    
    if mIndx == 1
        rowOffset = 0; % Monkey A = Rows 1-2
        currLMD = LMD_A; currXTicks = ticksA; xLimVal = [-45 120];
        lbl = 'Monkey A';
    else
        rowOffset = 4; % Monkey B = Rows 3-4
        currLMD = LMD_B; currXTicks = ticksB; xLimVal = [-10 70];
        lbl = 'Monkey B';
    end
    
    % --- ROW 1 (Top) - BLUE - vs Pre-B (Extensor) ---
    
    % Col 1: Synergy C vs Pre-B -> Pair [2, 3] (Pre 2=B, Post 3=C)
    ax1 = nexttile(1 + rowOffset);
    pairIdx = find(data_All.corrPair(:,1)==2 & data_All.corrPair(:,2)==3); 
    plotOverlay(ax1, data_All, data_Mean, pairIdx, dataBeh, colorTop, [-1 1], currLMD, currXTicks, xLimVal);
    if mIndx == 1, title(ax1, 'Synergy C', 'FontSize', 14, 'FontWeight', 'bold'); end
    text(ax1, 0.05, 0.9, lbl, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 12);
    text(ax1, 0.5, 0.2, 'vs Syn B (pre)', 'Units', 'normalized', 'Color', colorTop, 'FontSize', 9);
    
    % Col 2: Synergy D vs Pre-B -> Pair [2, 4] (Pre 2=B, Post 4=D)
    ax2 = nexttile(2 + rowOffset);
    pairIdx = find(data_All.corrPair(:,1)==2 & data_All.corrPair(:,2)==4); 
    plotOverlay(ax2, data_All, data_Mean, pairIdx, dataBeh, colorTop, [-1 1], currLMD, currXTicks, xLimVal);
    if mIndx == 1, title(ax2, 'Synergy D', 'FontSize', 14, 'FontWeight', 'bold'); end
    
    % --- ROW 2 (Bottom) - RED - vs Pre-A (Flexor) ---
    
    % Col 1: Synergy C vs Pre-A -> Pair [1, 3] (Pre 1=A, Post 3=C)
    ax3 = nexttile(3 + rowOffset);
    pairIdx = find(data_All.corrPair(:,1)==1 & data_All.corrPair(:,2)==3); 
    plotOverlay(ax3, data_All, data_Mean, pairIdx, dataBeh, colorBot, [-1 1], currLMD, currXTicks, xLimVal);
    text(ax3, 0.5, 0.2, 'vs Syn A (pre)', 'Units', 'normalized', 'Color', colorBot, 'FontSize', 9);
    
    % Col 2: Synergy D vs Pre-A -> Pair [1, 4] (Pre 1=A, Post 4=D)
    ax4 = nexttile(4 + rowOffset);
    pairIdx = find(data_All.corrPair(:,1)==1 & data_All.corrPair(:,2)==4); 
    plotOverlay(ax4, data_All, data_Mean, pairIdx, dataBeh, colorBot, [-1 1], currLMD, currXTicks, xLimVal);
    text(ax4, 0.5, 0.2, 'vs Syn A (pre)', 'Units', 'normalized', 'Color', colorBot, 'FontSize', 9);
    
    % Formatting
    ylabel(ax1, 'Correlation'); ylabel(ax3, 'Correlation'); 
    
    if mIndx == 2
        xlabel(ax3, 'Days Post-TT'); xlabel(ax4, 'Days Post-TT');
    else
        set(ax3, 'XTickLabel', []); set(ax4, 'XTickLabel', []);
    end
    set(ax1, 'XTickLabel', []); set(ax2, 'XTickLabel', []);
    set(ax2, 'YTickLabel', []); set(ax4, 'YTickLabel', []);
    
    yyaxis(ax2, 'right'); ylabel(ax2, 'Contact Time (ms)');
    yyaxis(ax4, 'right'); ylabel(ax4, 'Contact Time (ms)');
end

disp('Figure 10 generated successfully.');

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotOverlay(ax, data_All, data_Mean, idx, dataBeh, mainColor, yRange, landmarkDays, xTickValues, xLimVal)
    daysDiff = data_All.daysDiffList;
    y_line   = data_Mean.emgCorrMList{idx};
    y_std    = data_All.emgCorrSDList{idx};
    
    preEnd    = data_All.dayPreEnd;
    postStart = data_All.dayPostStart;
    
    % LEFT AXIS: CORRELATION
    yyaxis(ax, 'left'); hold(ax, 'on'); set(ax, 'YColor', 'k');
    preIdx = daysDiff <= preEnd; 
    postIdx = daysDiff >= postStart;
    
    plotTube(ax, daysDiff(preIdx), y_line(preIdx), y_std(preIdx), mainColor, 0.3);
    plotTube(ax, daysDiff(postIdx), y_line(postIdx), y_std(postIdx), mainColor, 0.3);
    
    plot(ax, daysDiff(preIdx), y_line(preIdx), 'Color', mainColor, 'LineWidth', 2);
    plot(ax, daysDiff(postIdx), y_line(postIdx), 'Color', mainColor, 'LineWidth', 2);
    
    ylim(ax, yRange); set(ax, 'YTick', -1:0.5:1);
    
    % RIGHT AXIS: BEHAVIOR
    yyaxis(ax, 'right'); set(ax, 'YColor', [0.2 0.2 0.2]);
    if ~isempty(dataBeh.days)
        bx = dataBeh.days; by = dataBeh.mean; bsd = dataBeh.std;
        preB = bx < 0; postB = bx >= 0;
        
        if any(postB), by(postB) = movmean(by(postB), 5); bsd(postB) = movmean(bsd(postB), 5); end
        
        if any(preB), plotShaded(ax, bx(preB), by(preB), bsd(preB), [0.1 0.1 0.1], 0.2); end
        if any(postB), plotShaded(ax, bx(postB), by(postB), bsd(postB), [0.1 0.1 0.1], 0.2); end
        
        maxY = max(by + bsd, [], 'all', 'omitnan'); if isnan(maxY)||maxY==0, maxY=100; end
        ylim(ax, [0, maxY*1.3]);
    else
        ylim(ax, [0 1]);
    end
    
    % Formatting
    xlim(ax, xLimVal);
    xline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    xline(ax, postStart, ':', 'Color', [0.6 0.6 0.6]);
    yline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    
    set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'Clipping', 'off');
    box(ax, 'off');
    
    % Landmarks
    yyaxis(ax, 'left');
    validLMD = landmarkDays(landmarkDays >= xLimVal(1) & landmarkDays <= xLimVal(2));
    if ~isempty(validLMD)
        plot(ax, validLMD, repmat(yRange(1), size(validLMD)), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'Clipping', 'off');
    end
    set(ax, 'XTick', xTickValues); hold(ax, 'off');
end

function plotTube(ax, x, y, sd, col, alphaVal)
    x = x(:)'; y = y(:)'; sd = sd(:)';
    fill(ax, [x, fliplr(x)], [y+sd, fliplr(y-sd)], col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function plotShaded(ax, x, y, sd, col, alphaVal)
    x = x(:)'; y = y(:)'; sd = sd(:)';
    fill(ax, [x, fliplr(x)], [y+sd, fliplr(y-sd)], col, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
    plot(ax, x, y, '--', 'Color', col, 'LineWidth', 1.5);
end

% --- DATA LOADING ---
function [dataMean, dataAll] = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    fName = sprintf('synergyData_%s_%s.mat', monkeyName, condition);
    fullPath = fullfile(matDir, fName);
    if ~exist(fullPath, 'file'), error('File not found: %s', fullPath); end
    
    D = load(fullPath, 'synergyData', 'percentTime', 'expDates', 'nameList');
    daysDiff = days(datetime(cellstr(D.expDates), 'InputFormat', 'yyMMdd') - ttDate);
    daysDiff = double(daysDiff(:)'); % Force Row Vector
    
    % Map Focus List Indices
    indices = zeros(1, length(focusList));
    for i=1:length(focusList), indices(i) = find(contains(D.nameList, focusList{i}), 1); end
    
    for i=1:length(focusList)
        idx = indices(i);
        nDays = length(daysDiff);
        meanMat = zeros(length(D.percentTime), nDays);
        allTrials = cell(1, nDays);
        
        for d=1:nDays
            raw = D.synergyData{idx, d};
            if ~isempty(raw)
                meanMat(:, d) = mean(raw, 1)';
                allTrials{d}  = raw';
            end
        end
        dataMean(i).data = meanMat; dataMean(i).daysDiff = daysDiff;
        dataAll(i).dataAll = allTrials; dataAll(i).daysDiff = daysDiff;
    end
end

function behData = loadBehavioralData(monkeyName, behDir)
    behData.days=[]; behData.mean=[]; behData.std=[];
    if ~exist(behDir, 'dir'), return; end
    
    dFile = dir(fullfile(behDir, 'days*.csv'));
    if isempty(dFile), return; end
    daysRaw = csvread(fullfile(behDir, dFile(1).name), 0, 0);
    daysRaw = daysRaw(:)';
    
    if strcmp(monkeyName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(fullfile(behDir, pat));
    
    [~, idx] = sort({files.name}); files = files(idx);
    c = {}; for k=1:length(files), c{1,k} = files(k).name(1:min(14, end)); end
    uFiles = unique(c, 'stable');
    
    meanC = nan(1, length(uFiles)); stdC = nan(1, length(uFiles));
    for i=1:length(uFiles)
        fIdx = find(startsWith({files.name}, uFiles{i}), 1);
        if isempty(fIdx), continue; end
        try
            dat = csvread(fullfile(behDir, files(fIdx).name), rOff, 1);
            if strcmp(monkeyName, 'Yachimun'), dat(dat==0)=NaN; end
            meanC(i) = nanmean(dat); stdC(i) = nanstd(dat);
        catch, end
    end
    
    minLen = min(length(meanC), length(daysRaw));
    [behData.days, sortIdx] = sort(daysRaw(1:minLen));
    behData.mean = meanC(sortIdx); behData.std = stdC(sortIdx);
end

% --- CALCULATIONS ---
function results = calcCorrAll(emgPreMean, emgPostMean, emgPostAll, emgPreAll, focusList)
    nFocus = length(focusList);
    nTotal = nFocus^2;
    interpN = 200; calcRange = round(0.35*interpN):round(0.65*interpN);
    
    emgCorrMList = cell(1, nTotal); emgCorrSDList = cell(1, nTotal); corrPair = zeros(nTotal, 2); 
    counter = 1;
    
    for postIdx = 1:nFocus
        for preIdx = 1:nFocus
            corrPair(counter,:) = [preIdx, postIdx];
            
            % Ref Signal (Pre Mean)
            preData = emgPreMean(preIdx).data;
            preInterp = interp1(1:size(preData,1), preData, linspace(1, size(preData,1), interpN), 'linear');
            refSig = mean(preInterp, 2); refSig = refSig(calcRange) - mean(refSig(calcRange));
            
            calcTrialCorrs = @(target) compute_trial_corrs(target, refSig, interpN, calcRange);
            [mPost, sPost] = calcTrialCorrs(emgPostAll(postIdx));
            [mPre, sPre]   = calcTrialCorrs(emgPreAll(postIdx));
            
            emgCorrMList{counter}  = [mPre, mPost];
            emgCorrSDList{counter} = [sPre, sPost];
            counter = counter + 1;
        end
    end
    
    d1 = emgPreMean(1).daysDiff; d2 = emgPostMean(1).daysDiff;
    results.daysDiffList = [d1(:)', d2(:)'];
    results.emgCorrMList = emgCorrMList; results.emgCorrSDList = emgCorrSDList;
    results.corrPair = corrPair;
    results.dayPreEnd = max(d1); results.dayPostStart = min(d2);
end

function [mVal, sVal] = compute_trial_corrs(targetStruct, refSig, interpN, calcRange)
    dataAll = targetStruct.dataAll; nDays = length(dataAll);
    corrMat = nan(nDays, 50);
    for d=1:nDays
        if isempty(dataAll{d}), continue; end
        trials = dataAll{d}; nTrials = size(trials, 2);
        lenT = size(trials, 1); tInterp = linspace(1, lenT, interpN);
        for t=1:nTrials
            sig = interp1(1:lenT, trials(:,t), tInterp, 'linear')';
            sig = sig(calcRange) - mean(sig(calcRange));
            [c, lags] = xcorr(refSig, sig, 0.15*interpN, 'coeff');
            corrMat(d, t) = c(lags==0);
        end
    end
    mVal = nanmean(corrMat, 2)'; sVal = nanstd(corrMat, 0, 2)';
end

function results = calcCorrMean_modified(emgPre, emgPost, focusList)
    nFocus = length(focusList); interpN = 200; calcRange = round(0.35*interpN):round(0.65*interpN);
    emgCorrList = cell(1, nFocus^2); counter = 1;
    
    for postIdx = 1:nFocus
        for preIdx = 1:nFocus
            preRef = interp1(1:size(emgPre(preIdx).data,1), emgPre(preIdx).data, linspace(1, size(emgPre(preIdx).data,1), interpN));
            refSig = mean(preRef, 2); refSig = refSig(calcRange) - mean(refSig(calcRange));
            
            doCorr = @(target) compute_mean_corr(target, refSig, interpN, calcRange);
            emgCorrList{counter} = [doCorr(emgPre(postIdx).data), doCorr(emgPost(postIdx).data)];
            counter = counter + 1;
        end
    end
    results.emgCorrMList = emgCorrList;
end

function corrs = compute_mean_corr(targetData, refSig, interpN, calcRange)
    if isempty(targetData), corrs = []; return; end
    len = size(targetData, 1); tInterp = interp1(1:len, targetData, linspace(1, len, interpN));
    corrs = zeros(1, size(targetData, 2));
    for i=1:size(targetData, 2)
        sig = tInterp(:, i); sig = sig(calcRange) - mean(sig(calcRange));
        [c, lags] = xcorr(refSig, sig, 0.15*interpN, 'coeff');
        corrs(i) = c(lags==0);
    end
end