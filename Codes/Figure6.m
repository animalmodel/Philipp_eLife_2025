% =========================================================================
% SCRIPT: Figure6.m
%
% PURPOSE: 
% Potting figure 6: EMG Profiles and Cross-Correlation Analysis
%
% UPDATED:
%   - Dynamic path detection for GitHub portability.
%   - Monkey B Cross-Corr: Plotted in Column 4 (Stacked).
%   - Post-Surgery EMG: Plots ONLY Landmark Days.
%   - Behavioral lines disconnected (Pre/Post).
%   - Y-Axis limits fixed (0 to Inf).
%
% AUTHORS: Roland Philipp 
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
% --- DYNAMIC PATH SETUP ---
% Automatically detects the folder this script is in (e.g., .../Codes)
% and sets the base directory to the parent folder (e.g., .../Philipp_eLife_2025)
scriptPath = fileparts(mfilename('fullpath'));
if isempty(scriptPath), scriptPath = pwd; end % Fallback for running sections
baseDir = fileparts(scriptPath); 

fprintf('Detected Base Directory: %s\n', baseDir);

% --- Paths ---
dirEMG     = fullfile(baseDir, 'Data', 'emg','emg_xls');
dirCorrMat = fullfile(baseDir, 'Data', 'emg','emg_mat');


dirBehA = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
dirBehB = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

outFigDir = fullfile(baseDir, 'outputFigures_Fig6');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% --- Settings ---
monkeyNameList = {'Yachimun', 'Seseki'};
calcType = 'EMG'; 

% Landmark Days (Used for both filtering EMG and plotting triangles)
LMD_A = [31, 64, 69, 79, 99];
LMD_B = [22, 36, 44, 48, 64];

%% 2. INITIALIZE FIGURE
% -------------------------------------------------------------------------
fig = figure('Name', 'Figure 6: Complete Analysis', 'Color', 'w', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.7 0.9]);

t = tiledlayout(6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 3. MAIN LOOP
% -------------------------------------------------------------------------
for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('--- Processing %s ---\n', monkeyName);
    
    if mIndx == 1
        % --- MONKEY A CONFIG ---
        emgList = {'EDCdist', 'FDSdist', 'FCR'};
        ttDate  = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        behPath = dirBehA;
        corrOrder = [1, 5, 2, 4]; 
        corrTitles = {'EDC pre/ FDS post', 'EDC pre/ EDC post', ...
                      'FDS pre/ FDS post', 'FDS pre/ EDC post'};
        currLMD = LMD_A;
        rowStart = 0; 
    else
        % --- MONKEY B CONFIG ---
        emgList = {'EDC', 'FDS', 'ECR'};
        ttDate  = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        behPath = dirBehB;
        corrOrder = [1, 3]; 
        corrTitles = {'EDC pre/ EDC post', 'FDS pre/ EDC post'};
        currLMD = LMD_B;
        rowStart = 12; 
    end
    
    % --- PART A: LEFT SIDE (EMG Profiles) ---
    [CPre, CPost] = loadEMGProfiles(monkeyName, dirEMG, emgList, ttDate);
    
    for i = 1:3
        r = i - 1; 
        
        % Pre-Surgery (Col 1)
        axPre = nexttile(rowStart + (r*4) + 1);
        plotEMG_Pre(axPre, CPre{i}, emgList{i});
        if i == 1 && mIndx == 1, title(axPre, 'Monkey A', 'FontSize', 12, 'FontWeight', 'bold'); end
        if i == 1 && mIndx == 2, title(axPre, 'Monkey B', 'FontSize', 12, 'FontWeight', 'bold'); end
        
        % Post-Surgery (Col 2)
        axPost = nexttile(rowStart + (r*4) + 2);
        
        % Skip FDS Post for Monkey B
        if mIndx == 2 && i == 2
            axis(axPost, 'off');
        else
            % Pass currLMD to filter plots
            plotEMG_Post(axPost, CPost{i}, emgList{i}, currLMD);
        end
    end
    
    % --- PART B: RIGHT SIDE (Cross Correlation) ---
    [dPreM, dPreAll]   = loadCorrData(calcType, monkeyName, 'PreAll', dirCorrMat, emgList, ttDate);
    [dPostM, dPostAll] = loadCorrData(calcType, monkeyName, 'PostAll', dirCorrMat, emgList, ttDate);
    
    resCorr_All  = calcCorrAll(dPreM, dPostM, dPostAll, dPreAll, emgList);
    resCorr_Mean = calcCorrMean(dPreM, dPostM, emgList);
    resBeh       = loadBehavioralData(monkeyName, behPath);
    
    nCorrPlots = length(corrOrder);
    for k = 1:nCorrPlots
        pIdx = corrOrder(k);
        
        if mIndx == 1
            % Monkey A: 2x2 Grid in Cols 3-4
            if k <= 2, tileIdx = rowStart + 2 + k; else, tileIdx = rowStart + 6 + (k-2); end
            axCorr = nexttile(tileIdx);
        else
            % Monkey B: Stacked in Col 4
            % Row 4 start=12 -> Tile 16 is R4,C4
            % Row 5 start=16 -> Tile 20 is R5,C4
            if k == 1
                tileIdx = 16; 
            else
                tileIdx = 20; 
            end
            axCorr = nexttile(tileIdx); 
        end
        
        plotOverlay(axCorr, resCorr_All.daysDiffList, ...
            resCorr_Mean.emgCorrMList{pIdx}, ...
            resCorr_All.emgCorrSDList{pIdx}, ...
            resBeh, ...
            corrTitles{k}, [-1 1], resCorr_All.dayPreEnd, resCorr_All.dayPostStart, currLMD);
            
        if mIndx == 1 && k <= 2, set(axCorr, 'XTickLabel', []); end
        if mIndx == 2 && k == 1, set(axCorr, 'XTickLabel', []); end
        
        % Y-Axis Labels logic
        if mIndx == 1 && mod(k,2) == 0
             yyaxis(axCorr, 'left'); set(axCorr, 'YTickLabel', []); ylabel('');
        end

    end
end

% Save
saveas(fig, fullfile(outFigDir, 'Figure6_Combined.fig'));
exportgraphics(fig, fullfile(outFigDir, 'Figure6_Combined.png'), 'Resolution', 300);
try
    print(fig, fullfile(outFigDir, 'Figure6_Combined.svg'), '-dsvg', '-painters');
    fprintf('Figure 6 saved to: %s\n', outFigDir);
catch
    fprintf('SVG save failed. Figure available in PNG.\n');
end


%% ========================================================================
%  HELPER FUNCTIONS: PLOTTING
% =========================================================================
function plotEMG_Pre(ax, dataStruct, mName)
    hold(ax, 'on');
    plot(ax, dataStruct.x, dataStruct.meanData, 'k-', 'LineWidth', 1.5); 
    xline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
    xlim(ax, [-15 15]); 
    ylim(ax, [0, inf]); 
    ylabel(ax, mName, 'FontWeight', 'bold', 'FontSize', 9);
    set(ax, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off');
end

function plotEMG_Post(ax, dataStruct, mName, LMD)
    hold(ax, 'on');
    
    % Filter for Landmark Days
    availableDays = str2double(dataStruct.days); % Convert string list to numbers
    
    % Find indices of available days that are in the Landmark list
    [~, idxInStruct, ~] = intersect(availableDays, LMD);
    
    if isempty(idxInStruct)
        warning('No matching Landmark Days found for %s', mName);
    end
    
    % Plot only selected days
    numSel = length(idxInStruct);
    colors = parula(numSel);
    
    for k = 1:numSel
        dIdx = idxInStruct(k);
        plot(ax, dataStruct.x, dataStruct.dailyMeans(:, dIdx), '-', 'Color', colors(k, :), 'LineWidth', 1.2);
    end
    
    xline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
    xlim(ax, [-15 15]); 
    ylim(ax, [0, inf]); 
    set(ax, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off', 'YTickLabel', []);
end

function plotOverlay(ax, daysDiff, meanLine, tubeSD, dataBeh, titleStr, yRange, preEnd, postStart, LMD)
    preIdx = daysDiff <= preEnd; postIdx = daysDiff >= postStart;
    
    % Left Axis (Correlation)
    yyaxis(ax, 'left'); hold(ax, 'on'); set(ax, 'YColor', 'k');
    plotTube(ax, daysDiff(preIdx), meanLine(preIdx), tubeSD(preIdx), [0.2 0.2 0.8], 0.3);
    plotTube(ax, daysDiff(postIdx), meanLine(postIdx), tubeSD(postIdx), [0.8 0.2 0.2], 0.3);
    plot(ax, daysDiff(preIdx), meanLine(preIdx), 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5);
    plot(ax, daysDiff(postIdx), meanLine(postIdx), 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
    ylim(ax, yRange); ylabel(ax, 'Corr Coeff', 'FontSize', 8);
    
    % Right Axis (Behavior)
    yyaxis(ax, 'right'); set(ax, 'YColor', [0.3 0.3 0.3]);
    
    bx = dataBeh.days(:)'; by = dataBeh.mean(:)'; bsd = dataBeh.std(:)';
    if ~isempty(bx)
        idxPreB  = bx < 0; idxPostB = bx >= 0;
        
        if any(idxPostB)
            by(idxPostB) = movmean(by(idxPostB), 5); 
            bsd(idxPostB) = movmean(bsd(idxPostB), 5); 
        end
        
        if any(idxPreB)
            fill(ax, [bx(idxPreB), fliplr(bx(idxPreB))], ...
                 [by(idxPreB)+bsd(idxPreB), fliplr(by(idxPreB)-bsd(idxPreB))], ...
                 [0.1 0.1 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax, bx(idxPreB), by(idxPreB), '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
        end
        
        if any(idxPostB)
            fill(ax, [bx(idxPostB), fliplr(bx(idxPostB))], ...
                 [by(idxPostB)+bsd(idxPostB), fliplr(by(idxPostB)-bsd(idxPostB))], ...
                 [0.1 0.1 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax, bx(idxPostB), by(idxPostB), '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1);
        end
        ylim(ax, [0 max(by+bsd, [], 'all')*1.3]); ylabel(ax, 'Time (s)', 'FontSize', 8);
    end
    
    xlim(ax, [-20 120]); title(ax, titleStr, 'FontSize', 10, 'FontWeight', 'normal');
    set(ax, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off');
    
    yyaxis(ax, 'left');
    validLMD = LMD(LMD >= -20 & LMD <= 120);
    plot(ax, validLMD, repmat(-1, size(validLMD)), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end

function plotTube(ax, x, y, sd, col, alpha)
    x=x(:)'; y=y(:)'; sd=sd(:)'; 
    fill(ax, [x, fliplr(x)], [y+sd, fliplr(y-sd)], col, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% --- LOADING FUNCTIONS ---
function [CPre, CPost] = loadEMGProfiles(monkeyName, dirPath, emgList, ttDate)
    nEMG = length(emgList);
    CPre = cell(nEMG, 1); CPost = cell(nEMG, 1);
    
    fname = fullfile(dirPath, ['emgData_' monkeyName '_PreAll.xlsx']);
    if ~exist(fname, 'file'), error('EMG File not found: %s', fname); end
    for i = 1:nEMG
        T = readtable(fname, 'Sheet', emgList{i}, 'VariableNamingRule', 'preserve');
        CPre{i}.x = T.TimeSample;
        CPre{i}.meanData = mean(T{:, 2:end}, 2, 'omitnan');
    end
    
    fname = fullfile(dirPath, ['emgData_' monkeyName '_PostAll.xlsx']);
    for i = 1:nEMG
        T = readtable(fname, 'Sheet', emgList{i}, 'VariableNamingRule', 'preserve');
        x = T.TimeSample; data = T{:, 2:end};
        expDays = cellstr(regexprep(string(T.Properties.VariableNames(2:end)), '^(\d{4})-(\d{2})(\d{2})$', '$1-$2-$3'));
        
        % Calculate Days relative to TTDate for matching with Landmarks
        dateObjs = datetime(expDays, 'InputFormat', 'yyyy-MM-dd');
        dayDiffs = days(dateObjs - ttDate);
        
        uDays = unique(dayDiffs, 'stable'); % Unique Numeric Days
        
        dMeans = zeros(length(x), length(uDays));
        for d = 1:length(uDays)
            % Match numeric days back to data columns
            dMeans(:, d) = mean(data(:, dayDiffs == uDays(d)), 2, 'omitnan');
        end
        CPost{i}.x = x; 
        CPost{i}.dailyMeans = dMeans;
        CPost{i}.days = string(uDays); % Store days as string list for easy matching
    end
end

function [dMean, dAll] = loadCorrData(calcType, mName, cond, dirPath, fList, ttDate)
    fName = fullfile(dirPath, ['emgData_' mName '_' cond '.mat']);
    if ~exist(fName, 'file'), error('Corr File missing: %s', fName); end
    D = load(fName); 
    indices = zeros(1, length(fList));
    for i=1:length(fList), indices(i) = find(contains(D.nameList, fList{i}), 1); end
    
    daysDiff = double(days(datetime(cellstr(D.expDates), 'InputFormat', 'yyMMdd') - ttDate));
    daysDiff = daysDiff(:)'; 
    
    for i=1:length(fList)
        dMean(i).data = []; dAll(i).dataAll = [];
        for d=1:length(daysDiff)
            raw = D.emgData{indices(i), d};
            if ~isempty(raw)
                dMean(i).data(:,d) = mean(raw, 1)';
                dAll(i).dataAll{d} = raw';
            end
        end
        dMean(i).daysDiff = daysDiff; dAll(i).daysDiff = daysDiff;
    end
end

function behData = loadBehavioralData(mName, dirPath)
    behData.days=[]; behData.mean=[]; behData.std=[];
    if ~exist(dirPath, 'dir'), return; end
    dFile = dir(fullfile(dirPath, 'days*.csv'));
    if isempty(dFile), return; end
    daysRaw = csvread(fullfile(dirPath, dFile(1).name), 0, 0);
    
    if strcmp(mName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(fullfile(dirPath, pat));
    [~, idx] = sort({files.name}); files = files(idx);
    
    c = {}; for k=1:length(files), c{1,k} = files(k).name(1:min(14,end)); end
    uFiles = unique(c, 'stable');
    
    meanC = []; stdC = [];
    for j=1:length(uFiles)
        fIdx = find(startsWith({files.name}, uFiles{j}), 1);
        if isempty(fIdx), continue; end
        try
            dat = csvread(fullfile(dirPath, files(fIdx).name), rOff, 1);
            if strcmp(mName, 'Yachimun'), dat(dat==0)=NaN; end
            meanC(j) = nanmean(dat(:)); stdC(j) = nanstd(dat(:));
        catch, end
    end
    
    if ~isempty(meanC)
        minLen = min(length(daysRaw), length(meanC));
        behData.days = daysRaw(1:minLen)'; 
        behData.mean = meanC(1:minLen);    
        behData.std  = stdC(1:minLen);     
    end
end

function results = calcCorrAll(dPreM, dPostM, dPostAll, dPreAll, fList)
    nF = length(fList); nTot = nF^2; cnt = 1; interpN = 200;
    cRange = round(0.35*interpN):round(0.65*interpN);
    
    for postI = 1:nF
        for preI = 1:nF
            preD = dPreM(preI).data;
            preInterp = interp1(1:size(preD,1), preD, linspace(1, size(preD,1), interpN));
            refSig = mean(preInterp, 2); refSig = refSig(cRange) - mean(refSig(cRange));
            
            calcStats = @(S) getCorrStats(S, refSig, interpN, cRange);
            [mPost, sPost] = calcStats(dPostAll(postI));
            [mPre, sPre]   = calcStats(dPreAll(postI));
            
            res.emgCorrMList{cnt} = [mPre, mPost];
            res.emgCorrSDList{cnt} = [sPre, sPost];
            res.corrPair(cnt,:) = [preI, postI];
            cnt = cnt + 1;
        end
    end
    
    d1 = dPreM(1).daysDiff; d2 = dPostM(1).daysDiff;
    res.daysDiffList = [d1(:)', d2(:)'];
    
    res.dayPreEnd = max(d1); res.dayPostStart = min(d2);
    results = res;
end

function [mVal, sVal] = getCorrStats(targetStruct, refSig, interpN, cRange)
    dAll = targetStruct.dataAll;
    vals = nan(length(dAll), 50);
    for d=1:length(dAll)
        if isempty(dAll{d}), continue; end
        trials = dAll{d};
        for t=1:size(trials,2)
            sig = interp1(1:size(trials,1), trials(:,t), linspace(1,size(trials,1),interpN))';
            sig = sig(cRange) - mean(sig(cRange));
            [c, lags] = xcorr(refSig, sig, 0.15*interpN, 'coeff');
            vals(d,t) = c(lags==0);
        end
    end
    mVal = nanmean(vals, 2)'; sVal = nanstd(vals, 0, 2)';
end

function results = calcCorrMean(dPreM, dPostM, fList)
    nF = length(fList); cnt = 1; interpN = 200; cRange = round(0.35*interpN):round(0.65*interpN);
    for postI = 1:nF
        for preI = 1:nF
            preD = dPreM(preI).data;
            refSig = mean(interp1(1:size(preD,1), preD, linspace(1, size(preD,1), interpN)), 2);
            refSig = refSig(cRange) - mean(refSig(cRange));
            
            doC = @(D) getMeanCorr(D, refSig, interpN, cRange);
            results.emgCorrMList{cnt} = [doC(dPreM(postI).data), doC(dPostM(postI).data)];
            cnt = cnt + 1;
        end
    end
end

function vals = getMeanCorr(data, refSig, interpN, cRange)
    if isempty(data), vals = []; return; end
    interpD = interp1(1:size(data,1), data, linspace(1, size(data,1), interpN));
    vals = zeros(1, size(data,2));
    for i=1:size(data,2)
        sig = interpD(:,i); sig = sig(cRange) - mean(sig(cRange));
        [c, lags] = xcorr(refSig, sig, 0.15*interpN, 'coeff');
        vals(i) = c(lags==0);
    end
end