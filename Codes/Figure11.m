% =========================================================================
% SCRIPT: Figure11.m (Plot_aaEMG_with_Bars)
%
% Description:
%   - Generates Figure 11: 4x4 Grid.
%   - Pairs 'aaEMG Time Course' with 'Statistical Bar Plots' for each synergy.
%   - Implements logic from 'AccuEMG_new.m' for bar graphs.
%
% AUTHOR: Roland Philipp
% =========================================================================
clear; clc; close all;

%% 1. Configuration
% -------------------------------------------------------------------------
% --- DYNAMIC PATH SETUP ---
scriptPath = fileparts(mfilename('fullpath'));
if isempty(scriptPath), scriptPath = pwd; end % Fallback for running sections
baseDir = fileparts(scriptPath); 

fprintf('Detected Base Directory: %s\n', baseDir);

% --- PATHS ---
matDirA = fullfile(baseDir, 'Data', 'emg', 'aggregated_EMG_data', 'M1');
matDirB = fullfile(baseDir, 'Data', 'emg', 'aggregated_EMG_data', 'M2');
behMainA = fullfile(baseDir, 'Data', 'behavior', 'data_M1');
behSubA  = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');
behMainB = fullfile(baseDir, 'Data', 'behavior', 'data_M2');
behSubB  = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');

% Verify Paths
if ~exist(matDirA, 'dir'), error('Data folder missing: %s', matDirA); end

% --- LANDMARK DAYS (Indices for Bar Plots) ---
% Derived from AccuEMG_new.m
LMD_Idx_A = [4, 5, 22, 25, 30, 42]; 
Labels_A  = {'Pre', '29', '64', '69', '79', '99'};

LMD_Idx_B = [3, 5, 12, 15, 17, 25];
Labels_B  = {'Pre', '22', '36', '44', '48', '64'};

% Bar Plot Colors
barColor = [0.6 0.6 0.6]; 

% Muscle Synergies
grpA = { [3, 4], [10, 11, 12], [5, 6], [1, 7, 8] }; 
grpB = { [8, 11], [1, 2, 3], [5, 12, 6], [4] };
synNames = {'Synergy A', 'Synergy B', 'Synergy C', 'Synergy D'};

%% 2. Load Data
% -------------------------------------------------------------------------
fprintf('--- Loading Data ---\n');

% Monkey A
daysA = loadTimeline(behMainA, 'daysoriginal.csv');
dataA = load_aaEMG('Yachimun', matDirA, grpA);
behA  = loadBehavioralData('Yachimun', behSubA, daysA);

% Monkey B
daysB = loadTimeline(behMainB, 'days.csv');
dataB = load_aaEMG('Seseki', matDirB, grpB);
behB  = loadBehavioralData('Seseki', behSubB, daysB);

%% 3. Plotting (4x4 Grid)
% -------------------------------------------------------------------------
fig = figure('Name', 'Figure 11: aaEMG + Bars', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Color', 'w');
t = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:4
    % --- MONKEY A (Left Side) ---
    
    % 1. Time Course
    ax1 = nexttile;
    plot_aaEMG_Panel(ax1, dataA, behA, i, ['Monkey A: ' synNames{i}], [-45 120], 1.0);
    
    % 2. Bar Plot
    ax2 = nexttile;
    plot_BarStats(ax2, dataA.synStats{i}, LMD_Idx_A, Labels_A, barColor);
    
    % --- MONKEY B (Right Side) ---
    
    % 3. Time Course
    ax3 = nexttile;
    plot_aaEMG_Panel(ax3, dataB, behB, i, ['Monkey B: ' synNames{i}], [-10 70], 0.6);
    
    % 4. Bar Plot
    ax4 = nexttile;
    plot_BarStats(ax4, dataB.synStats{i}, LMD_Idx_B, Labels_B, barColor);
    
    % --- Axis Formatting ---
    if i < 4
        xlabel(ax1,''); xlabel(ax2,''); xlabel(ax3,''); xlabel(ax4,'');
    end
end

disp('Figure Generated.');

%% ========================================================================
%  FUNCTIONS
% =========================================================================

function plot_BarStats(ax, statsStruct, lmdIndices, labels, barCol)
    % Extracts data for specific Landmark Days
    means = statsStruct.mean(lmdIndices);
    sems  = statsStruct.sem(lmdIndices);
    
    % Plot
    b = bar(ax, 1:length(means), means, 'FaceColor', barCol, 'EdgeColor', 'none');
    hold(ax, 'on');
    errorbar(ax, 1:length(means), means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Formatting
    set(ax, 'XTick', 1:length(means), 'XTickLabel', labels, 'TickDir', 'out', 'FontSize', 9);
    xtickangle(ax, 45);
    ylabel(ax, 'Mean EMG [\muV]');
    xlim(ax, [0.5 length(means)+0.5]);
    
    % Match Y-limits to reasonable range
    maxY = max(means + sems);
    if isnan(maxY) || maxY==0, maxY=1; end
    ylim(ax, [0 maxY*1.2]);
    box(ax, 'off');
end

function plot_aaEMG_Panel(ax, emgData, behData, synIdx, titleStr, xLimVal, widthFactor)
    % Left Axis (EMG)
    yyaxis(ax, 'left'); set(ax, 'YColor', 'k');
    x = emgData.days; 
    y = emgData.synStats{synIdx}.mean; % Use Mean
    
    plot(ax, x, y, 'k-', 'LineWidth', 1.5); hold(ax, 'on');
    scatter(ax, x, y, 20, 'k', 'filled');
    
    xlim(ax, xLimVal); 
    yMax = max(y, [], 'omitnan'); if isempty(yMax), yMax=1; end
    ylim(ax, [0 yMax*1.2]);
    
    % Right Axis (Behavior)
    yyaxis(ax, 'right'); set(ax, 'YColor', [0.5 0.5 0.5]);
    if ~isempty(behData.days)
        bx = behData.days; by = behData.mean;
        postIdx = bx >= 0;
        if any(postIdx), by(postIdx) = movmean(by(postIdx), 5); end
        plot(ax, bx, by, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        
        bMax = max(by, [], 'omitnan'); if isempty(bMax), bMax=1; end
        ylim(ax, [0 bMax*1.5]);
    end
    
    xline(ax, 0, 'k:');
    title(ax, titleStr, 'FontSize', 10);
    box(ax, 'off');
end

function data = load_aaEMG(monkeyName, dirPath, groups)
    % Robust loading logic
    if strcmpi(monkeyName, 'Yachimun')
        fname = fullfile(dirPath, 'alignedEMG_data.mat');
        varName = 'Ptrig2';
    else
        fname = fullfile(dirPath, 'alignedEMG_data_Seseki.mat');
        varName = 'Ptrig3';
    end
    
    if ~exist(fname, 'file'), error('File not found: %s', fname); end
    D = load(fname);
    P = D.(varName);
    
    % Use loaded x-axis or generate one
    if isfield(D, 'xpostdays'), days = D.xpostdays; else, days = 1:size(P.plotData_sel, 1); end
    
    % Define Window (-15 to 15)
    if isfield(P, 'x')
        wIdx = find(P.x >= -15 & P.x <= 15);
    else
        wIdx = 1:size(P.plotData_sel{1},2); 
    end
    
    synStats = cell(1, 4);
    for g=1:4
        muscles = groups{g};
        numSessions = size(P.plotData_sel, 1);
        
        dailyMean = nan(1, numSessions);
        dailySEM  = nan(1, numSessions);
        
        for s=1:numSessions
            if isempty(P.plotData_sel{s}), continue; end
            
            % 1. Extract Muscles
            M = P.plotData_sel{s}(muscles, :);
            
            % 2. Sum (Aggregated EMG)
            aggM = sum(M, 1);
            
            % 3. Window
            windowData = aggM(wIdx);
            
            % 4. Stats
            dailyMean(s) = mean(windowData);
            dailySEM(s)  = std(windowData) / sqrt(length(windowData));
        end
        
        synStats{g}.mean = dailyMean;
        synStats{g}.sem  = dailySEM;
    end
    
    data.days = days;
    data.synStats = synStats;
end

function days = loadTimeline(dirPath, fileName)
    fullPath = fullfile(dirPath, fileName);
    if exist(fullPath, 'file')
        days = csvread(fullPath, 0, 0);
    else
        d = dir(fullfile(dirPath, 'days*.csv'));
        if ~isempty(d), days = csvread(fullfile(dirPath, d(1).name), 0, 0);
        else, days = []; end
    end
end

function behData = loadBehavioralData(monkeyName, behDir, daysVector)
    currentDir = pwd;
    try cd(behDir); catch, behData.days=[]; behData.mean=[]; return; end
    
    if strcmpi(monkeyName, 'Yachimun'), pat='Ya*.csv'; rOff=3; else, pat='Se*.csv'; rOff=1; end
    files = dir(pat);
    
    % Sort & Unique
    [~, idx] = sort({files.name}); files = files(idx);
    c = {}; for k=1:length(files), c{1,k}=files(k).name(1:14); end
    uFiles = unique(c, 'stable');
    
    meanC = nan(1, length(uFiles));
    for j=1:length(uFiles)
        mIdx = find(startsWith({files.name}, uFiles{j}), 1);
        if isempty(mIdx), continue; end
        try
            dat = csvread(files(mIdx).name, rOff, 1);
            meanC(j) = nanmean(dat(:));
        catch, end
    end
    
    minLen = min(length(meanC), length(daysVector));
    behData.days = daysVector(1:minLen); behData.mean = meanC(1:minLen);
    
    % Force Row Vectors
    behData.days = behData.days(:)'; behData.mean = behData.mean(:)';
    
    cd(currentDir);
end