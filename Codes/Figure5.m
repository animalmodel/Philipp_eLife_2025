% =========================================================================
% SCRIPT: PlotBehavioralRecovery.m
% AUTHOR: Roland Philipp
% =========================================================================

clear; close all; clc;
baseDir = 'C:\Users\mypre\Documents\Manuscripts\Revision\post acceptance revision\Philipp_eLife_2025';

%% ========================================================================
%  PART 1: MONKEY A
% =========================================================================
disp('Processing Monkey A...');

% --- PATHS ---


dir_A_Main   = fullfile(baseDir, 'Data', 'behavior', 'data_M1');
dir_A_Ataxia = fullfile(baseDir, 'Data', 'behavior', 'data_M1', 'real_mov_time - contains ATAXIA data');

% 1. Load Main Timeline
cd(dir_A_Main);
daysFile = dir('daysoriginal.csv');
days_A = csvread(daysFile.name, 0, 0);

% 2. Load Contact/Pull Data
csvFiles = dir('Ya*.csv');
fileNames = {csvFiles.name};
tempList = {};
for i = 1:length(fileNames), tempList{1,i} = fileNames{1,i}(1:14); end
fileList = char(unique(tempList));

touchData_A = [];
pullData_A  = [];

for j = 1:size(fileList, 1)
    dat = csvread(strtrim(fileList(j,:)), 3, 1);
    touchData_A(1:length(dat)/6, j) = dat(6:6:end);
    pullData_A(1:length(dat)/6, j)  = dat(5:6:end);
end

% 3. Plot Monkey A Figures (Hardcoded Indices: Pre=1:5, Post=6:end)
meanTouch_A = nanmean(touchData_A(1:20, :))'; 
stdTouch_A  = nanstd(touchData_A(1:20, :))';
plotData(days_A, meanTouch_A, stdTouch_A, 5, 'Monkey A: Contact Times', 'Time [ms]');

meanPull_A  = nanmean(pullData_A(1:20, :))';
stdPull_A   = nanstd(pullData_A(1:20, :))';
plotData(days_A, meanPull_A, stdPull_A, 5, 'Monkey A: Pull Times', 'Time [ms]');


%% ========================================================================
%  PART 2: MONKEY A - ATAXIA
% =========================================================================
disp('Processing Monkey A (Ataxia)...');
cd(dir_A_Ataxia);

% 1. Load Data
csvFiles = dir('Ya*.csv');
fileNames = {csvFiles.name};
tempList = {};
for i = 1:length(fileNames), tempList{1,i} = fileNames{1,i}(1:14); end
fileList = char(unique(tempList));

ataxiaData_A = [];
ataxiaStd_A  = [];

for j = 1:size(fileList, 1)
    dat = csvread(strtrim(fileList(j,:)), 3, 1);
    ataxiaData_A(j,1) = nanmean(dat);
    ataxiaStd_A(j,1)  = nanstd(dat);
end

% 2. Load Timeline (Try local first, else fallback to Main)
if exist('days.csv', 'file')
    days_A_Ataxia = csvread('days.csv', 0, 0);
else
    days_A_Ataxia = days_A; 
end

% 3. Plot (Indices: Pre=1:5, Post=6:end)
plotData(days_A_Ataxia, ataxiaData_A, ataxiaStd_A, 5, 'Monkey A: Ataxia', 'Time [ms]');


%% ========================================================================
%  PART 3: MONKEY B (ALL METRICS)
% =========================================================================
disp('Processing Monkey B...');
clearvars -except days_A baseDir; 

% --- PATHS ---
dir_B_Main     = fullfile(baseDir, 'Data', 'behavior', 'data_M2');
dir_B_Plate    = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'plateTouch');
dir_B_Aperture = fullfile(baseDir, 'Data', 'behavior', 'data_M2', 'aperture');

% 1. Load MAIN Timeline (Used for ALL Monkey B plots)
cd(dir_B_Main);
daysFile = dir('days.csv');
days_B = csvread(daysFile.name, 0, 0);

% 2. Load Main Data (Contact)
csvFiles = dir('Se*.csv');
fileNames = {csvFiles.name};
tempList = {};
for i = 1:length(fileNames), tempList{1,i} = fileNames{1,i}(1:14); end
fileList = char(unique(tempList));

touchData_B = [];
for j = 1:size(fileList, 1)
    dat = csvread(strtrim(fileList(j,:)), 3, 1);
    touchData_B(1:length(dat)/6, j) = dat(3:6:end);
end

% Plot Contact (Indices: Pre=1:3, Post=4:end)
meanTouch_B = nanmean(touchData_B(1:20, :))';
stdTouch_B  = nanstd(touchData_B(1:20, :))';
plotData(days_B, meanTouch_B, stdTouch_B, 3, 'Monkey B: Contact Times', 'Time [ms]');


% 3. Load Plate Data (Ataxia)
cd(dir_B_Plate);
csvFiles = dir('Se*.csv');
fileNames = {csvFiles.name};
tempList = {};
for i = 1:length(fileNames), tempList{1,i} = fileNames{1,i}(1:14); end
fileList = char(unique(tempList));

plateData_B = []; plateStd_B = [];
for j = 1:size(fileList, 1)
    try dat = csvread(strtrim(fileList(j,:)), 1, 1); catch, dat = csvread(strtrim(fileList(j,:)), 3, 1); end
    plateData_B(j,1) = nanmean(dat);
    plateStd_B(j,1)  = nanstd(dat);
end

% Plot Plate (Uses days_B from Main)
plotData(days_B, plateData_B, plateStd_B, 3, 'Monkey B: Ataxia (Plate)', 'Time [ms]');


% 4. Load Aperture Data
cd(dir_B_Aperture);
csvFiles = dir('Se*.csv');
fileNames = {csvFiles.name};
tempList = {};
for i = 1:length(fileNames), tempList{1,i} = fileNames{1,i}(1:14); end
fileList = char(unique(tempList));

apertureData_B = []; apertureStd_B = [];
for j = 1:size(fileList, 1)
    try dat = csvread(strtrim(fileList(j,:)), 1, 1); catch, dat = csvread(strtrim(fileList(j,:)), 3, 1); end
    apertureData_B(j,1) = nanmean(dat);
    apertureStd_B(j,1)  = nanstd(dat);
end

% Plot Aperture (Uses days_B from Main)
plotData(days_B, apertureData_B, apertureStd_B, 3, 'Monkey B: Aperture', 'Aperture [mm]');

disp('SUCCESS: All 6 figures generated.');


%% ========================================================================
%  HELPER: PLOT DATA (Robust & Matches Original Logic)
% =========================================================================
function plotData(days, meanData, stdData, preIdxEnd, titleStr, yLabelStr)
    figure('Name', titleStr, 'Color', 'w');
    
    % Force Row Vectors
    days = days(:)'; meanData = meanData(:)'; stdData = stdData(:)';
    
    % Match Lengths (Fixes "Vector Mismatch" error)
    minLen = min([length(days), length(meanData)]);
    days = days(1:minLen);
    meanData = meanData(1:minLen);
    stdData = stdData(1:minLen);
    
    % --- Pre-Surgery Plot ---
    x = days(1:preIdxEnd); 
    y = meanData(1:preIdxEnd); 
    sd = stdData(1:preIdxEnd);
    
    fill([x, fliplr(x)], [y+sd, fliplr(y-sd)], 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none'); hold on;
    plot(x, y, 'k', 'LineWidth', 2);
    
    % --- Post-Surgery Plot ---
    if length(days) > preIdxEnd
        x = days(preIdxEnd+1:end); 
        y = meanData(preIdxEnd+1:end); 
        sd = stdData(preIdxEnd+1:end);
        
        % Smoothing (Only on Post, as in original)
        y = movmean(y, 5); 
        sd = movmean(sd, 5);
        
        fill([x, fliplr(x)], [y+sd, fliplr(y-sd)], 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        plot(x, y, 'k', 'LineWidth', 2);
    end
    
    title(titleStr); xlabel('Days relative to Surgery'); ylabel(yLabelStr);
    grid on; box off; set(gca, 'TickDir', 'out');
end