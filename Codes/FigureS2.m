% =========================================================================
% SCRIPT: FigureS2.m
%
% PURPOSE: 
%   Generates Supplementary Figure S2: EMG Patterns.
%   - Left Cols: Pre-surgery mean +/- SD (shaded tube).
%   - Right Cols: Post-surgery individual days (color-coded).
%
% UPDATES:
%   - Fixed Monkey A missing FDS/EDC (used 'FDSdist'/'EDCdist').
%   - Removed Post-Surgery panels for FDS, FCU, and FCR in Monkey B.
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
baseDir = fileparts(scriptPath); 

% Force calculation type to EMG
calcType = 'EMG';

% Set data directory
matDir = fullfile(baseDir, 'Data', 'emg', 'emg_mat');

% Set output directory
outFigDir = fullfile(scriptPath, 'outputFigures_S2');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify data directory exists
if ~exist(matDir, 'dir')
    error('Data directory not found: %s', matDir);
end

monkeyNameList = {'Yachimun', 'Seseki'};

% Define muscle lists
% NOTE: Monkey A uses 'dist' suffix internally for FDS/EDC
focusList_A = {'FDSdist', 'EDCdist', 'FDP', 'ED23', 'FCU', 'ECU', 'FCR', 'ECR', 'PL', 'BRD'};
focusList_B = {'FDS', 'EDC', 'FDP', 'ED23', 'FCU', 'ED45', 'FCR', 'ECU', 'PL', 'ECR', 'Deltoid', 'BRD'};

% Plotting colors
colPreLine = [0.2, 0.2, 0.6];  % Dark purple/blue for mean line
colPreFill = [0.4, 0.4, 0.8];  % Lighter purple/blue for shading

%% 2. MAIN LOOP (MONKEYS)
% -------------------------------------------------------------------------
for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);

    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        focusList = focusList_A;
        figTitle = 'Monkey A';
        cmapName = 'hot'; 
        cbarLabel = 'Post surgery days';
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        focusList = focusList_B;
        figTitle = 'Monkey B';
        cmapName = 'summer'; 
        cbarLabel = 'Post surgery days';
    end

    % --- Load Data (All Days) ---
    dataPre = loadMatData(calcType, monkeyName, 'PreAll', matDir, focusList, ttDate);
    dataPost = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);

    % --- Figure Setup ---
    nMuscles = length(focusList);
    nRows = ceil(nMuscles / 2);
    figHeight = nRows * 4; 

    f = figure('Name', ['Figure S2 - ' figTitle], 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 20, figHeight]);
    t = tiledlayout(nRows, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Colormap Setup for Post-Surgery ---
    allPostDays = [];
    for i = 1:nMuscles
        allPostDays = [allPostDays, dataPost(i).daysDiff];
    end
    if isempty(allPostDays)
        dayRange = [0 1]; 
    else
        dayRange = [min(allPostDays), max(allPostDays)];
    end
    colormap(cmapName);
    cMap = colormap;
    nColors = size(cMap, 1);

    % --- Plotting Loop (Muscle Pairs) ---
    for row = 1:nRows
        for colPair = 1:2
            muscIdx = (row - 1) * 2 + colPair;
            if muscIdx > nMuscles, continue; end
            
            % Get Muscle Name & Clean Label (remove 'dist')
            rawName = focusList{muscIdx};
            cleanName = strrep(rawName, 'dist', '');
            
            % Add asterisk to FDS and EDC
            if contains(cleanName, 'FDS') || strcmp(cleanName, 'EDC')
                plotLabel = ['*' cleanName];
            else
                plotLabel = cleanName;
            end

            % =============================================================
            % 1. Pre-Surgery Plot (Mean +/- SD)
            % =============================================================
            axPre = nexttile((row - 1) * 4 + colPair);
            hold(axPre, 'on');
            
            if ~isempty(dataPre(muscIdx).dataM)
                meanData = mean(dataPre(muscIdx).dataM, 2, 'omitnan');
                sdData = std(dataPre(muscIdx).dataM, 0, 2, 'omitnan');
                timeVec = dataPre(muscIdx).t;
                plotShadedError(axPre, timeVec, meanData, sdData, colPreLine, colPreFill);
            end
            
            % Formatting
            xline(axPre, 0, 'k--');
            xlim(axPre, [-15, 15]);
            ylabel(axPre, plotLabel, 'FontWeight', 'bold', 'FontSize', 9);
            if row < nRows, set(axPre, 'XTickLabel', []); end
            if colPair == 2, set(axPre, 'YTickLabel', []); ylabel(axPre, ''); title(axPre, plotLabel, 'FontWeight', 'bold', 'FontSize', 9); end
            if row == 1 && colPair == 1
                 text(axPre, 0.5, 1.15, figTitle, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            end

            % =============================================================
            % 2. Post-Surgery Plot (Individual Days)
            % =============================================================
            axPost = nexttile((row - 1) * 4 + colPair + 2);
            hold(axPost, 'on');

            % Monkey B Exception: Skip plotting FDS, FCU, and FCR in Post-Surgery
            skipPlot = false;
            if mIndx == 2
                if contains(cleanName, 'FDS') || contains(cleanName, 'FCU') || contains(cleanName, 'FCR')
                    skipPlot = true;
                end
            end

            if skipPlot
                % Hide axis completely for skipped plots
                axis(axPost, 'off');
            else
                if ~isempty(dataPost(muscIdx).dataM)
                    numDays = size(dataPost(muscIdx).dataM, 2);
                    for k = 1:numDays
                        dayVal = dataPost(muscIdx).daysDiff(k);
                        colorIdx = round(((dayVal - dayRange(1)) / (dayRange(2) - dayRange(1))) * (nColors - 1)) + 1;
                        colorIdx = max(1, min(nColors, colorIdx)); 
                        thisColor = cMap(colorIdx, :);
                        
                        plot(axPost, dataPost(muscIdx).t, dataPost(muscIdx).dataM(:, k), '-', 'Color', thisColor, 'LineWidth', 1);
                    end
                end
                 % Formatting
                xline(axPost, 0, 'k--');
                xlim(axPost, [-15, 15]);
                set(axPost, 'YTickLabel', []);
                if row < nRows, set(axPost, 'XTickLabel', []); end
                title(axPost, plotLabel, 'FontWeight', 'bold', 'FontSize', 9);
            end

            % Axis Labels (Bottom Row Only)
            if row == nRows && colPair == 1
                xlabel(axPre, 'Task range [%]', 'FontSize', 9);
                ylabel(axPre, {'Amplitude [µV]', plotLabel}, 'FontWeight', 'bold', 'FontSize', 9);
            end
            if row == nRows && colPair == 2
                 xlabel(axPost, 'Task range [%]', 'FontSize', 9);
            end

            % Colorbar (Last plot only)
            if muscIdx == nMuscles
                c = colorbar(axPost);
                c.Label.String = cbarLabel;
                cbarLimits = caxis(axPost); 
                caxis(axPost, dayRange); 
            end
        end
    end
    
    % Save figure
    exportgraphics(f, fullfile(outFigDir, [monkeyName '_FigureS2.png']), 'Resolution', 300);
    fprintf('Figure saved to: %s\n', fullfile(outFigDir, [monkeyName '_FigureS2.png']));
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotShadedError(ax, x, yMean, ySD, lineCol, fillCol)
    x = x(:)'; yMean = yMean(:)'; ySD = ySD(:)';
    yUpper = yMean + ySD;
    yLower = yMean - ySD;
    xFill = [x, fliplr(x)];
    yFill = [yUpper, fliplr(yLower)];
    fill(ax, xFill, yFill, fillCol, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(ax, x, yMean, '-', 'Color', lineCol, 'LineWidth', 1.5);
end

function dataAll = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    matPath = fullfile(matDir, ['emgData_' monkeyName '_' condition '.mat']);
    
    if ~exist(matPath, 'file')
         warning('Data file not found: %s', matPath);
         dataAll = struct('t', [], 'dataM', [], 'daysDiff', []);
         dataAll = repmat(dataAll, length(focusList), 1);
         return;
    end

    dataFile = load(matPath, 'emgData', 'percentTime', 'expDates', 'nameList');
    dataCell = dataFile.emgData;
    nList = length(focusList);
    dataAll = struct([]);

    for i = 1:nList
        % Robust Search: Match focusList name to nameList in file
        targetName = focusList{i};
        muscleIdx = find(strcmp(targetName, dataFile.nameList));
        
        if isempty(muscleIdx)
            % Try 'dist' fallback if not found (e.g. if user asked for EDC but file has EDCdist)
            muscleIdx = find(strcmp([targetName 'dist'], dataFile.nameList));
        end

        if isempty(muscleIdx)
            % Warning suppressed to avoid clutter if user intentionally requests missing muscle
            dataAll(i).t = dataFile.percentTime;
            dataAll(i).dataM = [];
            dataAll(i).daysDiff = [];
            continue;
        end

        t = dataFile.percentTime;
        expDayString = cellstr(dataFile.expDates);
        expDates = datetime(expDayString, 'InputFormat', 'yyMMdd');
        daysDiff = days(expDates - ttDate);

        numDays = length(expDates);
        muscleDataM = zeros(length(t), numDays);
        validDays = false(1, numDays);

        for dayI = 1:numDays
            dailyTrials = dataCell{muscleIdx, dayI};
            if ~isempty(dailyTrials)
                muscleDataM(:, dayI) = mean(dailyTrials, 1, 'omitnan')';
                validDays(dayI) = true;
            end
        end
        
        dataAll(i).t = t;
        dataAll(i).dataM = muscleDataM(:, validDays);
        dataAll(i).daysDiff = daysDiff(validDays);
    end
end