% =========================================================================
% SCRIPT: FigureS7.m
%
% PURPOSE: 
%   Generates Supplementary Figure S7: Time varying activation profiles 
%   of muscle synergies.
%
% LAYOUT (4 Rows x 4 Cols):
%   - Rows: Synergies A, B, C, D
%   - Col 1: Monkey A Pre (Mean +/- SD)
%   - Col 2: Monkey A Post (All days, color-coded)
%   - Col 3: Monkey B Pre (Mean +/- SD)
%   - Col 4: Monkey B Post (All days, color-coded)
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
baseDir = fileparts(scriptPath); 

% Force calculation type to Synergy
calcType = 'Synergy';

% Set data directory
matDir = fullfile(baseDir, 'Data', 'synergy');

% Set output directory
outFigDir = fullfile(scriptPath, 'outputFigures_S7');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify data directory exists
if ~exist(matDir, 'dir')
    error('Data directory not found: %s', matDir);
end

monkeyNameList = {'Yachimun', 'Seseki'};
focusList = {'A', 'B', 'C', 'D'}; % Synergies

% Plotting colors
colPreLine = [0.2, 0.2, 0.6];  % Dark purple/blue for mean line
colPreFill = [0.4, 0.4, 0.8];  % Lighter purple/blue for shading

% Figure Setup
figWidth_cm = 24; figHeight_cm = 20; 
f = figure('Name', 'Figure S7: Synergy Activation', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, figWidth_cm, figHeight_cm]);
T = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 2. MAIN LOOP
% -------------------------------------------------------------------------
for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('Processing %s...\n', monkeyName);

    if mIndx == 1
        ttDate = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        figTitle = 'Monkey A';
        % Custom Colormap: Red (early) to Black (late)
        cMapBase = [1 0 0]; % Start Red
        cMapEnd  = [0.2 0 0]; % End Dark Red/Black
        colOffset = 0; % Plot in cols 1 & 2
        cbarLabel = 'Post surgery days';
    else
        ttDate = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        figTitle = 'Monkey B';
        % Custom Colormap: Green (early) to Black (late)
        cMapBase = [0 1 0]; % Start Green
        cMapEnd  = [0 0.2 0]; % End Dark Green/Black
        colOffset = 2; % Plot in cols 3 & 4
        cbarLabel = 'Post surgery days';
    end

    % --- Load Data ---
    dataPre = loadMatData(calcType, monkeyName, 'PreAll', matDir, focusList, ttDate);
    dataPost = loadMatData(calcType, monkeyName, 'PostAll', matDir, focusList, ttDate);

    % --- Colormap Setup for Post-Surgery ---
    allPostDays = [];
    for i = 1:length(focusList)
        allPostDays = [allPostDays, dataPost(i).daysDiff];
    end
    if isempty(allPostDays), dayRange = [0 1]; else, dayRange = [min(allPostDays), max(allPostDays)]; end
    
    % Generate Gradient Colormap
    nColors = 64;
    cMap = [linspace(cMapBase(1), cMapEnd(1), nColors)', ...
            linspace(cMapBase(2), cMapEnd(2), nColors)', ...
            linspace(cMapBase(3), cMapEnd(3), nColors)'];

    % --- Plotting Loop (Rows = Synergies) ---
    for row = 1:4
        synIdx = row;
        synLabel = ['Synergy ' focusList{synIdx}];
        
        % =============================================================
        % 1. Pre-Surgery Plot (Mean +/- SD) - Column 1 (or 3)
        % =============================================================
        axPre = nexttile((row - 1) * 4 + colOffset + 1);
        hold(axPre, 'on');
        
        if ~isempty(dataPre(synIdx).dataM)
            meanData = mean(dataPre(synIdx).dataM, 2, 'omitnan');
            sdData = std(dataPre(synIdx).dataM, 0, 2, 'omitnan');
            timeVec = dataPre(synIdx).t;
            plotShadedError(axPre, timeVec, meanData, sdData, colPreLine, colPreFill);
        end
        
        % Formatting
        xline(axPre, 0, 'k--');
        xlim(axPre, [-15, 15]);
        ylabel(axPre, 'Amplitude [A.U.]', 'FontWeight', 'bold', 'FontSize', 9);
        
        % Text Label inside plot (Top Left)
        text(axPre, 0.05, 0.9, synLabel, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 9);

        if row < 4, set(axPre, 'XTickLabel', []); end
        if row == 1
             text(axPre, 0.5, 1.15, figTitle, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end

        % =============================================================
        % 2. Post-Surgery Plot (Individual Days) - Column 2 (or 4)
        % =============================================================
        axPost = nexttile((row - 1) * 4 + colOffset + 2);
        hold(axPost, 'on');

        if ~isempty(dataPost(synIdx).dataM)
            numDays = size(dataPost(synIdx).dataM, 2);
            for k = 1:numDays
                dayVal = dataPost(synIdx).daysDiff(k);
                colorIdx = round(((dayVal - dayRange(1)) / (dayRange(2) - dayRange(1))) * (nColors - 1)) + 1;
                colorIdx = max(1, min(nColors, colorIdx)); 
                thisColor = cMap(colorIdx, :);
                plot(axPost, dataPost(synIdx).t, dataPost(synIdx).dataM(:, k), '-', 'Color', thisColor, 'LineWidth', 1);
            end
        end
        
        % Formatting
        xline(axPost, 0, 'k--');
        xlim(axPost, [-15, 15]);
        set(axPost, 'YTickLabel', []);
        if row < 4, set(axPost, 'XTickLabel', []); end
        
        % Axis Labels (Bottom Row Only)
        if row == 4
            xlabel(axPre, 'task range[%]', 'FontSize', 9);
            xlabel(axPost, 'task range[%]', 'FontSize', 9);
        end

        % Colorbar (Last row, last column of the monkey block)
        if row == 4
            colormap(axPost, cMap); % Apply specific map to this axis
            c = colorbar(axPost);
            c.Label.String = cbarLabel;
            cbarLimits = caxis(axPost); 
            caxis(axPost, dayRange); 
        end
    end
end

% Save figure
exportgraphics(f, fullfile(outFigDir, 'FigureS7.png'), 'Resolution', 300);
fprintf('Figure S7 saved to: %s\n', fullfile(outFigDir, 'FigureS7.png'));


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
    matPath = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    
    if ~exist(matPath, 'file')
         warning('Data file not found: %s', matPath);
         dataAll = struct('t', [], 'dataM', [], 'daysDiff', []);
         dataAll = repmat(dataAll, length(focusList), 1);
         return;
    end

    dataFile = load(matPath, 'synergyData', 'percentTime', 'expDates', 'nameList');
    dataCell = dataFile.synergyData; % Synergy data is usually {synergy x day}
    
    nList = length(focusList);
    dataAll = struct([]);

    for i = 1:nList
        % Find index: In synergy files, nameList is usually {'A','B','C','D'}
        idx = find(strcmp(focusList{i}, dataFile.nameList));
        
        if isempty(idx)
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
        % Check data dimensions in first cell to init array
        if numDays > 0 && ~isempty(dataCell{idx,1})
             sz = size(dataCell{idx,1}, 2); % usually [1 x time] or [time x 1]
             % We need standard length matching t
             muscleDataM = zeros(length(t), numDays);
        else
             muscleDataM = [];
        end
        
        validDays = false(1, numDays);

        for dayI = 1:numDays
            dailyData = dataCell{idx, dayI};
            if ~isempty(dailyData)
                % dailyData might be [trials x time] or [1 x time]
                % We want mean across trials
                if size(dailyData, 2) == length(t)
                     % if [trials x time]
                     meanProfile = mean(dailyData, 1, 'omitnan')';
                elseif size(dailyData, 1) == length(t)
                     % if [time x trials]
                     meanProfile = mean(dailyData, 2, 'omitnan');
                else
                     % Try transposing
                     if size(dailyData, 1) == length(t)
                        meanProfile = mean(dailyData, 2, 'omitnan');
                     else
                        meanProfile = mean(dailyData, 1, 'omitnan')'; 
                     end
                end
                
                % Ensure length matches t
                if length(meanProfile) == length(t)
                    muscleDataM(:, dayI) = meanProfile;
                    validDays(dayI) = true;
                end
            end
        end
        
        dataAll(i).t = t;
        if ~isempty(muscleDataM)
            dataAll(i).dataM = muscleDataM(:, validDays);
            dataAll(i).daysDiff = daysDiff(validDays);
        else
            dataAll(i).dataM = [];
            dataAll(i).daysDiff = [];
        end
    end
end