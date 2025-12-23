% =========================================================================
% SCRIPT: FigureS3.m
%
% PURPOSE: 
%   Generates Supplementary Figure S3: EMG profiles for landmark days.
%
% LAYOUT OPTIMIZATION (8 Rows x 6 Cols):
%   - Monkey A Pre:  2 rows x 5 cols (10 muscles)
%   - Monkey A Post: 2 rows x 5 cols (10 muscles)
%   - Monkey B Pre:  2 rows x 6 cols (12 muscles)
%   - Monkey B Post: 2 rows (5 then 4) (9 muscles) -> Changed from 3x3 to maximize size
% 
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
baseDir = fileparts(scriptPath); 
matDir = fullfile(baseDir, 'Data', 'emg', 'emg_mat');
outFigDir = fullfile(scriptPath, 'outputFigures_S3');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

if ~exist(matDir, 'dir'), error('Data directory not found: %s', matDir); end

% --- Muscles & Landmark Days ---
% Monkey A
musclesA = {'EDCdist', 'ED23', 'ECU', 'ECR', 'FDSdist', 'FDP', 'FCU', 'PL', 'FCR', 'BRD'};
dispNamesA = strrep(musclesA, 'dist', ''); 
lmdDaysA = [31, 64, 69, 79, 99];
colorsA = [1 0 0; 1 0.5 0; 0 1 1; 0 0 1; 0 0 0]; 

% Monkey B
musclesB_Pre = {'EDC', 'ED23', 'ED45', 'ECU', 'ECR', 'Deltoid', 'FDS', 'FDP', 'FCU', 'PL', 'FCR', 'BRD'};
musclesB_Post = {'EDC', 'ED23', 'ED45', 'ECU', 'ECR', 'Deltoid', 'FDP', 'PL', 'BRD'}; 
lmdDaysB = [22, 36, 44, 48, 64];
colorsB = [1 0 0; 1 0.5 0; 0 1 1; 0 0 1; 0 0 0]; 

% Dates
ttDateA = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
ttDateB = datetime('2020-01-21','InputFormat','yyyy-MM-dd');

% Figure Setup
% 8 Rows (Total) x 6 Cols (Max Width)
% increased dimensions to ensure panels are physically large
figWidth_cm = 28; figHeight_cm = 37; 
f = figure('Name', 'Figure S3', 'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, figWidth_cm, figHeight_cm]);

% Use 'none' spacing to maximize internal area
T = tiledlayout(8, 6, 'TileSpacing', 'none', 'Padding', 'compact');


%% 2. MAIN PLOTTING LOOPS
% =========================================================================

% --- BLOCK 1: Monkey A Pre (Rows 1-2) ---
% 10 Muscles -> 2 Rows of 5. Occupies Cols 1-5.
fprintf('Processing Monkey A Pre...\n');
dataPreA = loadMatData('EMG', 'Yachimun', 'PreAll', matDir, musclesA, ttDateA);
plotBlock(T, dataPreA, dispNamesA, 1, 1, 2, 5, 6, 'Pre', [], [], 'Monkey A');

% --- BLOCK 2: Monkey A Post (Rows 3-4) ---
% 10 Muscles -> 2 Rows of 5. Occupies Cols 1-5.
fprintf('Processing Monkey A Post...\n');
dataPostA = loadMatData('EMG', 'Yachimun', 'PostAll', matDir, musclesA, ttDateA);
plotBlock(T, dataPostA, dispNamesA, 3, 1, 2, 5, 6, 'Post', lmdDaysA, colorsA, '');

% --- BLOCK 3: Monkey B Pre (Rows 5-6) ---
% 12 Muscles -> 2 Rows of 6. Occupies Cols 1-6.
fprintf('Processing Monkey B Pre...\n');
dataPreB = loadMatData('EMG', 'Seseki', 'PreAll', matDir, musclesB_Pre, ttDateB);
plotBlock(T, dataPreB, musclesB_Pre, 5, 1, 2, 6, 6, 'Pre', [], [], 'Monkey B');

% --- BLOCK 4: Monkey B Post (Rows 7-8) ---
% 9 Muscles -> Reorganized to 2 Rows (5 then 4) to allow bigger panels
% Row 7: 5 muscles (Cols 1-5)
% Row 8: 4 muscles (Cols 1-4)
fprintf('Processing Monkey B Post...\n');
dataPostB = loadMatData('EMG', 'Seseki', 'PostAll', matDir, musclesB_Post, ttDateB);

% Custom call for split layout
blockStartRow = 7;
blockMuscles = musclesB_Post;
nMusc = length(blockMuscles);
totalCols = 6; 
% Split 9 muscles into row of 5 and row of 4
for i = 1:nMusc
    if i <= 5
        r = blockStartRow; 
        c = i;
    else
        r = blockStartRow + 1;
        c = i - 5;
    end
    
    tileIdx = (r-1)*totalCols + c;
    ax = nexttile(T, tileIdx);
    
    % Plot Single Muscle
    plotSingleMuscle(ax, dataPostB(i), blockMuscles{i}, 'Post', lmdDaysB, colorsB);
    
    % Axes Labels (Only on bottom-left-most tile of this block)
    if i == 6 % First tile of bottom row
        xlabel(ax, 'Task range [%]', 'FontSize', 9);
        ylabel(ax, 'Amplitude [µV]', 'FontWeight', 'bold', 'FontSize', 9);
    else
        set(ax, 'XTickLabel', []);
        set(ax, 'YTickLabel', []);
    end
    
    % Legend on last tile
    if i == nMusc
         addLegend(ax, dataPostB(i), lmdDaysB, colorsB);
    end
end


% --- Finalize and Save ---
exportgraphics(f, fullfile(outFigDir, 'FigureS3.png'), 'Resolution', 300);
fprintf('Figure S3 saved to: %s\n', fullfile(outFigDir, 'FigureS3.png'));


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotBlock(T, dataStruct, muscleNames, startRow, startCol, nRows, nCols, totalCols, condition, lmdDays, colors, blockTitle)
    nMuscles = length(muscleNames);
    
    for i = 1:nMuscles
        % Calculate tile index
        row = startRow + ceil(i / nCols) - 1;
        col = startCol + mod(i - 1, nCols);
        tileIdx = (row - 1) * totalCols + col;
        ax = nexttile(T, tileIdx);

        plotSingleMuscle(ax, dataStruct(i), muscleNames{i}, condition, lmdDays, colors);
        
        % Block Title
        if i == 1 && ~isempty(blockTitle)
             text(ax, -0.1, 1.15, blockTitle, 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        end
        
        % Axis Labels (Bottom-left logic relative to this block)
        isBottomLeft = (ceil(i / nCols) == nRows) && (mod(i - 1, nCols) == 0);
        if isBottomLeft
            xlabel(ax, 'Task range [%]', 'FontSize', 9);
            ylabel(ax, 'Amplitude [µV]', 'FontWeight', 'bold', 'FontSize', 9);
        else
            set(ax, 'XTickLabel', []);
            set(ax, 'YTickLabel', []);
        end
        
        % Legend
        if strcmp(condition, 'Post') && i == nMuscles
             addLegend(ax, dataStruct(i), lmdDays, colors);
        end
    end
end

function plotSingleMuscle(ax, dStruct, name, condition, lmdDays, colors)
    hold(ax, 'on');
    xline(ax, 0, 'k--');

    if strcmp(condition, 'Pre')
        if ~isempty(dStruct.dataM)
            meanData = mean(dStruct.dataM, 2, 'omitnan');
            plot(ax, dStruct.t, meanData, 'k-', 'LineWidth', 1.5);
        end
    else
        if ~isempty(dStruct.dataM)
            for k = 1:length(lmdDays)
                dayIdx = find(dStruct.daysDiff == lmdDays(k));
                if ~isempty(dayIdx)
                    plot(ax, dStruct.t, dStruct.dataM(:, dayIdx), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
                end
            end
        end
    end

    title(ax, name, 'FontWeight', 'bold', 'FontSize', 10);
    xlim(ax, [-15, 15]);
    box(ax, 'off'); set(ax, 'TickDir', 'out', 'FontSize', 8);
    % axis(ax, 'square'); % Optional: Force exact square
end

function addLegend(ax, dStruct, lmdDays, colors)
    % Create dummy plots for legend
    L = []; txt = {};
    for k = 1:length(lmdDays)
        L(k) = plot(ax, nan, nan, '-', 'Color', colors(k,:), 'LineWidth', 1.5);
        txt{k} = sprintf('Day %d', lmdDays(k));
    end
    legend(ax, L, txt, 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 8);
end

% --- Data Loading Function ---
function dataAll = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)
    matPath = fullfile(matDir, ['emgData_' monkeyName '_' condition '.mat']);
    if ~exist(matPath, 'file')
         warning('Data file not found: %s', matPath);
         dataAll = struct('t', [], 'dataM', [], 'daysDiff', []); dataAll = repmat(dataAll, length(focusList), 1); return;
    end
    dataFile = load(matPath, 'emgData', 'percentTime', 'expDates', 'nameList');
    dataCell = dataFile.emgData; nList = length(focusList); dataAll = struct([]);
    for i = 1:nList
        muscleIdx = find(strcmp(focusList{i}, dataFile.nameList));
        if isempty(muscleIdx), muscleIdx = find(strcmp([focusList{i} 'dist'], dataFile.nameList)); end
        if isempty(muscleIdx)
            dataAll(i).t = dataFile.percentTime; dataAll(i).dataM = []; dataAll(i).daysDiff = []; continue;
        end
        t = dataFile.percentTime;
        expDates = datetime(cellstr(dataFile.expDates), 'InputFormat', 'yyMMdd');
        daysDiff = days(expDates - ttDate);
        numDays = length(expDates); muscleDataM = zeros(length(t), numDays); validDays = false(1, numDays);
        for dayI = 1:numDays
            dailyTrials = dataCell{muscleIdx, dayI};
            if ~isempty(dailyTrials)
                muscleDataM(:, dayI) = mean(dailyTrials, 1, 'omitnan')';
                validDays(dayI) = true;
            end
        end
        dataAll(i).t = t; dataAll(i).dataM = muscleDataM(:, validDays); dataAll(i).daysDiff = daysDiff(validDays);
    end
end