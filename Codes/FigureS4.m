% =========================================================================
% SCRIPT: FigureS4.m
%
% PURPOSE: 
%   Generates Supplementary Figure S4: Cross-correlation analysis.
%
% EXCLUSIONS (Monkey B): 
%   - FDS, FCU, FCR (Missing/Excluded post-surgery data)
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
baseDir = fileparts(scriptPath); 
matDir = fullfile(baseDir, 'Data', 'emg', 'emg_mat');
outFigDir = fullfile(scriptPath, 'outputFigures_S4');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

if ~exist(matDir, 'dir'), error('Data directory not found: %s', matDir); end

% --- Muscles Lists ---
% Monkey A (10 Muscles)
musclesA = {'FDSdist', 'EDCdist', 'PL', 'FDP', 'ED23', 'BRD', 'FCU', 'ECU', 'FCR', 'ECR'};
dispNamesA = strrep(musclesA, 'dist', ''); 

% Monkey B (7 Muscles) - FDS, FCU, FCR EXCLUDED
musclesB = {'FDP', 'EDC', 'ECR', 'PL', 'ED23', 'ECU', 'BRD', 'ED45', 'Deltoid'}; 
% Note: Although ED45 and Deltoid are in this list, previous specific requests 
% for S4 only mentioned excluding FDS/FCU/FCR. If ED45/Deltoid aren't needed 
% for S4 specifically, they can be removed, but I will keep them to fill the grid 
% unless instructed otherwise.

% Dates
ttDateA = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
ttDateB = datetime('2020-01-21','InputFormat','yyyy-MM-dd');

% Figure Setup
figWidth_cm = 22; figHeight_cm = 32; 
f = figure('Name', 'Figure S4: Cross-correlation', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, figWidth_cm, figHeight_cm]);
T = tiledlayout(7, 3, 'TileSpacing', 'compact', 'Padding', 'compact');


%% 2. MAIN LOOPS
% =========================================================================

fprintf('Processing Monkey A...\n');
processMonkey(T, 'Yachimun', matDir, musclesA, dispNamesA, ttDateA, 1, 'Monkey A');

fprintf('Processing Monkey B...\n');
% Pass musclesB twice (list and display names are same)
processMonkey(T, 'Seseki', matDir, musclesB, musclesB, ttDateB, 13, 'Monkey B'); 


% --- Finalize and Save ---
exportgraphics(f, fullfile(outFigDir, 'FigureS4.png'), 'Resolution', 300);
fprintf('Figure S4 saved to: %s\n', fullfile(outFigDir, 'FigureS4.png'));


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function processMonkey(T, monkeyName, matDir, muscleList, dispNames, ttDate, startTileIdx, blockLabel)
    
    % 1. Load Data
    dataPre  = loadMatData(monkeyName, 'PreAll', matDir, muscleList, ttDate);
    dataPost = loadMatData(monkeyName, 'PostAll', matDir, muscleList, ttDate);
    
    nMuscles = length(muscleList);
    currentTileIdx = startTileIdx;

    for i = 1:nMuscles
        ax = nexttile(T, currentTileIdx);
        hold(ax, 'on');
        
        % --- Calculate Correlations ---
        if isempty(dataPre(i).dataM)
            % If pre data missing, just skip plotting content
            warning('No pre-surgery data for %s', dispNames{i});
        else
            % 1. Reference: Mean Pre-surgery Profile
            refProfile = mean(dataPre(i).dataM, 2, 'omitnan');
            refProfile = refProfile(:); 
            nRef = length(refProfile);

            % 2. Calculate Correlation for each Post-surgery Day
            postDays = dataPost(i).daysDiff;
            postData = dataPost(i).dataM;
            dailyCorrs = [];
            validDays = [];

            if ~isempty(postData)
                for d = 1:length(postDays)
                    dailyProfile = postData(:, d);
                    dailyProfile = dailyProfile(:); 
                    
                    % Interpolation Fix
                    if length(dailyProfile) ~= nRef
                        if length(dailyProfile) > 1
                            x_old = linspace(0, 1, length(dailyProfile));
                            x_new = linspace(0, 1, nRef);
                            dailyProfile = interp1(x_old, dailyProfile, x_new, 'linear')';
                            dailyProfile = dailyProfile(:);
                        else
                            continue; 
                        end
                    end
                    
                    % Calculate Pearson Correlation
                    R = corrcoef(refProfile, dailyProfile, 'Rows', 'complete');
                    if numel(R) > 1
                        dailyCorrs = [dailyCorrs, R(1,2)]; %#ok<*AGROW>
                        validDays = [validDays, postDays(d)];
                    end
                end
            end
            
            % --- Plotting ---
            plot(ax, 0, 1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
            
            if ~isempty(validDays)
                % Gap Line
                plot(ax, [0, validDays(1)], [1, dailyCorrs(1)], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
                % Data Line
                plot(ax, validDays, dailyCorrs, 'k-', 'LineWidth', 1.5);
            end
        end
        
        yline(ax, 0, 'k:');

        % Styling
        title(ax, dispNames{i}, 'FontWeight', 'bold', 'FontSize', 10);
        ylim(ax, [-1.1, 1.1]);
        % Dynamic x-lim based on data present
        if exist('validDays','var') && ~isempty(validDays) 
            maxD = max(validDays); 
        else
            maxD = 100; 
        end
        xlim(ax, [-5, maxD+5]);
        
        box(ax, 'off'); set(ax, 'TickDir', 'out', 'FontSize', 8);
        set(ax, 'YTick', [-1, 0, 1]);
        
        if i == 1
             text(ax, -0.3, 1.2, blockLabel, 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        end

        % Axis Labels
        % Logic updated for the 7-row grid (Monkey A: Rows 1-4, Monkey B: Rows 5-7)
        % Monkey A bottom: items 10 (last col), 9 (mid col), 8 (first col?? No, layout varies)
        % Simpler approach: check if tile below is empty or out of bounds
        
        % Check if we are at bottom of a column within the block
        % Monkey A (10 items): 
        % Col 1: 1, 4, 7, 10
        % Col 2: 2, 5, 8
        % Col 3: 3, 6, 9
        isBottomA = strcmp(monkeyName, 'Yachimun') && (i==10 || i==8 || i==9);
        
        % Monkey B (7 items):
        % Col 1: 1, 4, 7
        % Col 2: 2, 5
        % Col 3: 3, 6
        isBottomB = strcmp(monkeyName, 'Seseki') && (i==7 || i==5 || i==6);

        if isBottomA || isBottomB
            xlabel(ax, 'Post surgery days', 'FontSize', 9);
        end
        
        isLeftPlot = mod(currentTileIdx-1, 3) == 0;
        if isLeftPlot
            ylabel(ax, 'Cross-corr coeff', 'FontSize', 9);
        else
            set(ax, 'YTickLabel', []);
        end

        currentTileIdx = currentTileIdx + 1;
    end
end

function dataAll = loadMatData(monkeyName, condition, matDir, focusList, ttDate)
    matPath = fullfile(matDir, ['emgData_' monkeyName '_' condition '.mat']);
    if ~exist(matPath, 'file')
         warning('Data file not found: %s', matPath);
         dataAll = struct('t', [], 'dataM', [], 'daysDiff', []); 
         dataAll = repmat(dataAll, length(focusList), 1); 
         return;
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
                if size(dailyTrials, 1) == length(t)
                     if size(dailyTrials, 2) > size(dailyTrials, 1), dailyMean = mean(dailyTrials, 1, 'omitnan')'; 
                     else, dailyMean = mean(dailyTrials, 2, 'omitnan'); end
                elseif size(dailyTrials, 2) == length(t)
                    dailyMean = mean(dailyTrials, 1, 'omitnan')';
                else
                    dailyMean = mean(dailyTrials, 1, 'omitnan')'; 
                end
                dailyMean = dailyMean(:);
                if length(dailyMean) == length(t)
                    muscleDataM(:, dayI) = dailyMean;
                    validDays(dayI) = true;
                else
                    x_in = linspace(0,1,length(dailyMean));
                    x_t  = linspace(0,1,length(t));
                    muscleDataM(:, dayI) = interp1(x_in, dailyMean, x_t)';
                    validDays(dayI) = true;
                end
            end
        end
        dataAll(i).t = t; dataAll(i).dataM = muscleDataM(:, validDays); dataAll(i).daysDiff = daysDiff(validDays);
    end
end