% =========================================================================
% SCRIPT: FigureS6.m
%
% PURPOSE: 
%   Generates Supplementary Figure S6: Synergy weights of all sessions.
%   - Displays weights for every recording day.
%   - Highlights main contributing muscles with red background boxes.
%
% UPDATE:
%   - Manually excluded red box for ECR in Monkey B, Synergy A (Pre).
%
% LAYOUT (8 Rows x 2 Cols):
%   - Rows 1-4: Monkey A (Synergies A-D)
%   - Rows 5-8: Monkey B (Synergies A-D)
%   - Left Col: Pre-surgery
%   - Right Col: Post-surgery
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
baseDir = fileparts(scriptPath); 
matDir = fullfile(baseDir, 'Data', 'synergy'); 
outFigDir = fullfile(scriptPath, 'outputFigures_S6');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

if ~exist(matDir, 'dir'), error('Data directory not found: %s', matDir); end

% --- Settings ---
monkeyNameList = {'Yachimun', 'Seseki'};
synergyList = {'A', 'B', 'C', 'D'};

% Muscle Lists
musclesA = {'BRD', 'ECR', 'ECU', 'ED23', 'EDC', 'FCR', 'FCU', 'FDP', 'FDS', 'PL'};
musclesB = {'BRD', 'Deltoid', 'ECR', 'ECU', 'ED23', 'ED45', 'EDC', 'FDP', 'PL'};

% Dates
ttDateA = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
ttDateB = datetime('2020-01-21','InputFormat','yyyy-MM-dd');

% Figure Setup
figWidth_cm = 24; figHeight_cm = 32; 
f = figure('Name', 'Figure S6: Synergy Weights', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, figWidth_cm, figHeight_cm]);
T = tiledlayout(8, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


%% 2. MAIN PLOTTING LOOPS
% =========================================================================

% --- Monkey A ---
fprintf('Processing Monkey A...\n');
dataPreA  = loadSynergyW('Yachimun', 'PreAll', matDir, synergyList, ttDateA);
dataPostA = loadSynergyW('Yachimun', 'PostAll', matDir, synergyList, ttDateA);

for i = 1:4
    % Pre (Left Column)
    axPre = nexttile(T, (i-1)*2 + 1);
    plotWeightBars(axPre, dataPreA(i).wAll, musclesA, ['Synergy ' synergyList{i}], i==1, i==4, {});
    if i==1, text(axPre, -0.15, 1.1, 'Monkey A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold'); end
    if i==1, title(axPre, 'Pre-surgery', 'FontSize', 11); end

    % Post (Right Column)
    axPost = nexttile(T, (i-1)*2 + 2);
    plotWeightBars(axPost, dataPostA(i).wAll, musclesA, '', false, i==4, {});
    if i==1, title(axPost, 'Post-surgery', 'FontSize', 11); end
end

% --- Monkey B ---
fprintf('Processing Monkey B...\n');
dataPreB  = loadSynergyW('Seseki', 'PreAll', matDir, synergyList, ttDateB);
dataPostB = loadSynergyW('Seseki', 'PostAll', matDir, synergyList, ttDateB);

startRow = 5;
for i = 1:4
    % Pre (Left Column)
    axPre = nexttile(T, (startRow+i-2)*2 + 1);
    
    % Override: Exclude ECR from red box for Synergy A (i=1) only
    if i == 1
        excludeList = {'ECR'};
    else
        excludeList = {};
    end
    
    plotWeightBars(axPre, dataPreB(i).wAll, musclesB, ['Synergy ' synergyList{i}], i==1, i==4, excludeList);
    
    if i==1, text(axPre, -0.15, 1.1, 'Monkey B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold'); end
    
    % Post (Right Column)
    axPost = nexttile(T, (startRow+i-2)*2 + 2);
    plotWeightBars(axPost, dataPostB(i).wAll, musclesB, '', false, i==4, {});
end

% --- Finalize and Save ---
exportgraphics(f, fullfile(outFigDir, 'FigureS6.png'), 'Resolution', 300);
fprintf('Figure S6 saved to: %s\n', fullfile(outFigDir, 'FigureS6.png'));


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotWeightBars(ax, wData, muscleLabels, yLabelStr, showYLabel, showXLabel, excludeMuscles)
    hold(ax, 'on');
    
    % 1. Extract Data Matrix
    if iscell(wData)
        nDays = length(wData);
        if nDays > 0
            nMuscles = length(wData{1});
            W = zeros(nMuscles, nDays);
            for d = 1:nDays
                if ~isempty(wData{d})
                    W(:,d) = wData{d};
                end
            end
        else
            W = [];
        end
    else
        W = wData;
    end

    if isempty(W)
        axis(ax, 'off'); return;
    end
    
    nMuscles = size(W, 1);
    
    % 2. Draw Red Background Boxes for High Weights
    meanW = mean(W, 2, 'omitnan');
    highIdx = find(meanW > 0.4);
    
    % Filter out manually excluded muscles
    if ~isempty(excludeMuscles)
        validIdx = [];
        for k = 1:length(highIdx)
            idx = highIdx(k);
            if ~ismember(muscleLabels{idx}, excludeMuscles)
                validIdx = [validIdx, idx]; %#ok<AGROW>
            end
        end
        highIdx = validIdx;
    end
    
    % Group consecutive indices
    if ~isempty(highIdx)
        groups = {}; currentGroup = highIdx(1);
        for k = 2:length(highIdx)
            if highIdx(k) == highIdx(k-1) + 1
                currentGroup = [currentGroup, highIdx(k)];
            else
                groups{end+1} = currentGroup;
                currentGroup = highIdx(k);
            end
        end
        groups{end+1} = currentGroup;
        
        for k = 1:length(groups)
            g = groups{k};
            xStart = g(1) - 0.45;
            xEnd   = g(end) + 0.45;
            fill(ax, [xStart xEnd xEnd xStart], [-0.1 -0.1 2.5 2.5], [1 0.85 0.85], 'EdgeColor', 'none');
        end
    end

    % 3. Plot Bars
    b = bar(ax, W, 'BarWidth', 1, 'EdgeColor', 'none');
    for k = 1:length(b), b(k).FaceColor = [0.2 0.3 0.7]; end
    
    % 4. Formatting
    xlim(ax, [0.5, nMuscles + 0.5]);
    ylim(ax, [0, 2]); 
    
    if showYLabel
        ylabel(ax, 'weight contribution [a.u.]', 'FontSize', 9);
        text(ax, -0.35, 0.5, yLabelStr, 'Units', 'normalized', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 10);
    else
        set(ax, 'YTickLabel', []);
    end
    
    if showXLabel
        set(ax, 'XTick', 1:nMuscles, 'XTickLabel', muscleLabels, 'XTickLabelRotation', 45, 'FontSize', 8);
    else
        set(ax, 'XTick', []);
    end
    
    box(ax, 'off');
    set(ax, 'TickDir', 'out');
end

% --- Data Loading Function ---
function dataAll = loadSynergyW(monkeyName, condition, matDir, focusList, ttDate)
    fname = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    dataAll = repmat(struct('wAll',[], 'daysDiff',[]), 1, length(focusList));
    
    if ~exist(fname, 'file'), warning('File not found: %s', fname); return; end
    
    try
        D = load(fname, 'synergyWData', 'expDates', 'nameList');
        expDates = datetime(cellstr(D.expDates), 'InputFormat','yyMMdd');
        daysDiff = days(expDates - ttDate);
        
        for i = 1:length(focusList)
            idx = find(strcmp(focusList{i}, D.nameList));
            if ~isempty(idx)
                dataAll(i).daysDiff = daysDiff;
                dataAll(i).wAll = D.synergyWData{idx}; 
            end
        end
    catch ME
        warning('Error loading %s: %s', fname, ME.message);
    end
end