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
%   - Col 2: Monkey A Post (All days, color-coded 'hot')
%   - Col 3: Monkey B Pre (Mean +/- SD)
%   - Col 4: Monkey B Post (All days, color-coded 'summer')
%
% UPDATES:
%   - Dynamic Y-axis limits (autoscales to max data value).
%   - Square aspect ratio for better visibility.
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

% Force calculation type to Synergy
calcType = 'Synergy';

% Set data directory
matDir = fullfile(baseDir, 'Data', 'synergy');

% Set output directory
outFigDir = fullfile(baseDir, 'outputFigures_S7');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify data directory exists
if ~exist(matDir, 'dir')
    error('Data directory not found: %s', matDir);
end

monkeyNameList = {'Yachimun', 'Seseki'};
synergyList = {'A', 'B', 'C', 'D'};

% Dates
ttDateA = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
ttDateB = datetime('2020-01-21','InputFormat','yyyy-MM-dd');

% Figure Setup
figWidth_cm = 28; figHeight_cm = 24; 
f = figure('Name', 'Figure S7: Synergy Profiles', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, figWidth_cm, figHeight_cm]);
T = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 2. MAIN LOOPS
% =========================================================================

% --- Load Data ---
fprintf('Loading Data...\n');
dataPreA  = loadMatData('Yachimun', 'PreAll', matDir, synergyList, ttDateA);
dataPostA = loadMatData('Yachimun', 'PostAll', matDir, synergyList, ttDateA);
dataPreB  = loadMatData('Seseki', 'PreAll', matDir, synergyList, ttDateB);
dataPostB = loadMatData('Seseki', 'PostAll', matDir, synergyList, ttDateB);

% --- Plotting Rows (Synergies) ---
for i = 1:4
    synLabel = synergyList{i};
    
    % --- COL 1: Monkey A Pre ---
    ax1 = nexttile(T, (i-1)*4 + 1);
    plotPre(ax1, dataPreA(i), ['Synergy ' synLabel], i==1, 'Monkey A', [0.2 0.2 0.6], [0.4 0.4 0.8]);
    
    % --- COL 2: Monkey A Post ---
    ax2 = nexttile(T, (i-1)*4 + 2);
    plotPost(ax2, dataPostA(i), '', i==1, 'Post-Surgery', 'hot');
    
    % --- COL 3: Monkey B Pre ---
    ax3 = nexttile(T, (i-1)*4 + 3);
    plotPre(ax3, dataPreB(i), ['Synergy ' synLabel], i==1, 'Monkey B', [0.2 0.2 0.6], [0.4 0.4 0.8]);
    
    % --- COL 4: Monkey B Post ---
    ax4 = nexttile(T, (i-1)*4 + 4);
    plotPost(ax4, dataPostB(i), '', i==1, 'Post-Surgery', 'summer');
    
    % Match Y-Limits across the row for comparison
    linkaxes([ax1, ax2, ax3, ax4], 'y');
    
    % Axis Labels (Bottom Row Only)
    if i == 4
        xlabel(ax1, '% Cycle'); xlabel(ax2, '% Cycle');
        xlabel(ax3, '% Cycle'); xlabel(ax4, '% Cycle');
    else
        set(ax1,'XTickLabel',[]); set(ax2,'XTickLabel',[]);
        set(ax3,'XTickLabel',[]); set(ax4,'XTickLabel',[]);
    end
end

% --- Finalize and Save ---
try
    exportgraphics(f, fullfile(outFigDir, 'FigureS7.png'), 'Resolution', 300);
    fprintf('Figure S7 saved to: %s\n', fullfile(outFigDir, 'FigureS7.png'));
catch
    fprintf('Warning: Could not save PNG. Check permissions.\n');
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotPre(ax, dStruct, rowTitle, isTop, colTitle, colLine, colFill)
    hold(ax, 'on');
    xline(ax, 0, 'k:');
    
    maxY = 1.0; % Default fallback
    
    if ~isempty(dStruct.dataM)
        meanData = mean(dStruct.dataM, 2, 'omitnan');
        sdData   = std(dStruct.dataM, 0, 2, 'omitnan');
        t        = dStruct.t;
        
        % Plot Tube
        x = t(:)'; y = meanData(:)'; sd = sdData(:)';
        fill(ax, [x, fliplr(x)], [y+sd, fliplr(y-sd)], colFill, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(ax, x, y, 'Color', colLine, 'LineWidth', 2);
        
        % Calculate Max for Scaling
        maxY = max(y+sd, [], 'all', 'omitnan');
    end
    
    % Styling
    xlim(ax, [-15 15]);
    
    % Dynamic Y-Lim (with a minimum of 1.0)
    finalY = max(1.2, maxY * 1.1);
    ylim(ax, [0 finalY]); 
    
    ylabel(ax, rowTitle, 'FontWeight', 'bold', 'FontSize', 10);
    
    if isTop
        title(ax, {colTitle, 'Pre-surgery'}, 'FontWeight', 'bold', 'FontSize', 11);
    end
    
    box(ax, 'off'); set(ax, 'TickDir', 'out'); pbaspect(ax, [1 1 1]);
end

function plotPost(ax, dStruct, rowTitle, isTop, colTitle, cmapName)
    hold(ax, 'on');
    xline(ax, 0, 'k:');
    
    maxY = 1.0; 
    
    if ~isempty(dStruct.dataM)
        nDays = size(dStruct.dataM, 2);
        
        % Colormap
        colormap(ax, cmapName);
        cMap = colormap(ax);
        nColors = size(cMap, 1);
        
        if nDays > 1
             dayRange = [min(dStruct.daysDiff), max(dStruct.daysDiff)];
        else
             dayRange = [0 1];
        end
        
        for k = 1:nDays
            dayVal = dStruct.daysDiff(k);
            % Map day to color
            if diff(dayRange) == 0, cIdx = 1; 
            else
                cIdx = round( (dayVal - dayRange(1)) / diff(dayRange) * (nColors-1) ) + 1;
            end
            cIdx = max(1, min(nColors, cIdx));
            
            plot(ax, dStruct.t, dStruct.dataM(:, k), 'Color', cMap(cIdx,:), 'LineWidth', 1);
        end
        
        % Colorbar (only for top plot to save space, or per plot)
        if isTop
             cb = colorbar(ax);
             cb.Label.String = 'Days Post';
             caxis(ax, dayRange);
        end
        
        maxY = max(dStruct.dataM, [], 'all', 'omitnan');
    end
    
    % Styling
    xlim(ax, [-15 15]);
    
    finalY = max(1.2, maxY * 1.1);
    ylim(ax, [0 finalY]);
    
    set(ax, 'YTickLabel', []);
    
    if isTop
        title(ax, colTitle, 'FontWeight', 'bold', 'FontSize', 11);
    end
    
    box(ax, 'off'); set(ax, 'TickDir', 'out'); pbaspect(ax, [1 1 1]);
end


function dataAll = loadMatData(monkeyName, condition, matDir, focusList, ttDate)
    fname = fullfile(matDir, ['synergyData_' monkeyName '_' condition '.mat']);
    dataAll = repmat(struct('t',[], 'dataM',[], 'daysDiff',[]), 1, length(focusList));
    
    if ~exist(fname, 'file'), warning('File not found: %s', fname); return; end
    
    try
        D = load(fname, 'synergyData', 'percentTime', 'expDates', 'nameList');
        expDates = datetime(cellstr(D.expDates), 'InputFormat','yyMMdd');
        daysDiff = days(expDates - ttDate);
        t = D.percentTime;
        
        for i = 1:length(focusList)
            idx = find(strcmp(focusList{i}, D.nameList));
            
            if ~isempty(idx)
                numDays = length(expDates);
                muscleDataM = zeros(length(t), numDays);
                validDays = false(1, numDays);
                
                for dayI = 1:numDays
                    dailyData = D.synergyData{idx, dayI};
                    if ~isempty(dailyData)
                        % Robust Mean Calculation
                        % dailyData might be [trials x time] or [time x trials]
                        if size(dailyData, 2) == length(t)
                             meanProfile = mean(dailyData, 1, 'omitnan')';
                        elseif size(dailyData, 1) == length(t)
                             meanProfile = mean(dailyData, 2, 'omitnan');
                        else
                             % Last resort: check elements
                             if numel(dailyData) == length(t)
                                 meanProfile = dailyData(:);
                             else
                                 meanProfile = [];
                             end
                        end
                        
                        if ~isempty(meanProfile)
                            muscleDataM(:, dayI) = meanProfile;
                            validDays(dayI) = true;
                        end
                    end
                end
                
                dataAll(i).t = t;
                dataAll(i).dataM = muscleDataM(:, validDays);
                dataAll(i).daysDiff = daysDiff(validDays);
            end
        end
    catch ME
        warning('Error loading %s: %s', fname, ME.message);
    end
end