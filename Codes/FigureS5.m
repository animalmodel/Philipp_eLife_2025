% =========================================================================
% SCRIPT: FigureS5.m
%
% PURPOSE: 
%   Generates Supplementary Figure S5: Variance Accounted For (VAF) Analysis.
%   - Searches for files containing 'vafDataList' (ignoring incomplete files).
%   - Plots 4 panels: Monkey A (Pre/Post) and Monkey B (Pre/Post).
%   - Compares Actual VAF (Blue) vs. Shuffled VAF (Red).
%
% AUTHORS: Roland Philipp (Adapted from Funato)
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
scriptPath = fileparts(mfilename('fullpath'));
if isempty(scriptPath), scriptPath = pwd; end
baseDir = fileparts(scriptPath); 

fprintf('Detected Base Directory: %s\n', baseDir);

% Output Directory
outFigDir = fullfile(baseDir, 'outputFigures_FigS5');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Search Paths (Will check content to find the correct file)
searchDirs = {
    fullfile(baseDir, 'matFiles'), ...          % Check user's folder first
    fullfile(baseDir, 'Data', 'matFiles'), ...
    fullfile(baseDir, 'Data', 'synergy')        % Check standard folder last
};

% Settings
monkeyNameList = {'Yachimun', 'Seseki'};
conditionList  = {'PreSelect', 'PostSelect'}; 

%% 2. MAIN PLOTTING
% -------------------------------------------------------------------------
fprintf('Generating Figure S5 (VAF Analysis)...\n');

figWidth_cm  = 20;  
figHeight_cm = 16; 
fig = figure('Name', 'Figure S5: VAF Analysis', 'Units', 'centimeters', ...
             'Position', [2, 2, figWidth_cm, figHeight_cm], 'Color', 'w');

t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    if mIndx == 1, lbl = 'Monkey A'; else, lbl = 'Monkey B'; end
    
    for cIndx = 1:2
        conditionName = conditionList{cIndx};
        if cIndx == 1, condLbl = 'Pre-surgery'; else, condLbl = 'Post-surgery (Landmark Days)'; end
        
        ax = nexttile;
        hold(ax, 'on');
        
        % --- SMART SEARCH ---
        % Find the file that actually contains the VAF variables
        fName = sprintf('synergyData_%s_%s.mat', monkeyName, conditionName);
        foundPath = '';
        vafDataList = [];
        vafShuffleDataList = [];
        
        for p = 1:length(searchDirs)
            testPath = fullfile(searchDirs{p}, fName);
            if exist(testPath, 'file')
                try
                    % Peek inside file variables
                    vars = whos('-file', testPath);
                    if ismember('vafDataList', {vars.name})
                        foundPath = testPath;
                        break; % Found the good file!
                    end
                catch
                    continue;
                end
            end
        end
        
        if isempty(foundPath)
            warning('Could not find a valid file with VAF data for: %s', fName);
            text(ax, 0.5, 0.5, 'Data Missing', 'HorizontalAlignment', 'center');
            continue;
        end
        
        % Load Data
        fprintf('Loading valid data from: %s\n', foundPath);
        D = load(foundPath, 'vafDataList', 'vafShuffleDataList');
        vafDataList = D.vafDataList;
        vafShuffleDataList = D.vafShuffleDataList;
        
        % --- PLOTTING ---
        synergyNum = (1:size(vafDataList, 1))';
        
        % 1. Shuffle (Chance) - Red Error Bars
        shuffle_mean = mean(vafShuffleDataList, 2);
        shuffle_std  = std(vafShuffleDataList, 0, 2);
        
        eb = errorbar(ax, synergyNum, shuffle_mean, shuffle_std, 'o-', ...
                      'LineWidth', 1, 'Color', [0.85, 0.33, 0.10], ...
                      'MarkerFaceColor', [0.85, 0.33, 0.10], 'MarkerSize', 4);
        
        % 2. Actual Data - Blue Lines per Day
        numDays = size(vafDataList, 2);
        for d = 1:numDays
            plot(ax, synergyNum, vafDataList(:, d), '-', 'Color', [0, 0.45, 0.74], 'LineWidth', 1);
            plot(ax, synergyNum, vafDataList(:, d), 'o', 'Color', [0, 0.45, 0.74], 'MarkerSize', 4);
        end
        
        % 3. Reference Line
        yline(ax, 0.8, 'k--', 'LineWidth', 1);
        
        % Formatting
        title(ax, {lbl, condLbl}, 'FontWeight', 'bold', 'FontSize', 10);
        xlabel(ax, 'Number of Synergies', 'FontSize', 9);
        ylabel(ax, 'VAF', 'FontSize', 9);
        
        xlim(ax, [0.5, size(vafDataList, 1) + 0.5]);
        ylim(ax, [0 1.05]); 
        set(ax, 'XTick', 1:size(vafDataList, 1));
        
        grid(ax, 'on'); box(ax, 'off'); set(ax, 'TickDir', 'out');
        
        % Legend (First panel only)
        if mIndx == 1 && cIndx == 1
            legend([eb, findobj(ax, 'Type','line','Color',[0, 0.45, 0.74],'Marker','o')], ...
                   {'Chance (Shuffle)', 'Actual Data'}, ...
                   'Location', 'southeast', 'Box', 'off', 'FontSize', 8);
        end
        
        hold(ax, 'off');
    end
end

% --- SAVE ---
set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto', 'PaperSize', [figWidth_cm figHeight_cm]);
try
    print(fig, fullfile(outFigDir, 'FigureS5_VAF.svg'), '-dsvg', '-painters');
    fprintf('Figure S5 saved to: %s\n', outFigDir);
catch
    exportgraphics(fig, fullfile(outFigDir, 'FigureS5_VAF.png'), 'Resolution', 300);
    fprintf('SVG failed. Saved PNG to: %s\n', outFigDir);
end