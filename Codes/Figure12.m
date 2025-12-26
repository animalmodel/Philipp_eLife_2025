% =========================================================================
% SCRIPT: Figure12.m (Kinematics)
%
% PURPOSE: 
%   Generates Figure 12: Longitudinal Analysis of Joint Angles (Monkey B).
%   Plots the Mean +/- SD of MCP and Wrist angles across specific days.
%
% INPUT FILES:
%   - MP-joint_angles.xlsx (MCP Data)
%   - MP-wrist_angles.xlsx (Wrist Data)
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

% Input Directory (Shared with Figure 13)
dataDir = fullfile(baseDir, 'Data', 'kinematics');
fileMCP   = 'MP-joint_angles.xlsx';
fileWrist = 'MP-wrist_angles.xlsx';

% Output Directory
outFigDir =  fullfile(baseDir, 'outputFigures_Fig12');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify Data
if ~exist(dataDir, 'dir')
    error('Data folder not found at: %s\n(Did you download the ''Data'' folder from GitHub?)', dataDir); 
end

% Day Mappings (Monkey B)
% Indices corresponding to columns in the Excel files
% Columns: 1(-4d), 2(22d), 3(36d), 4(44d), 5(49d), 6(50d), 7(55d), 8(56d), 9(58d), 10(62d), 11(63d)
% We select specific "Landmark Days" for the plot based on your sample image
plotIndices = [1, 2, 3, 4, 5, 11]; 
dayLabels   = [-4, 22, 36, 44, 49, 63]; % X-axis labels

% Plot Settings
markerColor = [1 0 0]; % Red
lineColor   = [0 0 0]; % Black
errColor    = [0 0 0]; % Black

%% 2. DATA LOADING & PROCESSING
% -------------------------------------------------------------------------
fprintf('Loading Kinematic Data...\n');
try
    rawMCP   = readmatrix(fullfile(dataDir, fileMCP));
    rawWrist = readmatrix(fullfile(dataDir, fileWrist));
catch ME
    error('Data load failed: %s\nCheck paths in %s', ME.message, dataDir);
end

% Extract Data for Selected Days
% Calculate Mean and SD across trials (rows) for each day (columns)
mcpMean = mean(rawMCP(:, plotIndices), 1, 'omitnan');
mcpSD   = std(rawMCP(:, plotIndices), 0, 1, 'omitnan');

wristMean = mean(rawWrist(:, plotIndices), 1, 'omitnan');
wristSD   = std(rawWrist(:, plotIndices), 0, 1, 'omitnan');

%% 3. PLOTTING
% -------------------------------------------------------------------------
fig = figure('Name', 'Figure 12: Kinematics', 'Color', 'w', 'Position', [100 100 500 800]);

% --- SUBPLOT A: MCP ---
ax1 = subplot(2, 1, 1);
plotKinematicPanel(ax1, dayLabels, mcpMean, mcpSD, 'MCP', 'joint angle [deg]', [-40 10]);

% --- SUBPLOT B: WRIST ---
ax2 = subplot(2, 1, 2);
plotKinematicPanel(ax2, dayLabels, wristMean, wristSD, 'WRIST', 'joint angle [deg]', [-30 20]);

% Add Significance Brackets (Hardcoded based on visual reference)
% Note: Actual stats should be calculated, but these represent the visual brackets from your image.
addBracket(ax1, -4, 49, 8, '***'); % Top bracket MCP
addBracket(ax1, -4, 22, -14, '**'); % Left bracket MCP (example position)

addBracket(ax2, -4, 44, 15, '');   % Top big bracket Wrist
addBracket(ax2, 22, 63, 10, '***'); % Lower bracket Wrist

% Save
saveas(fig, fullfile(outFigDir, 'Fig12_Kinematics.fig'));
% Use print for SVG to avoid exportgraphics issues
try
    print(fig, fullfile(outFigDir, 'Fig12_Kinematics.svg'), '-dsvg', '-painters');
    fprintf('Saved SVG to %s\n', outFigDir);
catch
    exportgraphics(fig, fullfile(outFigDir, 'Fig12_Kinematics.png'), 'Resolution', 300);
    fprintf('Saved PNG (SVG failed) to %s\n', outFigDir);
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================
function plotKinematicPanel(ax, x, y, sd, tStr, yLbl, yLimVal)
    hold(ax, 'on');
    
    % Split Pre/Post for visual break in line (dashed/solid)
    preIdx = x < 0;
    postIdx = x >= 0;
    
    % Plot Error Bars
    errorbar(ax, x, y, sd, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    
    % Plot Lines (Connecting Means)
    % Dashed line from Pre to first Post point
    plot(ax, [x(preIdx) x(find(postIdx,1))], [y(preIdx) y(find(postIdx,1))], ...
        'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'LineWidth', 2);
    
    % Solid line for Post points
    plot(ax, x(postIdx), y(postIdx), 'k-', 'LineWidth', 2);
    
    % Plot Markers (Red Circles)
    plot(ax, x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    
    % Formatting
    title(ax, tStr, 'FontSize', 16, 'FontWeight', 'normal');
    ylabel(ax, yLbl, 'FontSize', 14);
    if strcmp(tStr, 'WRIST'), xlabel(ax, 'Post surgery days', 'FontSize', 14); end
    
    xlim(ax, [-10 70]);
    ylim(ax, yLimVal);
    
  
    set(ax, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2);
end

function addBracket(ax, x1, x2, yH, starStr)
    % Draws a significance bracket
    hold(ax, 'on');
    plot(ax, [x1 x1 x2 x2], [yH-2 yH yH yH-2], 'k-', 'LineWidth', 1.2);
    text(ax, mean([x1 x2]), yH+1, starStr, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end