% =========================================================================
% SCRIPT: Figure13.m (Tenodesis Coupling)
%
% PURPOSE: 
%   Generates Figure 13: Tenodesis Coupling Analysis (Wrist vs MCP Angles).
%   Analyze kinematic coupling refinement in Monkey B over 11 sessions.
%
% LAYOUT (AllDays Option):
%   - Facet Grid: 11 individual day plots arranged in rows.
%   - Combined Plot: One large plot aggregating all days (Bottom Right).
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

% Input Directory
dataDir = fullfile(baseDir, 'Data', 'kinematics');
fileMCP   = 'MP-joint_angles.xlsx';
fileWrist = 'MP-wrist_angles.xlsx';

% Output Directory
outFigDir =  fullfile(baseDir, 'outputFigures_Fig13');
if ~exist(outFigDir, 'dir'), mkdir(outFigDir); end

% Verify Data
if ~exist(dataDir, 'dir')
    error('Data folder not found at: %s\n(Did you download the ''Data'' folder from GitHub?)', dataDir); 
end

% Analysis Settings
plotOption = 'AllDays'; % 'LandmarkDays' or 'AllDays'
saveFig    = true;

% Day Mappings
allDayLabels = [-4, 22, 36, 44, 49, 50, 55, 56, 58, 62, 63];
colors = parula(length(allDayLabels)); 

%% 2. DATA LOADING
% -------------------------------------------------------------------------
fprintf('Loading Kinematic Data...\n');
try
    mcpData   = readmatrix(fullfile(dataDir, fileMCP));
    wristData = readmatrix(fullfile(dataDir, fileWrist));
catch ME
    error('Data load failed: %s\nCheck paths in %s', ME.message, dataDir);
end

% Check Dimensions (Expect 20 trials x 11 days)
if ~isequal(size(mcpData), [20, 11])
    warning('Unexpected data dimensions. Expected 20x11, got %dx%d.', size(mcpData));
end

%% 3. SETUP PLOT LAYOUT
% -------------------------------------------------------------------------
if strcmp(plotOption, 'AllDays')
    dayIndices = 1:11;
    nPlots = 11;
    % Layout: 3 Rows x 5 Cols. 
    % Rows 1-2 fill normally. Row 3 fills 3 slots. Combined plot takes 2x2 at bottom right.
    nRows = 3; nCols = 5;
    combIdx = [9, 10, 14, 15]; % Combined plot spans these subplot indices
    
else % LandmarkDays
    dayIndices = [1, 2, 3, 4, 5, 11];
    nPlots = 6;
    nRows = 2; nCols = 3;
    combIdx = [];
end

% Determine Global Axis Limits (for consistent scaling)
selWrist = wristData(:, dayIndices);
selMCP   = mcpData(:, dayIndices);

wRange = range(selWrist(:)); mRange = range(selMCP(:));
limWrist = [min(selWrist(:)) - 0.1*wRange, max(selWrist(:)) + 0.1*wRange];
limMCP   = [min(selMCP(:))   - 0.1*mRange, max(selMCP(:))   + 0.1*mRange];

%% 4. GENERATE FIGURE
% -------------------------------------------------------------------------
fig = figure('Name', ['Tenodesis Coupling: ' plotOption], 'Color', 'w', 'WindowState', 'maximized');
sgtitle(['Refinement of Tenodesis Coupling (' plotOption ')'], 'FontSize', 14, 'FontWeight', 'bold');

% --- A. Plot Individual Days (Facets) ---
for i = 1:nPlots
    colIdx = dayIndices(i);
    
    % Determine Subplot Position (Custom 5+3+3 logic for AllDays)
    if strcmp(plotOption, 'AllDays')
        if i <= 5,     spIdx = i;      % Row 1
        elseif i <= 8, spIdx = i;      % Row 2 (Left)
        else,          spIdx = i + 2;  % Row 3 (Left, skip combined slots)
        end
    else
        spIdx = i;
    end
    
    ax = subplot(nRows, nCols, spIdx);
    
    % Call Helper to Plot Scatter & Fit
    plotDayFit(ax, wristData(:, colIdx), mcpData(:, colIdx), ...
               colors(i,:), allDayLabels(i), limWrist, limMCP);
end

% --- B. Plot Combined Summary (Bottom Right) ---
if strcmp(plotOption, 'AllDays') && ~isempty(combIdx)
    axComb = subplot(nRows, nCols, combIdx);
    hold(axComb, 'on');
    
    for i = 1:nPlots
        scatter(axComb, wristData(:, dayIndices(i)), mcpData(:, dayIndices(i)), ...
                30, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    title(axComb, 'All Days Combined', 'FontSize', 11);
    xlabel(axComb, 'Wrist Angle (deg)', 'FontSize', 9);
    ylabel(axComb, 'MCP Angle (deg)', 'FontSize', 9);
    xlim(axComb, limWrist); ylim(axComb, limMCP);
    grid(axComb, 'off'); set(axComb, 'TickDir', 'out', 'FontSize', 9);
    
    % Colorbar Setup
    colormap(axComb, colors);
    cb = colorbar(axComb);
    cb.Label.String = 'Days Relative to Surgery';
    cb.Ticks = linspace(1/(2*nPlots), 1 - 1/(2*nPlots), nPlots);
    cb.TickLabels = allDayLabels(dayIndices);
    cb.FontSize = 9;
end

%% 5. SAVE FIGURE (Updated for SVG Compatibility)
% -------------------------------------------------------------------------
if saveFig
    fNameBase = fullfile(outFigDir, ['Tenodesis_' plotOption]);
    
    % A4 Dimensions (cm)
    A4_W = 29.7; A4_H = 21.0; Margin = 1.5;
    
    set(fig, 'PaperUnits', 'centimeters', 'PaperOrientation', 'landscape');
    set(fig, 'PaperSize', [A4_W, A4_H]);
    
    % Calc optimal position keeping screen aspect ratio
    scrPos = get(fig, 'Position'); % Pixels
    aspRatio = scrPos(3) / scrPos(4);
    
    availW = A4_W - 2*Margin;
    availH = A4_H - 2*Margin;
    
    if (availW / availH) > aspRatio
        figH = availH; figW = figH * aspRatio;
    else
        figW = availW; figH = figW / aspRatio;
    end
    
    leftPos = (A4_W - figW) / 2;
    botPos  = (A4_H - figH) / 2;
    
    set(fig, 'PaperPosition', [leftPos, botPos, figW, figH]);
    
    % Save .fig (Editable)
    saveas(fig, [fNameBase '.fig']);
    
    % Save .png (High Res Raster)
    exportgraphics(fig, [fNameBase '.png'], 'Resolution', 300);
    
    % Save Vector (Robust Method)
    try
        print(fig, [fNameBase '.svg'], '-dsvg', '-painters');
        fprintf('Saved SVG: %s\n', [fNameBase '.svg']);
    catch
        exportgraphics(fig, [fNameBase '.pdf'], 'ContentType', 'vector');
        fprintf('SVG failed. Saved PDF instead: %s\n', [fNameBase '.pdf']);
    end
    
    fprintf('Figure saved to: %s\n', outFigDir);
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================
function plotDayFit(ax, xData, yData, col, dayLabel, xLim, yLim)
    % Plots Scatter points and Linear Regression line
    
    scatter(ax, xData, yData, 50, col, 'filled', 'MarkerFaceAlpha', 0.7);
    hold(ax, 'on');
    
    % Linear Regression (Fit)
    [p, ~] = polyfit(xData, yData, 1);
    xFit = linspace(min(xData), max(xData), 10);
    yFit = polyval(p, xFit);
    
    plot(ax, xFit, yFit, 'k-', 'LineWidth', 1.5);
    
    % Calc R-squared
    corrMat = corrcoef(xData, yData);
    if numel(corrMat) > 1 && ~any(isnan(corrMat(:)))
        r2 = corrMat(1,2)^2;
    else
        r2 = NaN;
    end
    
    % Formatting
    title(ax, sprintf('Day %d', dayLabel), 'FontSize', 11);
    xlabel(ax, 'Wrist (deg)', 'FontSize', 9);
    ylabel(ax, 'MCP (deg)', 'FontSize', 9);
    xlim(ax, xLim); ylim(ax, yLim);
    
    text(ax, 0.05, 0.9, sprintf('R^2 = %.2f', r2), ...
        'Units', 'normalized', 'FontSize', 9, 'FontWeight', 'bold');
    
    grid(ax, 'off'); 
    set(ax, 'FontSize', 9, 'TickDir', 'out');
    hold(ax, 'off');
end