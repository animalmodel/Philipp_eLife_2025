%% Ultrasound Data Visualization
% Refactored for Roland Philipp (modified from Funato)
% Date: 2026-01-05
% Description: Loads ultrasound stimulation data for Start Frame, EDC, and FDS
% trials and generates publication-ready comparison plots.

clear; close all; clc;

% =========================================================================
% Configuration & Path Setup
% =========================================================================
% Define base directory for data (Update this path if the folder moves)
baseDataDir = fullfile(pwd, 'Ultrasound_data');

% Define plot settings for publication quality
lineColors = lines(7); % Generate distinct colors
lineWidth = 1.5;
fontSize  = 12;
fontName  = 'Arial';

%% =========================================================================
% Section 1: Load and Plot Stimulus Start Frame Data
% =========================================================================
disp('Processing Stimulus Start Frame Data...');

startFrameFile = fullfile(baseDataDir, 'stim_start_frame.mat');

if exist(startFrameFile, 'file')
    startFrameData = load(startFrameFile, 'coordinate_list');
    coordinateList = startFrameData.coordinate_list;
    
    fig1 = figure('Name', 'Index Finger Z-Displacement', 'Color', 'w');
    tiledlayout(1, 1, 'Padding', 'compact');
    nexttile;
    hold on;
    
    % Iterate through each day in the coordinate list
    legendEntries = {};
    for dayIndex = 1:length(coordinateList)
        dataDay = coordinateList{dayIndex};
        
        % Plot Index Finger Z (Column 3)
        plot(dataDay(:, 3), ...
            'LineWidth', lineWidth, ...
            'Color', lineColors(dayIndex, :));
        
        legendEntries{end+1} = sprintf('Day %d', dayIndex);
    end
    
    % Formatting
    xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel('Displacement (mm)', 'FontSize', fontSize, 'FontName', fontName);
    title('Index Finger Z-Displacement', 'FontSize', fontSize+2, 'FontName', fontName);
    legend(legendEntries, 'Location', 'best', 'Box', 'off');
    set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', fontSize, 'FontName', fontName);
    hold off;
else
    warning('File not found: %s', startFrameFile);
end

%% =========================================================================
% Section 2: Load and Plot EDC Stimulation Data
% =========================================================================
disp('Processing EDC Stimulation Data...');

edcFile = fullfile(baseDataDir, 'radial_trimmed_stim_data(-60_to_60 frame_5days).mat');

if exist(edcFile, 'file')
    edcData = load(edcFile, 'stim_average_data', 'days_list');
    edcAverageData = edcData.stim_average_data;
    
    % Funato's original code selected specific indices [1, 2, 4]. 
    % Keeping this logic but making it explicit.
    selectedDayIndices = [1, 2, 4]; 
    edcDaysList = edcData.days_list(selectedDayIndices);
    
    fig2 = figure('Name', 'EDC Stimulation Responses', 'Color', 'w');
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'EDC Stimulation', 'FontSize', fontSize+4, 'FontName', fontName);
    
    % Loop through the 3 conditions/components stored in stim_average_data
    for compIndex = 1:3
        nexttile;
        hold on;
        
        legendEntries = {};
        for k = 1:length(edcDaysList)
            % Extract data for the specific component
            currentData = edcAverageData{compIndex};
            
            % Plot the specific day (row k corresponds to days in selectedDayIndices)
            plot(currentData(k, :), ...
                'LineWidth', lineWidth, ...
                'Color', lineColors(k, :));
            
            legendEntries{end+1} = char(edcDaysList(k));
        end
        
        % Formatting per tile
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', fontSize, 'FontName', fontName);
        ylabel('Amplitude (a.u.)'); % Assuming arbitrary units, update if known (e.g., uV)
        
        % Only add legend to the first plot to avoid clutter
        if compIndex == 1
            legend(legendEntries, 'Location', 'bestoutside', 'Box', 'off');
        end
        
        % Only add x-label to the bottom plot
        if compIndex == 3
            xlabel('Elapsed Time from Stimulus (frames)', 'FontSize', fontSize, 'FontName', fontName);
        end
        hold off;
    end
else
    warning('File not found: %s', edcFile);
end

%% =========================================================================
% Section 3: Load and Plot FDS Stimulation Data
% =========================================================================
disp('Processing FDS Stimulation Data...');

fdsFile = fullfile(baseDataDir, 'ulnar_trimmed_stim_data(-60_to_60 frame_5days).mat');

if exist(fdsFile, 'file')
    fdsData = load(fdsFile, 'stim_average_data', 'days_list');
    fdsAverageData = fdsData.stim_average_data;
    fdsDaysList = fdsData.days_list;
    
    fig3 = figure('Name', 'FDS Stimulation Responses', 'Color', 'w');
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'FDS Stimulation', 'FontSize', fontSize+4, 'FontName', fontName);
    
    % Loop through the 3 conditions/components
    for compIndex = 1:3
        nexttile;
        hold on;
        
        legendEntries = {};
        for k = 1:length(fdsDaysList)
            currentData = fdsAverageData{compIndex};
            
            plot(currentData(k, :), ...
                'LineWidth', lineWidth, ...
                'Color', lineColors(k, :));
            
            legendEntries{end+1} = char(fdsDaysList(k));
        end
        
        % Formatting per tile
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', fontSize, 'FontName', fontName);
        ylabel('Amplitude (a.u.)');
        
        if compIndex == 1
            legend(legendEntries, 'Location', 'bestoutside', 'Box', 'off');
        end
        
        if compIndex == 3
            xlabel('Elapsed Time from Stimulus (frames)', 'FontSize', fontSize, 'FontName', fontName);
        end
        hold off;
    end
else
    warning('File not found: %s', fdsFile);
end

disp('Plotting complete.');