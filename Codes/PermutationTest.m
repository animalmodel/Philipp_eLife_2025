% =========================================================================
% SCRIPT: PermutationTest.m
%
% PURPOSE: 
%   Performs a permutation test (Monte Carlo) to statistically compare
%   synergy weights between two conditions (e.g., Pre vs. Post).
%   
%   - Loads synergy data (W matrices).
%   - Calculates the Cosine Similarity of the actual data.
%   - Shuffles the muscle weights randomly (10,000 times) to build a null distribution.
%   - Calculates a p-value based on where the actual similarity falls in the distribution.
%
% AUTHORS: Roland Philipp
% =========================================================================

clear; clc; close all;

%% 1. CONFIGURATION & PATHS
% -------------------------------------------------------------------------
% --- DYNAMIC PATH SETUP ---
scriptPath = fileparts(mfilename('fullpath'));
if isempty(scriptPath), scriptPath = pwd; end % Fallback for running sections
baseDir = fileparts(scriptPath); 

fprintf('Detected Base Directory: %s\n', baseDir);

% Input Directory (Data/synergy)
matDir = fullfile(baseDir, 'Data', 'synergy'); 

% Verify Data Path Exists
if ~exist(matDir, 'dir')
    error('Data folder not found at: %s\nPlease check directory structure.', matDir);
end

% Settings
monkeyNameList = {'Yachimun', 'Seseki'};
synergyLabels  = {'A', 'B', 'C', 'D'};
nPermutations  = 10000; % Number of shuffles

%% 2. MAIN LOOP
% -------------------------------------------------------------------------
fprintf('Running Permutation Tests (%d iterations)...\n', nPermutations);

for mIndx = 1:2
    monkeyName = monkeyNameList{mIndx};
    fprintf('\n--- %s ---\n', monkeyName);
    
    % Load Data
    dataPre  = loadSynergyW(monkeyName, 'PreAll', matDir, synergyLabels);
    dataPost = loadSynergyW(monkeyName, 'PostAll', matDir, synergyLabels);
    
    for s = 1:4
        % --- ROBUST DATA EXTRACTION (PRE) ---
        rawPre = dataPre(s).wAll;
        if iscell(rawPre)
             % Cell array case: Remove empty cells and concat
             rawPre = rawPre(~cellfun('isempty', rawPre));
             if isempty(rawPre)
                 wPreAll = [];
             else
                 wPreAll = [rawPre{:}];
             end
        else
             % Matrix case: Use directly
             wPreAll = rawPre;
        end
        
        if isempty(wPreAll)
            fprintf('Synergy %s: No Pre-data found.\n', synergyLabels{s});
            continue;
        end
        
        wRef = mean(wPreAll, 2);
        
        % --- ROBUST DATA EXTRACTION (POST) ---
        rawPost = dataPost(s).wAll;
        wPostAll = [];
        
        if iscell(rawPost)
            % Cell array case
            rawPost = rawPost(~cellfun('isempty', rawPost));
            nPost = length(rawPost);
            if nPost > 0
                if nPost > 5
                    wPostAll = [rawPost{nPost-4:nPost}]; % Last 5 sessions
                else
                    wPostAll = [rawPost{:}];
                end
            end
        else
            % Matrix case
            nPost = size(rawPost, 2);
            if nPost > 0
                if nPost > 5
                    wPostAll = rawPost(:, nPost-4:nPost); % Last 5 cols
                else
                    wPostAll = rawPost;
                end
            end
        end
        
        if isempty(wPostAll)
            fprintf('Synergy %s: No Post-data found.\n', synergyLabels{s});
            continue;
        end
        
        wTest = mean(wPostAll, 2);
        
        % 1. Calculate Actual Similarity
        actualSim = getCosineSim(wRef, wTest);
        
        % 2. Permutation Test
        nullDist = zeros(1, nPermutations);
        nMuscles = length(wRef);
        
        for k = 1:nPermutations
            % Shuffle muscle weights in wTest
            wShuffled = wTest(randperm(nMuscles));
            nullDist(k) = getCosineSim(wRef, wShuffled);
        end
        
        % 3. Calculate P-Value
        % Null Hypothesis: The similarity is random.
        % We check how many random shuffles produced a similarity HIGHER than actual.
        % (One-tailed test: Is the structure significantly conserved?)
        p_random = (sum(nullDist >= actualSim) + 1) / (nPermutations + 1);
        
        fprintf('Synergy %s: Similarity = %.3f, p(random) = %.4f\n', ...
                synergyLabels{s}, actualSim, p_random);
    end
end

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function sim = getCosineSim(v1, v2)
    if norm(v1)==0 || norm(v2)==0, sim=0; return; end
    sim = dot(v1, v2) / (norm(v1) * norm(v2));
end

function dataAll = loadSynergyW(monkeyName, condition, matDir, focusList)
    fName = ['synergyData_' monkeyName '_' condition '.mat'];
    fullPath = fullfile(matDir, fName);
    
    dataAll = repmat(struct('wAll',[]), 1, length(focusList));
    
    if ~exist(fullPath, 'file')
        warning('File not found: %s', fullPath);
        return;
    end
    
    try
        D = load(fullPath, 'synergyWData', 'nameList');
        for i = 1:length(focusList)
            idx = find(strcmp(focusList{i}, D.nameList));
            if ~isempty(idx)
                dataAll(i).wAll = D.synergyWData{idx};
            end
        end
    catch
    end
end