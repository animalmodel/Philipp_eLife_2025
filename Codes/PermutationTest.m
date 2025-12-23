clear all;
close all

baseDir = 'C:\Users\mypre\Documents\Manuscripts\Revision\post acceptance revision\Philipp_eLife_2025';


matDir    =  fullfile(baseDir, 'Data', 'synergy');
outFigDir = 'outputFigures';
calcType = 'Synergy';

monkeyNameList = {'Yachimun', 'Seseki'};
for mIndx = 1:2
% for mIndx = 1
    monkeyName = monkeyNameList{mIndx};
    if(mIndx==1)
        ttDate   = datetime('2017-05-30','InputFormat','yyyy-MM-dd');
        focusList = {'A', 'B', 'C', 'D'};
        preDayNum = 1;
    else
        ttDate   = datetime('2020-01-21','InputFormat','yyyy-MM-dd');
        focusList = {'A', 'B', 'C', 'D'};
        preDayNum = 1;
    end
    % ===== Load Mat Data  (Selected Days for Plot)=====
    dataSelectPre  = loadMatData(calcType, monkeyName, 'PreSelect',  matDir, focusList, ttDate);
    dataSelectPost = loadMatData(calcType, monkeyName, 'PostSelect',  matDir, focusList, ttDate);
    %==%

    % ===== Synergy Plot & calc Correlation=====
    figure;
    tiledlayout(length(focusList),2,'TileSpacing','compact','Padding','compact');
    for i = 1:length(focusList)
        %- Pre -%
        synergyData = dataSelectPre(i).dataAll;
        synergyDataTrialPreMeanList = zeros(length(dataSelectPre(i).t), length(synergyData));
        dayI = 1;
        synergyDataPre = synergyData{dayI};
        nexttile;
        synergyPreTrialMean = mean(synergyDataPre,2);
        synergyDataTrialPreMeanList(:,dayI) = synergyPreTrialMean;
        synergyDataTrialSD   = std(synergyDataPre, 0, 2);
        hold on;
        plot(dataSelectPre(i).t, synergyPreTrialMean, 'b-');
        plot(dataSelectPre(i).t, synergyPreTrialMean+synergyDataTrialSD, 'c-');
        plot(dataSelectPre(i).t, synergyPreTrialMean-synergyDataTrialSD, 'c-');
        hold off;
        xlim([-15 15]);
        if(dayI==1), ylabel(focusList{i}); end
        if(i==1), title('Pre'); end
        if(i==length(focusList)), xlabel('time %'); end
        %
        %- Post -%
        synergyData = dataSelectPost(i).dataAll;
        synergyDataTrialPostMeanList = zeros(length(dataSelectPost(i).t), length(synergyData));
        % dayI = 1;
        dayI = length(synergyData);
        synergyDataPost = synergyData{dayI};
        nexttile;
        synergyPostTrialMean = mean(synergyDataPost,2);
        synergyDataTrialPostMeanList(:,dayI) = synergyPostTrialMean;
        synergyDataTrialSD   = std(synergyDataPost, 0, 2);
        hold on;
        plot(dataSelectPost(i).t, synergyPostTrialMean, 'b-');
        plot(dataSelectPost(i).t, synergyPostTrialMean+synergyDataTrialSD, 'c-');
        plot(dataSelectPost(i).t, synergyPostTrialMean-synergyDataTrialSD, 'c-');
        hold off;
        xlim([-15 15]);
        if(dayI==1), ylabel(focusList{i}); end
        if(i==1), title('Post (last day)'); end
        if(i==length(focusList)), xlabel('time %'); end
        
        % size(synergyDataPre)
        % size(synergyDataPost)

        %- Statistics -%
        corr = cosine_corr_interp(synergyPreTrialMean, synergyPostTrialMean);        
        pWil = ranksum(synergyPreTrialMean, synergyPostTrialMean); % Wilcoxon rank-sum test
        pPerm = PermTestTimeSeriesInterp(synergyDataPre, synergyDataPost);
        disp(['Pre-Post (last day) ' monkeyNameList{mIndx} ', Synergy ' focusList{i} ', Cosine Corr= ' num2str(corr), ', Wilcoxon P= ' num2str(pWil) ', Permutation test P=' num2str(pPerm)]);
        %--%
    end
end

%%
function c = cosine_corr_interp(x, y)
t1 = linspace(0,1,length(x));
t2 = linspace(0,1,length(y));

y2 = interp1(t2, y, t1);   % interp y so that its length becomes the same length with x
c = (x(:)' * y2(:)) / (norm(x) * norm(y2));
end

%%
function p = PermTestTimeSeriesInterp(x, y)
    nPerm = 10000;  % number of permulation
    targetT = 300;  % number of interpolated time points

    % ----------------------------------------------------
    % Step 1: Interpolate both datasets to 300 time points
    % ----------------------------------------------------
    [Tx1, nX] = size(x);
    [Ty1, nY] = size(y);
    % Original time axes
    tx = linspace(1, Tx1, Tx1);
    ty = linspace(1, Ty1, Ty1);
    % Interpolated time axes
    t_new_x = linspace(1, Tx1, targetT);
    t_new_y = linspace(1, Ty1, targetT);
    % Interpolate x
    x_interp = zeros(targetT, nX);
    for i = 1:nX
        x_interp(:, i) = interp1(tx, x(:, i), t_new_x, 'linear');
    end
    % Interpolate y
    y_interp = zeros(targetT, nY);
    for j = 1:nY
        y_interp(:, j) = interp1(ty, y(:, j), t_new_y, 'linear');
    end

    % ----------------------------------------------------
    % Step 2: Compute observed statistic
    %         (norm of difference between mean waveforms)
    % ----------------------------------------------------
    meanX = mean(x_interp, 2);  % 300 × 1
    meanY = mean(y_interp, 2);  % 300 × 1
    obsStat = norm(meanX - meanY);  % scalar

    % ----------------------------------------------------
    % Step 3: Permutation test
    % ----------------------------------------------------
    allData = [x_interp, y_interp];   % 300 × (nX + nY)
    nTotal = nX + nY;
    permStats = zeros(nPerm, 1);
    for k = 1:nPerm
        idx = randperm(nTotal);
        Xp = allData(:, idx(1:nX));        % permuted X group
        Yp = allData(:, idx(nX+1:end));    % permuted Y group
        meanXp = mean(Xp, 2);
        meanYp = mean(Yp, 2);
        permStats(k) = norm(meanXp - meanYp);  % statistic for this permutation
    end
    % p-value = proportion of permutations with stat >= observed stat
    p = mean(permStats >= obsStat);
end


%% Load Data from Mat files
function dataAll = loadMatData(calcType, monkeyName, condition, matDir, focusList, ttDate)

matPath  = [matDir '/synergyData_' monkeyName '_' condition '.mat'];
dataFile = load(matPath, 'synergyData', 'synergyWData', 'percentTime', 'expDates','nameList');
dataCell = dataFile.synergyData;
dataSynergyWCell = dataFile.synergyWData;

nList = length(focusList);
focusListNum = zeros(1,nList);
for i=1:nList
    focusListNum(i) = find(strcmp(focusList{i}, dataFile.nameList));
end
% dataMean = struct([]);
dataAll  = struct([]);
for i = 1:nList
    t = dataFile.percentTime;
    expDayString = cellstr(dataFile.expDates);
    expDates = datetime(expDayString, 'InputFormat','yyMMdd');
    daysDiff = days(expDates - ttDate);
    %
    meanData = zeros(length(t),length(expDates));
    allData  = cell(1, length(expDates));
    for dayI = 1:length(expDates)
        meanData(:,dayI) = mean(dataCell{focusListNum(i),dayI},1)';
        allData{dayI}    = (dataCell{focusListNum(i),dayI})';
    end
    dataAll(i).t         = t;
    dataAll(i).dataAll   = allData;
    dataAll(i).dataM     = meanData;
    dataAll(i).daysDiff  = daysDiff;
    %
    dataAll(i).wAll      = dataFile.synergyWData{i};
end
end
