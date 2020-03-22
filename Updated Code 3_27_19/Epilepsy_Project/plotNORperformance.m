%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot NOR (and Barnes) performance
% NOTE: Modify heavily if use!
% 
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Get Data

% Animal information excel file
excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Grabs the animals of interest
%fileInfo = readNORcodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(excelFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},{'LV'},1,[]);
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Indices of animals belong to a spike group
spikeRateThresh = 2.5;
spikeGroups = {'High Spike Rate';'Low Spike Rate'};
spikeTable = plotNORspikeTotals(excelFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
  if spikeTable{iSpike,end} >= spikeRateThresh % High spike rate group
      spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
  elseif spikeTable{iSpike,end} < spikeRateThresh % Low spike rate group
    spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
  end
end

% Read in NOR performance from excel file
[performanceNOR,performanceBarnes]  = readPerformance(excelFile,5,7,6,5:11,analysisAnimals);

% Remove any animals that did not interact with an object or interacted
% with only 1 object.
performanceNOR(performanceNOR == 0) = NaN;
performanceNOR(performanceNOR == 1) = NaN;

%% Plot Performance and Calcuate Stats

% Save directory for plots
saveDir = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Performance/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir) % Makes saveDirectory if it doesn't exist
end

% Set graphing UDIs as empty
UDIs.groups = {};
UDIs.groupColor = {};
UDIs.XAXIS = [];
UDIs.YAXIS = [];
UDIs.COLORMAP = '';
UDIs.XLABEL = {''};
UDIs.YLABEL = '';
UDIs.TITLE = {''};
UDIs.COLORBARTITLE = '';
UDIs.SAVENAME = '';

% Group parameters for plots
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};

%%% 4 Groups NOR %%%
dataCell = {performanceNOR(groupIdx{1})',performanceNOR(groupIdx{2})',performanceNOR(groupIdx{3})',performanceNOR(groupIdx{4})'};

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.YLABEL = 'NOR Performance';

[pval,~,~] = anova1(performanceNOR,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
%[results,means] = multcompare(stats,'CType','bonferroni');
UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDir 'performanceNOR_box_4groups'];

UDIs = plotAnalyses(dataCell,NaN,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 2 Groups NOR %%%
dataCell = {performanceNOR(spikeGroupIdx{1})',performanceNOR(spikeGroupIdx{2})'};

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.YLABEL = 'NOR Performance';

[~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDir 'performanceNOR_box_2groups'];

UDIs = plotAnalyses(dataCell,NaN,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 4 Groups Barnes %%%
dataCell = {performanceBarnes(groupIdx{1})',performanceBarnes(groupIdx{2})',performanceBarnes(groupIdx{3})',performanceBarnes(groupIdx{4})'};

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.YLABEL = '';

[pval,~,~] = anova1(performanceBarnes,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
%[results,means] = multcompare(stats,'CType','bonferroni');
UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDir 'performanceBarnes_box_4groups'];

UDIs = plotAnalyses(dataCell,NaN,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 2 Groups NOR %%%
dataCell = {performanceBarnes(spikeGroupIdx{1})',performanceBarnes(spikeGroupIdx{2})'};

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.YLABEL = '';

[~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDir 'performanceBarnes_box_2groups'];

UDIs = plotAnalyses(dataCell,NaN,UDIs,'boxPlotScatter');
clear data*
%%%

% Outlier detection
% mu = cellfun(@nanmean,dataCell);
% std = cellfun(@nanstd,dataCell);
% 
% for i = 1:4
%     group = (dataCell{i} - mu(i))/std(i)
%     disp(find(group > 2.5 | group < -2.5))
% end
