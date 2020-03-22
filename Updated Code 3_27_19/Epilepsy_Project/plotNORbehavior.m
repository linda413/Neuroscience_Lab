%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to quantify time spent in a specific behavior during NOR
% NOTE: Modify heavily if use!
% 
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Find % time in each behavior

videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
[num,txt,raw] = xlsread(videoAnalysisFile,1,'C3:DS102'); % Read in specific rows and columns

animalTask = txt(:,1);
animals = cellfun(@(x) strtok(x,'_'),animalTask,'UniformOutput',false);
novel1Idx = find(~cellfun(@isempty,strfind(animalTask,'novel1')));
novel2Idx = find(cellfun(@isempty,strfind(animalTask,'novel1')));

behaviors = txt(:,2:end);

noData = sum(strcmp(behaviors,'NO DATA'),2);

Mbehav = sum((strcmp(behaviors,{'M'}) + strcmp(behaviors,{'ATO'}) + strcmp(behaviors,{'ABO'}) + strcmp(behaviors,{'LTO'}) + strcmp(behaviors,{'LBO'})),2);

NMbehav = sum(strcmp(behaviors,'NM'),2);

Sbehav = sum((strcmp(behaviors,'S') + strcmp(behaviors,{'SNOT'}) + strcmp(behaviors,{'SNOB'})),2);

Ibehav = sum(contains(behaviors,'I'),2);

totalTime = 120 - noData; % Data for only first 120 seconds

moving = Mbehav./totalTime;

notMovingSitting = (NMbehav + Sbehav)./totalTime;

interacting = Ibehav./totalTime;

%% Stats and plots

% Grabs the animals of interest
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(codesFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},{'LV'},1,[]);

% Indices of animals belonging to a treatment group
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Remove animals that did not interact with an object or interacted with only 1 object
performance  = readNORperformance(codesFile,5,fileInfo(:,1));
removeIdx = find(performance == 0 | performance == 1 | isnan(performance));

newRemoveIdx = [];
for iRemove = 1:length(removeIdx)
   newRemoveIdx =  [newRemoveIdx; find(strcmp(animals,fileInfo{removeIdx(iRemove),1}))];
end

moving(newRemoveIdx) = NaN; notMovingSitting(newRemoveIdx) = NaN; interacting(newRemoveIdx) = NaN;

% ANOVAs
task = cell(length(animals),1);
task(novel1Idx) = {'novel1'}; task(novel2Idx) = {'novel2'};

group = cell(length(animals)/2,2);
group(groupIdx{1},:) = {'SHAM'}; group(groupIdx{2},:) = {'PILO'}; group(groupIdx{3},:) = {'STIM'}; group(groupIdx{4},:) = {'BURST'};
group = group'; group = group(:);

% Save information
% Makes saveDirectory if it doesn't exist
saveDir = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Behavior/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%%% Moving
[p,tbl,stats] = anovan(moving,{group,task});
% Note: The results variable returns the means of the groups collapsed
% across task AND the standard errors (not the means of the groups for each
% task)!
[results,means] = multcompare(stats,'CType','bonferroni')

moving1 = moving(novel1Idx); moving2 = moving(novel2Idx);

for iGroup = 1:4 % 2 groups of 4 parameters
    
    data(1,iGroup) = nanmean(moving1(groupIdx{iGroup}));
    err(1,iGroup) = nanstd(moving1(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    data(2,iGroup) = nanmean(moving2(groupIdx{iGroup}));
    err(2,iGroup) = nanstd(moving2(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    
end

figure('Position',[100 100 1000 750])
h = barwitherr(err,data);
set(gca,'XTickLabel',{'Novel 1','Novel 2'})
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend('SHAM','PILO','STIM','BURST')
ylabel('Proportion of Time')
set(h(1),'FaceColor',[0 0 1],'EdgeColor',[0 0 1])
set(h(2),'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
set(h(3),'FaceColor',[0 1 0],'EdgeColor',[0 1 0])
set(h(4),'FaceColor',[1 0 1],'EdgeColor',[1 0 1])
title('Moving: Main Effects of Group (p = 0.012 PILO-STIM) and Session (p = 0.0017)')
save2pdf([saveDir 'behavior_moving_4groups'])

%%% NotMoving/Sitting
[p,tbl,stats] = anovan(notMovingSitting,{group,task});
[results,means] = multcompare(stats,'CType','bonferroni')

notMovingSitting1 = notMovingSitting(novel1Idx); notMovingSitting2 = notMovingSitting(novel2Idx);
for iGroup = 1:4 % 2 groups of 4 parameters
    
    data(1,iGroup) = nanmean(notMovingSitting1(groupIdx{iGroup}));
    err(1,iGroup) = nanstd(notMovingSitting1(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    data(2,iGroup) = nanmean(notMovingSitting2(groupIdx{iGroup}));
    err(2,iGroup) = nanstd(notMovingSitting2(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    
end

figure('Position',[100 100 1000 750])
h = barwitherr(err,data);
set(gca,'XTickLabel',{'Novel 1','Novel 2'})
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend('SHAM','PILO','STIM','BURST')
ylabel('Proportion of Time')
set(h(1),'FaceColor',[0 0 1],'EdgeColor',[0 0 1])
set(h(2),'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
set(h(3),'FaceColor',[0 1 0],'EdgeColor',[0 1 0])
set(h(4),'FaceColor',[1 0 1],'EdgeColor',[1 0 1])
title('Not Moving/Sitting: No Main Effects')
save2pdf([saveDir 'behavior_notMovingSitting_4groups'])

%%% Interacting
[p,tbl,stats] = anovan(interacting,{task,group});
[results,means] = multcompare(stats,'CType','bonferroni')

interacting1 = interacting(novel1Idx); interacting2 = interacting(novel2Idx);

for iGroup = 1:4 % 2 groups of 4 parameters
    
    data(1,iGroup) = nanmean(interacting1(groupIdx{iGroup}));
    err(1,iGroup) = nanstd(interacting1(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    data(2,iGroup) = nanmean(interacting2(groupIdx{iGroup}));
    err(2,iGroup) = nanstd(interacting2(groupIdx{iGroup}))/sqrt(numel(groupIdx{iGroup}));
    
end

figure('Position',[100 100 1000 750])
h = barwitherr(err,data);
set(gca,'XTickLabel',{'Novel 1','Novel 2'})
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend('SHAM','PILO','STIM','BURST')
ylabel('Proportion of Time')
set(h(1),'FaceColor',[0 0 1],'EdgeColor',[0 0 1])
set(h(2),'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
set(h(3),'FaceColor',[0 1 0],'EdgeColor',[0 1 0])
set(h(4),'FaceColor',[1 0 1],'EdgeColor',[1 0 1])
title('Interacting: Main Effect of Session (p = 0.005)')
save2pdf([saveDir 'behavior_interacting_4groups'])