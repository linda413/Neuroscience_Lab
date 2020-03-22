%% Separates the animals into groups by spike rate rather than treatment

clc; clear; close all;

% Grabs the animals of interest
%fileInfo = readAnimalcodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
fileInfo = readAnimalCodes(excelFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},{'LV'},1,[]);
analysisAnimals = fileInfo(:,1);

% Gets the bad channels
channels = readBadChannels(fileInfo,excelFile,2,'novel1');

% Makes saveDirectory if it doesn't exist
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Spiking/'];
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Indices of animals belonging to a group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Group colors
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};

%oddCohort = {'A53';'A54';'A55';'A56';'A57';'A58'};
%oddCohortIdx = ismember(analysisAnimals,oddCohort);

faceColors = cell(length(analysisAnimals),1);
faceColors(groupIdx{1}) = groupColor(1); faceColors(groupIdx{2}) = groupColor(2); faceColors(groupIdx{3}) = groupColor(3); faceColors(groupIdx{4}) = groupColor(4);
edgeColors = faceColors;
%edgeColors(oddCohortIdx) = {[0 0 0]};

% Spikes
%spikeTable = plotNORspikeTotals(excelFile,spikeDir,region,plotIT,saveDir)
spikeTable = plotNORspikeTotals(excelFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/','LV',0,[]);
spikeRate = spikeTable(:,6);
emptyIdx = cellfun(@isempty,spikeRate);
spikeRate(emptyIdx) = {NaN};

% P-episode
PepisodeAll = nan(length(analysisAnimals),30);
%avePepisode = nan(length(analysisAnimals),1);
for iLoad = 1:length(analysisAnimals)
    fileName = ['/Users/amberschedlbauer/Desktop/Gurkoff Lab/Pepisode/Pepisode_' analysisAnimals{iLoad} '_novel1_' channels{iLoad} '.mat'];
    if exist(fileName,'file') == 2
        load(fileName,'Pepisode','freqTimeWin','backgroundFrequencies','frequenciesOfInterest')
        %[~,idxF1] = min(abs(backgroundFrequencies - frequenciesOfInterest(1)));
        %[~,idxF1] = min(abs(backgroundFrequencies - 7));
        %[~,idxF2] = min(abs(backgroundFrequencies - frequenciesOfInterest(2)));
        %[~,idxF2] = min(abs(backgroundFrequencies - 8));
        %temp = nanmean(sum(freqTimeWin(idxF1:idxF2,:),1) > 0);
        %avePepisode(iLoad) = temp*100;
        PepisodeAll(iLoad,:) = Pepisode;
        clear Pepisode
    end
end
[~,idxF1] = min(abs(backgroundFrequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(backgroundFrequencies - frequenciesOfInterest(2)));
avePepisode = nanmean(PepisodeAll(:,idxF1:idxF2),2);

% Performance
[performance,~]  = readPerformance(excelFile,5,7,6,5:11,analysisAnimals);
performance(performance == 0) = NaN;
performance(performance == 1) = NaN;

%% Plot 2D scatter plots of spike rate, Pepisode, and performance

%figure('Position',[100 100 2000 500])

allVars = [[spikeRate{:}]' avePepisode performance];

% Spike Rate vs. Pepisode
figure('Position',[100 100 500 500])
%subplot(131)
for iGroup = 1:size(allVars,1)
    scatter(allVars(iGroup,1),allVars(iGroup,2),60,'filled','MarkerFaceColor',faceColors{iGroup},'MarkerEdgeColor',edgeColors{iGroup})
    hold on
end
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Spikes/Min')
ylabel('Pepisode (%)')
axis([0 max(allVars(:,1))+5 0 max(allVars(:,2))+5])
title('Spike Rate vs. Pepisode')
hold off
save2pdf([saveDirectory 'spikeCluster_2Dscatter1_LV'])

% Spike Rate vs. Performance
figure('Position',[100 100 500 500])
%subplot(132)
for iGroup = 1:size(allVars,1)
    scatter(allVars(iGroup,1),allVars(iGroup,3),60,'filled','MarkerFaceColor',faceColors{iGroup},'MarkerEdgeColor',edgeColors{iGroup})
    hold on
end
hline(0.5,'--k')
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Spikes/Min')
ylabel('Performance')
axis([0 max(allVars(:,1))+5 0 1])
title('Spike Rate vs. Performance')
hold off
save2pdf([saveDirectory 'spikeCluster_2Dscatter2_LV'])

% Pepisode vs. Performance
figure('Position',[100 100 500 500])
%subplot(133)
for iGroup = 1:size(allVars,1)
    scatter(allVars(iGroup,2),allVars(iGroup,3),60,'filled','MarkerFaceColor',faceColors{iGroup},'MarkerEdgeColor',edgeColors{iGroup})
    hold on
end
hline(0.5,'--k')
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Pepisode (%)')
ylabel('Performance')
axis([0 max(allVars(:,2))+5 0 1])
title('Pepisode vs. Performance')
hold off
save2pdf([saveDirectory 'spikeCluster_2Dscatter3_LV'])

%% Use K means to group high and low spike rate groups

X = allVars(:,1:2);
[kIdx,C] = kmeans(X,2,'Distance','cosine');

figure('Position',[100 100 500 500])
hold on
plot(X(kIdx==1,1),X(kIdx==1,2),'.','MarkerSize',20,'Color',[0 0 0])
plot(X(kIdx==2,1),X(kIdx==2,2),'.','MarkerSize',20,'Color',[0.5 0.5 0.5])
vline(2.5,'--r')
title('K Means Clusters')
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Spikes/Min')
ylabel('Pepisode(% time)')
axis([0 max(allVars(:,1))+5 0 max(allVars(:,2))+5])
box off
save2pdf([saveDirectory 'spikeCluster_KmeansClusters_LV'])

%% Scatter plots separated by K means group

%figure('Position',[100 100 2000 1000])

allXlabels = {'Spikes/Min';'Spikes/Min';'Pepisode (%)';'Spikes/Min';'Spikes/Min';'Pepisode (%)'};
allYlabels = {'Pepisode (%)';'Performance';'Performance';'Pepisode (%)';'Performance';'Performance'};
variablePair = {[1 2];[1 3];[2 3];[1 2];[1 3];[2 3]};

for iPlot = 1:6
    
    if iPlot < 4
        Kgroup = 1;
    else
        Kgroup = 2;
    end
    
    figure('Position',[100 100 500 500])
    %subplot(2,3,iPlot)
    hold on
    for iGroup = find(kIdx == Kgroup)'
        scatter(allVars(iGroup,variablePair{iPlot}(1)),allVars(iGroup,variablePair{iPlot}(2)),60,'filled','MarkerFaceColor',faceColors{iGroup},'MarkerEdgeColor',edgeColors{iGroup})
    end
    if ismember(iPlot,[2 3 5 6])
        hline(0.5,'--k')
    end
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    xlabel(allXlabels{iPlot})
    ylabel(allYlabels{iPlot})
    hold off
    
end

save2pdf([saveDirectory 'spikeCluster_KmeansClustersGroups_LV'])

%% Additional calculations

% Correlation
% a = allVars(kIdx == 1,2:3);
% [remove,~] = find(isnan(a));
% a(remove,:) = [];
% [rho,pval] = corr(a(:,1),a(:,2))

% Regression
% x = a(:,1); y = a(:,2);
% bls = regress(y,[ones(length(x),1) x]);
% brob = robustfit(x,y);
% scatter(x,y,'filled'); grid on; hold on
% plot(x,bls(1)+bls(2)*x,'k','LineWidth',2);
% plot(x,brob(1)+brob(2)*x,'Color',[0.5 0.5 0.5],'LineWidth',2)
