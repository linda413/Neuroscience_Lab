%% Longitudinal Changes in Pepisode

% Model: postpilostim04 pepisode and power

clc; clear; close all;

%% Get data

fileInfo = readAnimalCodes('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx',[],{'SHAM';'PILO'},{'NO'; 'STIM';'BURST'},{'LV'},1,[]);
analysisAnimals = fileInfo(:,1);

channels = readBadChannels(fileInfo,'/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx',2,session);

session = {'prepilostim01';'prepilostim02';'prepilostim03';'postpilostim04';'postpilostim16';'postpilostim31';'novel1';'novel2'};

frequencies = logspace(log10(2),log10(30),30);

PepisodeAll = nan(length(frequencies),size(fileInfo,1),length(session));

for iTime = 1:length(session)
    for iAnimal = 1:length(analysisAnimals)
        pEpisodeDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Pepisode/';
        pEpisodeFile = [pEpisodeDir 'Pepisode_' analysisAnimals{iAnimal} '_' session{iTime} '_' channels{iAnimal} '.mat'];
        if exist(pEpisodeFile,'file') == 2
            load(pEpisodeFile,'Pepisode')
            PepisodeAll(:,iAnimal,iTime) = Pepisode;
        end
    end
end

PepisodeAll(:,37,:) = NaN; % NaN out that one outlier PILO animal (A54)

%% Plot Pepisode

groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};

[~,idxF1] = min(abs(frequencies - 6));
[~,idxF2] = min(abs(frequencies - 10));

%%% All sessions, groups, individuals, frequencies
ctr = 1;
figure('Position',[100 100 2000 1000])
for iSession = 1:length(session)
    for iGroup = 1:length(groups)
        subplot(8,4,ctr)
        hold on
        plot(frequencies,PepisodeAll(:,groupIdx{iGroup},iSession),'Color',groupColor{iGroup},'LineWidth',1);
        plot(frequencies,nanmean(PepisodeAll(:,groupIdx{iGroup},iSession),2),'Color',groupColor{iGroup},'LineWidth',4);
        ylim([0 100])
        
        if iSession == 1
            title(groups{iGroup})
        end
        
        if iGroup == 1
            ylabel(session{iSession})
        end
        
        ctr = ctr + 1; 
    end
end
%%%

%%% Raw - Average, Peak, or 7.4 Theta Pepisode All Sessions Line Plot
aveTheta = squeeze(nanmean(PepisodeAll(idxF1:idxF2,:,:)));
peakTheta = squeeze(max(PepisodeAll(idxF1:idxF2,:,:)));
singleTheta = squeeze(PepisodeAll(15,:,:));
theta = aveTheta;

% Line
figure('Position',[100 100 2000 1000])
subplot(121)
hold on
for iGroup = 1:length(groups)
    plot(1:length(session),theta(groupIdx{iGroup},:),'Color',[groupColor{iGroup} 0.25],'LineWidth',1);
    plot(1:length(session),nanmean(theta(groupIdx{iGroup},:)),'Color',groupColor{iGroup},'LineWidth',4);
    ylim([0 50])
end
set(gca,'XTick',1:length(session),'XTickLabel',session)
ylabel('Average Pepisode')
grid on

%%%

%%% Compared to Pre - Average, Peak, or 7.4 Theta Pepisode Line Plot
newPepisode = cat(3,nanmean(PepisodeAll(:,:,1:3),3),PepisodeAll(:,:,4:end));
thetaPre1 = newPepisode./repmat(newPepisode(:,:,1),[1 1 6]);

aveTheta = squeeze(nanmean(thetaPre1(idxF1:idxF2,:,:)));
peakTheta = squeeze(max(thetaPre1(idxF1:idxF2,:,:)));
singleTheta = squeeze(thetaPre1(15,:,:));
theta = aveTheta;

% Line
%figure('Position',[100 100 1000 1000])
subplot(122)
hold on
for iGroup = 1:length(groups)
    plot(1:6,theta(groupIdx{iGroup},:),'Color',[groupColor{iGroup} 0.25],'LineWidth',1);
    plot(1:6,nanmean(theta(groupIdx{iGroup},:)),'Color',groupColor{iGroup},'LineWidth',4);
    ylim([0 3])
end
set(gca,'XTick',1:6,'XTickLabel',['pre'; session(4:end)])
ylabel('Average Pepisode')
grid on

%%%

%%% Scatter plot for each session
UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = [];
UDIs.YAXIS = [];
UDIs.COLORMAP = '';
UDIs.XLABEL = {''};
UDIs.YLABEL = 'Average Theta Pepisode';
UDIs.TITLE = '';
UDIs.COLORBARTITLE = '';
UDIs.SAVENAME = '';

dataCell = {aveTheta(groupIdx{1},1)',aveTheta(groupIdx{2},1)',aveTheta(groupIdx{3},1)',aveTheta(groupIdx{4},1)'};
[pval,~,~] = anova1(aveTheta(:,1),groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
UDIs.TITLE = ['1-Way ANOVA: p = ' num2str(pval)];
UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%% Read in NOR performance and spike rate from excel file
performance  = readPerformance('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx',5,7,analysisAnimals);
performance(performance == 0) = NaN;
performance(performance == 1) = NaN;

spikeTable = plotNORspikeTotals('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx','/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/','LV',0,[]);
spikeRate = spikeTable(:,6);
emptyIdx = cellfun(@isempty,spikeRate);
spikeRate(emptyIdx) = {NaN};
%%%
