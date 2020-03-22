%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spike rate over timepoints across epilepsy experiment
%
% Amber Schedlbauer
% 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% UDIs

% Animal info -> Can either be excel file or filesToAnalyze_date.mat file
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Enter the path of where all the data is located
mainDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/';

% Read in excel file for bad channels
sheetNum = 2;

% Save directory for figures
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Spiking/'];

%% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},{'LV'},1,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
animals = fileInfo(:,1);
   
%% Get spikes for "pre" period and average across all 3 days

tasks = {'prepilostim01';'prepilostim02';'prepilostim03'};

for iTask = 1:length(tasks)
    
    channels = readBadChannels(fileInfo,codesFile,sheetNum,tasks{iTask});
    
    for iAnimal = 1:length(animals)
        
        if isempty(channels{iAnimal})
            spikeTotalPre(iAnimal,iTask) = NaN;
            spikeTimePre(iAnimal,iTask) = NaN;
        else
            load([mainDir 'spikeTotals_' animals{iAnimal} '_' tasks{iTask} '_' channels{iAnimal} '.mat'],'artifactReject','spikeTotal');
            spikeTotalPre(iAnimal,iTask) = spikeTotal;
            spikeTimePre(iAnimal,iTask) = sum(~artifactReject);
        end
        
        clear artifactReject spikeTotal
        
    end
    
    clear currentCol channels
    
end

spikeRatePre = nansum(spikeTotalPre,2)./(nansum(spikeTimePre,2)/60000);
%spikeRatePre = nan(length(animals),1);

%% Get spikes for each postpilostim days
tasks = {'postpilostim04';'postpilostim08';'postpilostim12';'postpilostim16';'postpilostim30';'postpilostim31'};

for iTask = 1:length(tasks)
    
    channels = readBadChannels(fileInfo,codesFile,sheetNum,tasks{iTask});

    for iAnimal = 1:length(animals)
        
        if isempty(channels{iAnimal})
            spikeRatePost(iAnimal,iTask) = NaN;
        else
            load([mainDir 'spikeTotals_' animals{iAnimal} '_' tasks{iTask} '_' channels{iAnimal} '.mat'],'spikeRate');
            spikeRatePost(iAnimal,iTask) = spikeRate;
        end
        
    end
    
    clear currentCol channels
    
end

%% Get spikes for NOR

%spikeTable = plotNORspikeTotals(excelFile,spikeDir,region,plotIT,saveDir)
spikeTable = plotNORspikeTotals(codesFile,mainDir,'LV',0,[]);
spikeRateNOR = spikeTable(:,6);
spikeRateNOR(cellfun(@isempty,spikeRateNOR)) = {NaN};
spikeRateNOR = cell2mat(spikeRateNOR);

return

%% Plot

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Indices of animals belong to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};

% Indices of animals belong to a spike group
spikeGroups = {'HSR';'LSR'};
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
  if spikeTable{iSpike,end} >= 2.5 % High spike rate group
      spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
  elseif spikeTable{iSpike,end} < 2.5 % Low spike rate group
    spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
  end
end
spikeGroupColor = {[0 0 0];[0 1 1]};

allSpikeRates = [spikeRatePre spikeRatePost spikeRateNOR];
%baselineNorm = (allSpikeRates - allSpikeRates(:,2))./allSpikeRates(:,2)*100;
%baselineNorm(allSpikeRates(:,1) == 0,:) = allSpikeRates(allSpikeRates(:,1) == 0,:);
baselineNorm = allSpikeRates;

% Individual animals
% hCmap = colormap(jet(length(animals)));
% figure('Position',[100 100 1500 500])
% for iPlot = 1:length(animals)
%     plot(baselineNorm(iPlot,:),'Color',hCmap(iPlot,:),'LineWidth',4)
%     hold on
%     pause
% end
% hold off
% set(gca,'XTick',1:(2+length(tasks)),'XTickLabel',[{'prepilostim'}; tasks; {'NOR'}],'FontSize',20)
% %ylabel('% Change from Prestim')
% save2pdf([saveDirectory 'spikeRateAllTimes_individualAnimals'])

% 4 treatment groups
figure('Position',[100 100 1500 500])
for iPlot = 1:length(groups)
    %plot(allSpikeRates(groupIdx{iPlot},:)','.-','Color',groupColor{iPlot},'MarkerSize',20)
    plot(baselineNorm(groupIdx{iPlot},:)','.-','Color',groupColor{iPlot},'MarkerSize',20)
    hold on
end
set(gca,'XTick',1:(2+length(tasks)),'XTickLabel',[{'prepilostim'}; tasks; {'NOR'}])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Spike Rate (spikes/min)')
%ylabel('% Change from Prestim')
save2pdf([saveDirectory 'spikeRateAllTimes_scatter_4groups'])

% 4 treatment group means
figure('Position',[100 100 1500 500])
for iPlot = 1:length(groups)
    %plot(allSpikeRates(groupIdx{iPlot},:)','.','Color',groupColor{iPlot},'MarkerSize',20)
    plot(nanmean(allSpikeRates(groupIdx{iPlot},:))','.-','Color',groupColor{iPlot},'MarkerSize',20,'LineWidth',4)
    hold on
end
set(gca,'XTick',1:(2+length(tasks)),'XTickLabel',[{'prepilostim'}; tasks; {'NOR'}])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Mean Spike Rate (spikes/min)')
save2pdf([saveDirectory 'spikeRateAllTimes_scatter_4groups_means'])

% 2 treatment groups
figure('Position',[100 100 1500 500])
for iPlot = 1:length(spikeGroups)
    %plot(allSpikeRates(groupIdx{iPlot},:)','.','Color',groupColor{iPlot},'MarkerSize',20)
    %plot(allSpikeRates(spikeGroupIdx{iPlot},:)','.-','Color',spikeGroupColor{iPlot},'MarkerSize',20)
    plot(baselineNorm(spikeGroupIdx{iPlot},:)','.-','Color',spikeGroupColor{iPlot},'MarkerSize',20)
    hold on
end
set(gca,'XTick',1:(2+length(tasks)),'XTickLabel',[{'prepilostim'}; tasks; {'NOR'}])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Spike Rate (spikes/min)')
%ylabel('% Change from Prestim')
save2pdf([saveDirectory 'spikeRateAllTimes_scatter_2groups'])

% % Change plot
% change = allSpikeRates - allSpikeRates(:,2);
% figure('Position',[100 100 500 1000])
% hold on
% for iPlot = 1:length(groups)
%     plot(ones(1,length(groupIdx{iPlot})),change(groupIdx{iPlot},3)','.','Color',groupColor{iPlot},'MarkerSize',20)
% end
% for iPlot = 1:length(spikeGroups)
%      plot(2*ones(1,length(spikeGroupIdx{iPlot})),change(spikeGroupIdx{iPlot},3)','.','Color',spikeGroupColor{iPlot},'MarkerSize',20)
% end
% hline(0,'k','')
% xlim([0 3])
% set(gca,'XTick',0:3,'XTickLabel',{' ','NOR','NOR',' '})
% set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
% %ylabel('Spike Rate (spikes/min)')
% ylabel('Difference in spike rate from postpilo04')
% save2pdf([saveDirectory 'spikeRateAllTimes_changeFrompostPilopilostim04'])

% Plots spike rates over time versus postpilostim04
figure('Position',[100 100 2000 1000])
hold on
colors = {'r','y','g','c','b','m'};
for iPlot = 1:6
    
    plot(spikeRatePost(:,iPlot),spikeRateNOR,'.','Color',colors{iPlot},'MarkerSize',40)
    text(spikeRatePost(:,iPlot),spikeRateNOR,animals)
    plot(0:30,0:30,'r')
    pause
end
xlabel('Spike Rate PPID04')
ylabel('Spike Rate NOR')

vars = [spikeRatePost(:,1) spikeRateNOR];
[rho,pval] = corr(vars(:,1),vars(:,2))