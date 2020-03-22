function spikeTable = plotNORspikeTotals(excelFile,spikeDir,region,plotIT,saveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% spikeTable = plotNORspikeTotals(excelFile,spikeDir,region,plotIT,saveDir)
%
% Displays the spike counts for the groups
%
% Inputs:
%   excelFile - File that contains the NOR animal codes
%   Ex: excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%
%   spikeDir - Location of the spike total files
%   Ex: spikeDir = '/home/amschedl/Desktop/NOR/Spikes/';
%       spikeDir = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Spikes/';
%       spikeDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/';
%
%   region
%   Ex: region = 'LD';
%
%   plotIT - Enter 1 to output scatter box plots or 0 for no plots
%   Ex: plotIT = 1;
%
%   saveDir - Enter in the directory to save the figures if plotIT = 1
%   Ex: saveDir = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Spiking/'];
%
% Output:
%   spikeTable - table of animal and spike rate information
%
% Amber Schedlbauer - 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify the files of interest
fileInfo = readAnimalCodes(excelFile,[],{'SHAM';'PILO'},{'NO'; 'STIM';'BURST'},region,1,[]);

% Identifies which channels for the identified animals
channels1 = readBadChannels(fileInfo,excelFile,2,'novel1');
channels2 = readBadChannels(fileInfo,excelFile,2,'novel2');

% Allocation
spikeTable = cell(size(fileInfo,1),6);

for iAnimal = 1:size(fileInfo,1)
    
    fileName1 = [spikeDir 'spikeTotals_' fileInfo{iAnimal,1} '_novel1_' channels1{iAnimal} '.mat'];
    fileName2 = [spikeDir 'spikeTotals_' fileInfo{iAnimal,1} '_novel2_' channels2{iAnimal} '.mat'];
    
    normal1 = NaN; normal2 = NaN;
    supra1 = NaN; supra2 = NaN;
    both1 = NaN; both2 = NaN;
    totalTime1 = NaN; totalTime2 = NaN;
    
    if exist(fileName1,'file') == 2
        
        load(fileName1)
        
        normal1 = spikeTotal;
        supra1 = supraInChan;
        both1 = normal1 + supra1;
        totalTime1 = sum(~artifactReject)/60000;
        
        clear artifactReject supraInChan spikeTotal spikeRate
        
    end
    
    if exist(fileName2,'file') == 2
        
        load(fileName2)
        
        normal2 = spikeTotal;
        supra2 = supraInChan;
        both2 = normal2 + supra2;
        totalTime2 = sum(~artifactReject)/60000;
        
        clear artifactReject supraInChan spikeTotal spikeRate
        
    end
    
    if isnan(normal1) && isnan(normal2)
        continue
    end
    
    spikeRate = nansum([normal1 normal2]) / nansum([totalTime1 totalTime2]);
    supraSpikeRate = nansum([supra1 supra2]) / nansum([totalTime1 totalTime2]);
    bothRate = nansum([both1 both2]) / nansum([totalTime1 totalTime2]);
    spikeTable(iAnimal,:) = [fileInfo(iAnimal,1:3) spikeRate supraSpikeRate bothRate];
    
end


% Remove empty rows
%spikeTable(~sum(~cellfun(@isempty,spikeTable),2),:) = [];

% Sort table to group animals from same group and treatment
%[~,sortIdx] = sortrows(spikeTable(:,2:3));
%spikeTable = spikeTable(sortIdx,:);

%%% Plot results %%%
if plotIT
    
    % Makes saveDirectory if it doesn't exist
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    
    %% Scatter plots %%
    
    %     % OLD CODE that doesn't use plotBoxScatter -> Modify heavily if used
    %     figure('Position',[100 100 500 500]);
    %     scatter(ones(1,length(shamIdx)),[spikeTable{shamIdx,4}],'b','filled')
    %     hold on
    %     scatter(2*ones(1,length(piloIdx)),[spikeTable{piloIdx,4}],'r','filled')
    %     scatter(3*ones(1,length(stimIdx)),[spikeTable{stimIdx,4}],'g','filled')
    %     scatter(4*ones(1,length(burstIdx)),[spikeTable{burstIdx,4}],'m','filled')
    %     scatter(1,nanmean([spikeTable{shamIdx,4}]),100,'k','filled')
    %     scatter(2,nanmean([spikeTable{piloIdx,4}]),100,'k','filled')
    %     scatter(3,nanmean([spikeTable{stimIdx,4}]),100,'k','filled')
    %     scatter(4,nanmean([spikeTable{burstIdx,4}]),100,'k','filled')
    %     set(gca,'XTick',1:6,'XTickLabel',{'SHAM','PILO','STIM','BURST',''},'FontSize',20,'FontWeight','bold')
    %     xlim([0 5])
    %     ylabel('Spikes/Min')
    %     title('Normal Threshold Spikes')
    %     save2pdf([saveDir 'NormalSpikeCount_scatter'])
    %
    %     figure('Position',[100 100 500 500]);
    %     scatter(ones(1,length(shamIdx)),[spikeTable{shamIdx,5}],'b','filled')
    %     hold on
    %     scatter(2*ones(1,length(piloIdx)),[spikeTable{piloIdx,5}],'r','filled')
    %     scatter(3*ones(1,length(stimIdx)),[spikeTable{stimIdx,5}],'g','filled')
    %     scatter(4*ones(1,length(burstIdx)),[spikeTable{burstIdx,5}],'m','filled')
    %     scatter(1,nanmean([spikeTable{shamIdx,5}]),80,'k','filled')
    %     scatter(2,nanmean([spikeTable{piloIdx,5}]),80,'k','filled')
    %     scatter(3,nanmean([spikeTable{stimIdx,5}]),80,'k','filled')
    %     scatter(4,nanmean([spikeTable{burstIdx,5}]),80,'k','filled')
    %     set(gca,'XTick',1:6,'XTickLabel',{'SHAM','PILO','STIM','BURST',''},'FontSize',20,'FontWeight','bold')
    %     xlim([0 5])
    %     ylabel('Spikes/Min')
    %     title('Supra Threshold Spikes')
    %     save2pdf([saveDir 'SupraSpikeCount_scatter'])
    %
    %     figure('Position',[100 100 500 500]);
    %     scatter(ones(1,length(shamIdx)),[spikeTable{shamIdx,6}],'b','filled')
    %     hold on
    %     scatter(2*ones(1,length(piloIdx)),[spikeTable{piloIdx,6}],'r','filled')
    %     scatter(3*ones(1,length(stimIdx)),[spikeTable{stimIdx,6}],'g','filled')
    %     scatter(4*ones(1,length(burstIdx)),[spikeTable{burstIdx,6}],'m','filled')
    %     scatter(1,nanmean([spikeTable{shamIdx,6}]),80,'k','filled')
    %     scatter(2,nanmean([spikeTable{piloIdx,6}]),80,'k','filled')
    %     scatter(3,nanmean([spikeTable{stimIdx,6}]),80,'k','filled')
    %     scatter(4,nanmean([spikeTable{burstIdx,6}]),80,'k','filled')
    %     set(gca,'XTick',1:6,'XTickLabel',{'SHAM','PILO','STIM','BURST',''},'FontSize',20,'FontWeight','bold')
    %     xlim([0 5])
    %     ylabel('Spikes/Min')
    %     title('Normal + Supra Threshold Spikes')
    %     save2pdf([saveDir 'BothSpikeCount_scatter'])
    
    %% Box plots - 4 groups %%
    
    shamIdx = find(strcmp(spikeTable(:,2),'SHAM'));
    piloIdx = find(strcmp(spikeTable(:,2),'PILO') & strcmp(spikeTable(:,3),'NO'));
    stimIdx = find(strcmp(spikeTable(:,2),'PILO') & strcmp(spikeTable(:,3),'STIM'));
    burstIdx = find(strcmp(spikeTable(:,2),'PILO') & strcmp(spikeTable(:,3),'BURST'));
    groups = {'SHAM','PILO','STIM','BURST'};
    groupColor = {'b';'r';'g';'m'};
    
    %     spikeCell = {[spikeTable{shamIdx,4}],[spikeTable{piloIdx,4}],[spikeTable{stimIdx,4}],[spikeTable{burstIdx,4}]};
    %     plotBoxScatter(spikeCell,groups,groupColor)
    %     set(gca,'FontSize',20,'FontWeight','bold')
    %     ylabel('Spikes/Min')
    %     set(gca,'FontSize',20,'FontWeight','bold')
    %     title('Normal Spikes')
    %     %save2pdf([saveDir 'NormalSpikeCount_box_' region])
    
    %     spikeCell = {[spikeTable{shamIdx,5}],[spikeTable{piloIdx,5}],[spikeTable{stimIdx,5}],[spikeTable{burstIdx,5}]};
    %     plotBoxScatter(spikeCell,groups,groupColor)
    %     set(gca,'FontSize',20,'FontWeight','bold')
    %     ylabel('Spikes/Min')
    %     set(gca,'FontSize',20,'FontWeight','bold')
    %     title('Supra Threshold Spikes')
    %     %save2pdf([saveDir 'SupraSpikeCount_box_' region])
    
    figure('Position',[100 100 750 500])
    spikeCell = {[spikeTable{shamIdx,6}],[spikeTable{piloIdx,6}],[spikeTable{stimIdx,6}],[spikeTable{burstIdx,6}]};
    plotBoxScatter(spikeCell,groups,groupColor)
    set(gca,'FontSize',20,'FontWeight','bold')
    ylabel('Spikes/Min')
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    %title('Total Number of Spikes')
    
    % Unbalanced 1-way ANOVA for significance for total spike counts
    groupLabel = cell(1,size(fileInfo,1));
    groupLabel(shamIdx) = {'SHAM'}; groupLabel(piloIdx) = {'PILO'}; groupLabel(stimIdx) = {'STIM'}; groupLabel(burstIdx) = {'BURST'};
    totalCounts = spikeTable(:,6);
    totalCounts(cellfun(@isempty,spikeTable(:,6))) = {NaN};
    [p,~,~] = anova1(cell2mat(totalCounts),groupLabel,'off');
    title(['1-Way ANOVA: p = ' num2str(p)])
    
    %[results,means] = multcompare(stats,'CType','bonferroni')
    
    save2pdf([saveDir 'spikeRate_box_4groups_' region])
    
    %     % DUH YOU DONT NEED THIS
    %     %% Box plots - 2 groups %%
    %
    %     % Indices of animals belong to a spike group
    %     spikeGroups = {'High Spike Rate';'Low Spike Rate'};
    %     spikeGroupIdx = cell(1,2);
    %     for iSpike = 1:numel(spikeTable(:,end))
    %         if spikeTable{iSpike,end} >= 2.5 % High spike rate group
    %             spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    %         elseif spikeTable{iSpike,end} < 2.5 % Low spike rate group
    %             spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    %         end
    %     end
    %
    %     spikeCell = {[spikeTable{spikeGroupIdx{1},6}],[spikeTable{spikeGroupIdx{2},6}]};
    %     plotBoxScatter(spikeCell,spikeGroups,[])
    %     set(gca,'FontSize',20,'FontWeight','bold')
    %     ylabel('Spikes/Min')
    %     %title('Total Number of Spikes')
    %     save2pdf([saveDir 'spikeRate_box_2groups' region])
    
end

end

