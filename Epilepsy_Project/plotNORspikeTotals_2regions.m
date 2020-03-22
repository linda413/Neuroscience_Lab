%%% Script to plot NOR spike rates for RV and LV

clc; clear; close all;

excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
spikeDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/';
saveDir = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Spiking/'];

fileInfo = readAnimalCodes(excelFile,[],{'SHAM';'PILO'},{'NO'; 'STIM';'BURST'},{'RD'},1,[]);

spikeTableLV = plotNORspikeTotals(excelFile,spikeDir,'LV',0,saveDir);
spikeTableRV = plotNORspikeTotals(excelFile,spikeDir,'RV',0,saveDir);
spikeTableLD = plotNORspikeTotals(excelFile,spikeDir,'LD',0,saveDir);
spikeTableRD = plotNORspikeTotals(excelFile,spikeDir,'RD',0,saveDir);

spikeTableLV(cellfun(@isempty,spikeTableLV(:,6)),6) = {NaN};
spikeTableRV(cellfun(@isempty,spikeTableRV(:,6)),6) = {NaN};
spikeTableLD(cellfun(@isempty,spikeTableLD(:,6)),6) = {NaN};
spikeTableRD(cellfun(@isempty,spikeTableRD(:,6)),6) = {NaN};

data = [cell2mat(spikeTableLV(:,6)) cell2mat(spikeTableRV(:,6)) cell2mat(spikeTableRD(:,6)) cell2mat(spikeTableLD(:,6))];

figure('Position',[100 100 1000 1000])
for iPlot = 1:size(data,1) % 27
    
    if isempty(spikeTableLV{iPlot,1}) && isempty(spikeTableLD{iPlot,1}) && isempty(spikeTableRV{iPlot,1}) && isempty(spikeTableRD{iPlot,1})
       continue
    end
    
    % Animals that were originally labeled as high and low spikers
    if spikeTableLV{iPlot,6} > 2.5
        color = [0 0 0]; %[0.5 0.5 0.5 0.5]; % 
    elseif isnan(spikeTableLV{iPlot,6})
       continue
    else
        color = [0.5 0.5 0.5 0.5]; % [0 0 1]; %
    end
    
%     % Animals that had seizures in the long-term recordings
%     if ismember(fileInfo{iPlot,1},{'A42','A49','A50','A58','A72','A74'})
%        color = 'k';
%     end
    
    
    hP = plot(1:4,data(iPlot,:),'.-','Color',color,'LineWidth',4,'MarkerSize',30);
    hold on
    
end

hline(2.5,'--m')
%xlim([0.75 2.25])
set(gca,'XTick',1:4,'XTickLabel',{'LV','RV','RD','LD'})
ylabel('Spike Rate (spikes/min)')
xlabel('Region')
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')

% Makes saveDirectory if it doesn't exist
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
save2pdf([saveDir 'spikeRate_2regions_subset'])