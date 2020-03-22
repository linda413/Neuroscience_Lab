%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Peri-spike histogram analyses for Gurkoff Lab
%
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input into pEpisode)
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% NOR tasks
session = {'novel1';'novel2'};

% Group to analyze
groupInput = {'SHAM';'PILO'};

% Treatment to analyze
treatmentInput = {'NO'; 'STIM';'BURST'};

% Region to analyze (enter only 1)
electrode = {'LV'};

% Enter the spike interval surrounding the peak identified in countInterictalSpikes_AS.m
spikeInterval = -100:250;

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/Spiking/'];

%% Loads spike files and grabs the marked spikes

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,1,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Indices of animals belong to a spike group
spikeGroups = {'High Spike Rate';'Low Spike Rate'};
spikeThreshold = 2.5;
spikeTable = plotNORspikeTotals(codesFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
    if spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
        spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    elseif spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
        spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    end
end

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session{1});

% Loop through the animals and create a peri-spike histogram for each
for iAnimal = 1:size(analysisAnimals)
    
    figure('Position',[100 100 2000 500])
    spikeSegNeg1 = []; spikeSegPos1 = []; spikeSegNeg2 = []; spikeSegPos2 = [];
    
    if ~isempty(channels{iAnimal})
        
        subplot(131)
        
        % Load file
        fileName = rdir([directoryPath 'spikeTotals_' analysisAnimals{iAnimal} '_' session{1} '_' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            
            load(fileName.name);
            spikeSegNeg1 = nan(1,length(spikeInterval));
            
            if ~isempty(spikeIntervalsNeg)
                for iSpike = 1:size(spikeIntervalsNeg,1)
                    temp = values(spikeIntervalsNeg(iSpike,1):spikeIntervalsNeg(iSpike,2));
                    temp(length(spikeInterval)) = NaN;
                    spikeSegNeg1 = [spikeSegNeg1; temp];
                end
                plot(spikeInterval,spikeSegNeg1')
                hold on
                plot(spikeInterval,nanmean(spikeSegNeg1,1),'Color','k','LineWidth',4)
            end
            
            spikeSegPos1 = nan(1,length(spikeInterval));
            if ~isempty(spikeIntervalsPos)
                for iSpike = 1:size(spikeIntervalsPos,1)
                    temp = values(spikeIntervalsPos(iSpike,1):spikeIntervalsPos(iSpike,2));
                    temp(length(spikeInterval)) = NaN;
                    spikeSegPos1 = [spikeSegPos1; temp];
                end
                
                plot(spikeInterval,spikeSegPos1')
                plot(spikeInterval,nanmean(spikeSegPos1,1),'Color','r','LineWidth',4)
            end
            
            clear artifactReject spikeIntervals* spikeRate spikeTotal supraInChan values temp
            
            ylim([-3 3])
            title('Novel 1')
            
        end
        
        subplot(132)
        
        % Load file
        fileName = rdir([directoryPath 'spikeTotals_' analysisAnimals{iAnimal} '_' session{2} '_' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            
            load(fileName.name);
            spikeSegNeg2 = nan(1,length(spikeInterval));
            
            if ~isempty(spikeIntervalsNeg)
                for iSpike = 1:size(spikeIntervalsNeg,1)
                    temp = values(spikeIntervalsNeg(iSpike,1):spikeIntervalsNeg(iSpike,2));
                    temp(351) = NaN;
                    spikeSegNeg2 = [spikeSegNeg2; temp];
                end
                plot(spikeInterval,spikeSegNeg2')
                hold on
                plot(spikeInterval,nanmean(spikeSegNeg2,1),'Color','k','LineWidth',4)
            end
            
            spikeSegPos2 = nan(1,length(spikeInterval));
            if ~isempty(spikeIntervalsPos)
                for iSpike = 1:size(spikeIntervalsPos,1)
                    temp = values(spikeIntervalsPos(iSpike,1):spikeIntervalsPos(iSpike,2));
                    temp(length(spikeInterval)) = NaN;
                    spikeSegPos2 = [spikeSegPos2; temp];
                end
                
                plot(spikeInterval,spikeSegPos2')
                plot(spikeInterval,nanmean(spikeSegPos2,1),'Color','r','LineWidth',4)
            end
            
            clear artifactReject spikeIntervals* spikeRate spikeTotal supraInChan values
            
            ylim([-3 3])
            title('Novel 2')
            
        end
        
        subplot(133)
        plot(spikeInterval,nanmean([spikeSegNeg1; spikeSegNeg2],1),'Color','k','LineWidth',4)
        hold on
        plot(spikeInterval,nanmean([spikeSegPos1; spikeSegPos2],1),'Color','r','LineWidth',4)
        ylim([-3 3])
        title('Novel 1 + Novel 2')
        
    end
    
    spikeNegAll{iAnimal} = [spikeSegNeg1; spikeSegNeg2];
    spikePosAll{iAnimal} = [spikeSegPos1; spikeSegPos2];
    
    %pause
    
    % Save the output plots
    if ~exist(saveDirectory,'dir')
        mkdir(saveDirectory)
    end
    %save2pdf([saveDirectory 'periSpike_' analysisAnimals{iAnimal}])
    
    close(gcf)
    
end

%% Plot Group Peri-spike histograms

% Colors for treatment groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {'r';'b'};

spikeNegAll = cellfun(@nanmean,spikeNegAll,'UniformOutput',false);
spikePosAll = cellfun(@nanmean,spikePosAll,'UniformOutput',false);

figure('Position',[100 100 1500 1000])
hold on
for iGroup = 1:length(groupIdx)
    data = spikeNegAll(groupIdx{iGroup})';
    data(cellfun(@numel,data) == 1) = [];
    data = cell2mat(data);
    plot(spikeInterval,data,'Color',groupColor{iGroup})
    hL(iGroup) = plot(spikeInterval,nanmean(data),'LineWidth',4,'Color',groupColor{iGroup});
    data = spikePosAll(groupIdx{iGroup})';
    data(cellfun(@numel,data) == 1) = [];
    data = cell2mat(data);
    plot(spikeInterval,data,'Color',groupColor{iGroup})
    plot(spikeInterval,nanmean(data),'LineWidth',4,'Color',groupColor{iGroup})
    %set(gca,'YTick',1:7:length(frequencies),'YTickLabel',frequencies(1:7:end),'FontSize',18,'FontWeight','bold')
    xlabel('Time (s)')
    ylabel('mV')
    title('Peri-spike Histogram')
    pause
end
hL = legend(hL,groups);
save2pdf([saveDirectory 'periSpikeHistogram'])
