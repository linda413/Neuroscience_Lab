%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Comparison of NOR1 and NOR2

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Animal information excel file
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Behavior from video excel file
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';

% Pepisode path
PepisodePath = '/Volumes/SP PHD U3/AmberS/Pepisode/';

% Region to analyze (enter only 1)
electrode = 'LV';

% Enter a vector of frequencies (log-space)
frequencies = logspace(log10(2),log10(30),50);

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious
windowOfInterest = [0 120];

% Sampling rate
sRate = 1000;

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/'];

%% Analysis

% Grabs the animals of interest
%fileInfo = readNORcodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(codesFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},{electrode},1,[]);
analysisAnimals = fileInfo(:,1);

% Read bad channels
channels1 = readBadChannels(fileInfo,codesFile,2,'novel1');
channels2 = readBadChannels(fileInfo,codesFile,2,'novel2');

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(frequencies - frequenciesOfInterest(2)));

% Allocations
[numWinNov1ITO,numWinNov1IBO,numWinNov2ITO,numWinNov2IBO] = deal(zeros(size(fileInfo,1),1));
[powerNov1ITO,powerNov1IBO,powerNov2ITO,powerNov2IBO] = deal(cell(1,length(analysisAnimals)));
[PepisodeNov1ITO,PepisodeNov1IBO,PepisodeNov2ITO,PepisodeNov2IBO] = deal(cell(1,length(analysisAnimals)));

% Loop through animals
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels1{iAnimal})  && ~isempty(channels2{iAnimal}) % Chcecks to make sure both channels are good
        
        % Load Pepisode -MAT file
        PepisodeMat1 = load([PepisodePath 'Pepisode_' analysisAnimals{iAnimal} '_novel1_' channels1{iAnimal} '.mat']);
        PepisodeMat2 = load([PepisodePath 'Pepisode_' analysisAnimals{iAnimal} '_novel2_' channels2{iAnimal} '.mat']);
        
        % Identify windows of interaction for Novel 1 & 2
        [~,timesNov1ITO] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{'novel1'},{'ITO'});
        [~,timesNov1IBO] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{'novel1'},{'IBO'});
        [~,timesNov2ITO] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{'novel2'},{'ITO'});
        [~,timesNov2IBO] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{'novel2'},{'IBO'});
        
        % Novel 1 ITO
        if ~isempty(timesNov1ITO{1})
            numWinNov1ITO(iAnimal) = size(timesNov1ITO{1},1);
            for iTime = 1:numWinNov1ITO(iAnimal)
                startTime = timesNov1ITO{1}(iTime,1)*sRate - 999;
                endTime = timesNov1ITO{1}(iTime,2)*sRate;
                powerNov1ITO{iAnimal}(:,:,iTime) = PepisodeMat1.power(:,startTime:endTime);
                PepisodeNov1ITO{iAnimal}(:,:,iTime) = PepisodeMat1.freqTimeWin(:,startTime:endTime);
            end
        else
            powerNov1ITO{iAnimal} = NaN(length(frequencies),3000);
            PepisodeNov1ITO{iAnimal} = NaN(length(frequencies),3000);
        end
        
        % Novel 1 IBO
        if ~isempty(timesNov1IBO{1})
            numWinNov1IBO(iAnimal) = size(timesNov1IBO{1},1);
            for iTime = 1:numWinNov1IBO(iAnimal)
                startTime = timesNov1IBO{1}(iTime,1)*sRate - 999;
                endTime = timesNov1IBO{1}(iTime,2)*sRate;
                powerNov1IBO{iAnimal}(:,:,iTime) = PepisodeMat1.power(:,startTime:endTime);
                PepisodeNov1IBO{iAnimal}(:,:,iTime) = PepisodeMat1.freqTimeWin(:,startTime:endTime);
            end
        else
            powerNov1IBO{iAnimal} = NaN(length(frequencies),3000);
            PepisodeNov1IBO{iAnimal} = NaN(length(frequencies),3000);
        end
        
        % Novel 2 ITO
        if ~isempty(timesNov2ITO{1})
            numWinNov2ITO(iAnimal) = size(timesNov2ITO{1},1);
            for iTime = 1:numWinNov2ITO(iAnimal)
                startTime = timesNov2ITO{1}(iTime,1)*sRate - 999;
                endTime = timesNov2ITO{1}(iTime,2)*sRate;
                powerNov2ITO{iAnimal}(:,:,iTime) = PepisodeMat2.power(:,startTime:endTime);
                PepisodeNov2ITO{iAnimal}(:,:,iTime) = PepisodeMat2.freqTimeWin(:,startTime:endTime);
            end
        else
            powerNov2ITO{iAnimal} = NaN(length(frequencies),3000);
            PepisodeNov2ITO{iAnimal} = NaN(length(frequencies),3000);
        end
        
        % Novel 2 IBO
        if ~isempty(timesNov2IBO{1})
            numWinNov2IBO(iAnimal) = size(timesNov2IBO{1},1);
            for iTime = 1:numWinNov2IBO(iAnimal)
                startTime = timesNov2IBO{1}(iTime,1)*sRate - 999;
                endTime = timesNov2IBO{1}(iTime,2)*sRate;
                powerNov2IBO{iAnimal}(:,:,iTime) = PepisodeMat2.power(:,startTime:endTime);
                PepisodeNov2IBO{iAnimal}(:,:,iTime) = PepisodeMat2.freqTimeWin(:,startTime:endTime);
            end
        else
            powerNov2IBO{iAnimal} = NaN(length(frequencies),3000);
            PepisodeNov2IBO{iAnimal} = NaN(length(frequencies),3000);
        end
        
    else
        % Need to put NaNs in cells of "bad" animals so can perform cellfun
        [powerNov1ITO{iAnimal},powerNov1IBO{iAnimal},powerNov2ITO{iAnimal},powerNov2IBO{iAnimal}] = deal(NaN(length(frequencies),3000));
        [PepisodeNov1ITO{iAnimal},PepisodeNov1IBO{iAnimal},PepisodeNov2ITO{iAnimal},PepisodeNov2IBO{iAnimal}] = deal(NaN(length(frequencies),3000));
    end
    
end

% Identify which object has been replaced in Novel 2
[~,novelObject,~] = xlsread(codesFile,4);
topAnimals = novelObject(strcmp(novelObject(:,5),'TOP'),1); 
topIdx = ismember(analysisAnimals,topAnimals);
bottomAnimals = novelObject(strcmp(novelObject(:,5),'BOTTOM'),1);
bottomIdx = ismember(analysisAnimals,bottomAnimals);

% Power Analysis
avePowerNov1ITO = cellfun(@(x) nanmean(x,3),powerNov1ITO,'UniformOutput',false);
avePowerNov1IBO = cellfun(@(x) nanmean(x,3),powerNov1IBO,'UniformOutput',false);
avePowerNov2ITO = cellfun(@(x) nanmean(x,3),powerNov2ITO,'UniformOutput',false);
avePowerNov2IBO = cellfun(@(x) nanmean(x,3),powerNov2IBO,'UniformOutput',false);

avePowerNov1ITO = reshape([avePowerNov1ITO{:}],[length(frequencies) 3000 length(analysisAnimals)]);
avePowerNov1IBO = reshape([avePowerNov1IBO{:}],[length(frequencies) 3000 length(analysisAnimals)]);
avePowerNov2ITO = reshape([avePowerNov2ITO{:}],[length(frequencies) 3000 length(analysisAnimals)]);
avePowerNov2IBO = reshape([avePowerNov2IBO{:}],[length(frequencies) 3000 length(analysisAnimals)]);

avePowerITO = 10*log10(avePowerNov2ITO ./ avePowerNov1ITO);
avePowerIBO = 10*log10(avePowerNov2IBO ./ avePowerNov1IBO);

avePowerNew = NaN(length(frequencies),3000,length(analysisAnimals));
avePowerNew(:,:,topIdx) = avePowerITO(:,:,topIdx);
avePowerNew(:,:,bottomIdx) = avePowerIBO(:,:,bottomIdx);

avePowerOld = NaN(length(frequencies),3000,length(analysisAnimals));
avePowerOld(:,:,topIdx) = avePowerIBO(:,:,topIdx);
avePowerOld(:,:,bottomIdx) = avePowerITO(:,:,bottomIdx);

% Pepisode Analysis
avePepisodeNov1ITO = cellfun(@(x) nanmean(x,3),PepisodeNov1ITO,'UniformOutput',false);
avePepisodeNov1IBO = cellfun(@(x) nanmean(x,3),PepisodeNov1IBO,'UniformOutput',false);
avePepisodeNov2ITO = cellfun(@(x) nanmean(x,3),PepisodeNov2ITO,'UniformOutput',false);
avePepisodeNov2IBO = cellfun(@(x) nanmean(x,3),PepisodeNov2IBO,'UniformOutput',false);

% Indices of animals belong to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Indices of animals belong to a spike group
% Indices of animals belong to a spike group
spikeGroups = {'High Spike Rate';'Low Spike Rate'};
spikeThreshold = 2.5;
spikeTable = plotNORspikeTotals(codesFile,'/Volumes/SP PHD U3/AmberS/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
  if spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
      spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
  elseif spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
    spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
  end
end

%% PLOT

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
windowTime = -1:(1/sRate):1.999;

%%% New Objects %%%

%%% Average spectrogram for 4 groups
data = cat(3,nanmean(avePowerNew(:,:,groupIdx{1}),3),nanmean(avePowerNew(:,:,groupIdx{2}),3),nanmean(avePowerNew(:,:,groupIdx{3}),3),nanmean(avePowerNew(:,:,groupIdx{4}),3));
colorbarLims = [min(data(:)) max(data(:))];
figure('Position',[100 100 1500 1000])
for iGroup = 1:length(groupIdx)
   ax1(iGroup) = subplot(2,2,iGroup);
   imagesc(windowTime,1:length(frequencies),data(:,:,iGroup))
   hline([idxF1 idxF2],{'--k','--k'})
   axis xy
   set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
   cH = colorbar;
   cH.Label.String = 'dB';
   caxis(colorbarLims)
   xlabel('Time (s)')
   ylabel('Frequencies (Hz)')
   title(['Power ' groups{iGroup}])
end
linkaxes(ax1)
save2pdf([saveDirectory 'spectrogram_newInteracting_4groups_' electrode])
%%%

%%% Average spectrogram for 2 groups
data = cat(3,nanmean(avePowerNew(:,:,spikeGroupIdx{1}),3),nanmean(avePowerNew(:,:,spikeGroupIdx{2}),3));
colorbarLims = [min(data(:)) max(data(:))];
figure('Position',[100 100 1500 500])
for iGroup = 1:length(spikeGroupIdx)
   ax2(iGroup) = subplot(1,2,iGroup);
   imagesc(windowTime,1:length(frequencies),data(:,:,iGroup))
   hline([idxF1 idxF2],{'--k','--k'})
   axis xy
   set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
   cH = colorbar;
   cH.Label.String = 'dB';
   caxis(colorbarLims)
   xlabel('Time (s)')
   ylabel('Frequencies (Hz)')
   title(['Power ' spikeGroups{iGroup}])
end
linkaxes(ax2)
save2pdf([saveDirectory 'spectrogram_newInteracting_2groups_' electrode])
%%%

%%% Old Objects %%%

%%% Average spectrogram for 4 groups
data = cat(3,nanmean(avePowerOld(:,:,groupIdx{1}),3),nanmean(avePowerOld(:,:,groupIdx{2}),3),nanmean(avePowerOld(:,:,groupIdx{3}),3),nanmean(avePowerOld(:,:,groupIdx{4}),3));
colorbarLims = [min(data(:)) max(data(:))];
figure('Position',[100 100 1500 1000])
for iGroup = 1:length(groupIdx)
   ax1(iGroup) = subplot(2,2,iGroup);
   imagesc(windowTime,1:length(frequencies),data(:,:,iGroup))
   hline([idxF1 idxF2],{'--k','--k'})
   axis xy
   set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
   cH = colorbar;
   cH.Label.String = 'dB';
   caxis(colorbarLims)
   xlabel('Time (s)')
   ylabel('Frequencies (Hz)')
   title(['Power ' groups{iGroup}])
end
linkaxes(ax1)
save2pdf([saveDirectory 'spectrogram_oldInteracting_4groups_' electrode])
%%%

%%% Average spectrogram for 2 groups
data = cat(3,nanmean(avePowerOld(:,:,spikeGroupIdx{1}),3),nanmean(avePowerOld(:,:,spikeGroupIdx{2}),3));
colorbarLims = [min(data(:)) max(data(:))];
figure('Position',[100 100 1500 500])
for iGroup = 1:length(spikeGroupIdx)
   ax2(iGroup) = subplot(1,2,iGroup);
   imagesc(windowTime,1:length(frequencies),data(:,:,iGroup))
   hline([idxF1 idxF2],{'--k','--k'})
   axis xy
   set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
   cH = colorbar;
   cH.Label.String = 'dB';
   caxis(colorbarLims)
   xlabel('Time (s)')
   ylabel('Frequencies (Hz)')
   title(['Power ' spikeGroups{iGroup}])
end
linkaxes(ax2)
save2pdf([saveDirectory 'spectrogram_oldInteracting_2groups_' electrode])
%%%