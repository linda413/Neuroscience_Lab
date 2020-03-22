%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coherence Analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Animal information excel file
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'novel1';

% Group to analyze
groupInput = {'SHAM';'PILO'};

% Treatment to analyze
treatmentInput = {'NO'; 'STIM';'BURST'};

% Regions to analyze (enter only 2)
electrodes = {'LV','LPFC'};

% Enter the frequency range for coherence anlaysis
freqRange = [2 30];

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine pEpisode
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
%saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/'];
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrodes{1} '/'];

% Calculate power for specific behaviors
addSpecificBehaviors = 1;
    videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
    %behavior = {'NM';'S'};
    %behavior = {'M'};
    behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate coherence

% Grabs the animals of interest
%fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{1},1,[]);
channels1 = readBadChannels(fileInfo,codesFile,2,session);
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{2},1,[]);
channels2 = readBadChannels(fileInfo,codesFile,2,session);
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Indices of animals belong to a spike group (Include only for epilepsy)
% Indices of animals belong to a spike group
spikeGroups = {'Low Spike Rate';'High Spike Rate'};
spikeThreshold = 2.5;
spikeTable = plotNORspikeTotals(codesFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
    if spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
        spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    elseif spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
        spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    end
end

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Allocations to save memory
coherenceBehavior = cell(1,length(analysisAnimals));
numWindows = zeros(1,length(analysisAnimals));
ctr = 1;
aveCoherence = nan(513,4,length(analysisAnimals));

%% Loop through animals to calculate coherence
tic
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels1{iAnimal}) && ~isempty(channels2{iAnimal})

        % Load FILES
        fileName1 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels1{iAnimal} '.mat']);
        LFP1 = load(fileName1.name);
        
        fileName2 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels2{iAnimal} '.mat']);
        LFP2 = load(fileName2.name);
        
        % Get time/freq info
        sRate = LFP1.sFreq;
        timeStamps = round(LFP1.timestamps - LFP1.timestamps(1),3);
        interestWindowIdx = find(timeStamps >= windowOfInterest(1) & timeStamps < windowOfInterest(2));
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP1,'spike_intervals') % Check to see if this field exists
            LFP1.bad_intervals = [LFP1.bad_intervals; LFP1.spike_intervals];
        end
        
        if isfield(LFP2,'spike_intervals') % Check to see if this field exists
            LFP2.bad_intervals = [LFP2.bad_intervals; LFP2.spike_intervals];
        end
        
        % Filter raw signal for 60 Hz line noise
        d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
        LFP1.values = filtfilt(d_60,LFP1.values);
        LFP2.values = filtfilt(d_60,LFP2.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%% COHERENCE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%
        [interpolatedData1,badIdx1] = interpolateBadData(LFP1); LFP1.values = interpolatedData1;
        [interpolatedData2,badIdx2] = interpolateBadData(LFP2); LFP2.values = interpolatedData2;
        
        %[coherence,freqs] = coherence_AS(LFP1,LFP2,freqRange,timeWindow,stepWindowSize,plotIT)
        [coherence,freqs,timeBins] = coherence_AS(LFP1,LFP2,freqRange,[windowOfInterest; windowOfInterest],4,0);
        timeBins = timeBins + 2; % Step window size = 2 (see above), so the timing of the coherence estimate is shifted ~1 sec.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Frequencies of interest
        freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);
        
        % Find timeBins with artifact and remove
        artifactReject = zeros(1,length(interpolatedData1));
        artifactReject(unique([badIdx1 badIdx2])) = 1;
        artifactReject = logical(artifactReject(interestWindowIdx));
        newTime = 0:1/LFP1.sFreq:(windowOfInterest(2) - windowOfInterest(1))-(1/LFP1.sFreq);
        artifactTimes = newTime(artifactReject);
        for iReshape = 1:(length(timeBins)-1)
           timeBinsReshaped(iReshape,:) = timeBins(iReshape):(1/LFP1.sFreq):(timeBins(iReshape+1) - 1/LFP1.sFreq); 
        end
        removeIdx = find(sum(ismember(timeBinsReshaped,artifactTimes),2));
        coherence(:,removeIdx) = NaN;
        
        % Save calculated information for all animals
        if ctr == 1
            coherenceAll = nan([size(coherence) length(analysisAnimals)]);
            coherenceAll(:,:,iAnimal) = coherence;
            ctr = 2;
        else
            coherenceAll(:,:,iAnimal) = coherence;
        end
        
%         % Sanity Checks: Plot individual raw data and coherence
%         figure('Position',[100 100 1500 1000])
%         colormap(jet)
%         ax(1) = subplot(211);
%         plot(timeStamps(interestWindowIdx),[interpolatedData1(interestWindowIdx)' interpolatedData2(interestWindowIdx)'])
%         legend(electrodes)
%         ax(2) = subplot(212);
%         imagesc(timeBins,1:sum(freqIdx),coherence(freqIdx,:))
%         axis xy
%         set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx))%,'XTick',1:5:length(binSteps),'XTickLabel',binSteps(1:5:end))
%         linkaxes(ax,'x')
        
    else
        continue
    end
    
    if addSpecificBehaviors
        
        [fileNames,timeWindow] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{session},behavior);
        
        % Checks to ensure there are not multiple files
        if isempty(fileNames) || isempty(timeWindow{1})
            continue
        else
       
            numWindows(iAnimal) = size(timeWindow{1},1);
            ctr = 1;
            for iTime = 1:numWindows(iAnimal)          
                startTime = timeWindow{1}(iTime,1) - 1;
                endTime = timeWindow{1}(iTime,2);
                if startTime < min(timeBins) || endTime > max(timeBins)
                    continue
                else
                    timeIdx = timeBins >= startTime & timeBins <= endTime;
                    coherenceBehavior{iAnimal}(:,:,ctr) = coherenceAll(:,timeIdx,iAnimal);
                    ctr = ctr + 1;
                end
            end
            
            aveCoherence(:,:,iAnimal) = nanmean(coherenceBehavior{iAnimal},3);
            
        end
    end
    
    clear analysisChan* fileName* LFP* interpolatedData* timeBinsReshaped timeWindow
    
end
toc

return

%% Plot

% Colors for treatment groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};

freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);

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

%%% Aligned Coherence Matrix - 4 Groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
data = cat(3,nanmean(coherenceAll(:,:,groupIdx{1}),3),nanmean(coherenceAll(:,:,groupIdx{2}),3),nanmean(coherenceAll(:,:,groupIdx{3}),3),nanmean(coherenceAll(:,:,groupIdx{4}),3));
dataCell{1} = data(:,:,2) - data(:,:,1); dataCell{2} = data(:,:,3) - data(:,:,1); dataCell{3} = data(:,:,4) - data(:,:,1);

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = 1:size(data,2);
UDIs.YAXIS = freqs(freqIdx);
UDIs.COLORMAP = 'jet';
UDIs.XLABEL = {'','Time Bins',''};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = {'PILO-SHAM';'STIM-SHAM';'BURST-SHAM'};
UDIs.COLORBARTITLE = '';
UDIs.SAVENAME = [saveDirectory 'coherence_matrix_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
clear data*
%%%

%%% 2 Groups %%%
data = cat(3,nanmean(coherenceAll(:,:,spikeGroupIdx{1}),3),nanmean(coherenceAll(:,:,spikeGroupIdx{2}),3));
dataCell{1} = data(:,:,2) - data(:,:,1);

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = 1:size(data,2);
UDIs.YAXIS = freqs(freqIdx);
UDIs.COLORMAP = 'jet';
UDIs.XLABEL = {'Time Bins'};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = {[spikeGroups{2} ' - ' spikeGroups{1}]};
UDIs.COLORBARTITLE = 'dB';
UDIs.SAVENAME = [saveDirectory 'coherence_matrix_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specific Behaviors
if addSpecificBehaviors
    
    %%% Aligned Coherence Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% 4 Groups %%%
    data = cat(3,nanmean(aveCoherence(freqIdx,:,groupIdx{1}),3),nanmean(aveCoherence(freqIdx,:,groupIdx{2}),3),nanmean(aveCoherence(freqIdx,:,groupIdx{3}),3),nanmean(aveCoherence(freqIdx,:,groupIdx{4}),3));
    dataCell{1} = data(:,:,2) - data(:,:,1); dataCell{2} = data(:,:,3) - data(:,:,1); dataCell{3} = data(:,:,4) - data(:,:,1);

    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.XAXIS = -1:2; %1:size(data,2);
    UDIs.YAXIS = freqs(freqIdx);
    UDIs.COLORMAP = 'jet';
    UDIs.XLABEL = {'','Time (s)',''};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = {'PILO-SHAM';'STIM-SHAM';'BURST-SHAM'};
    UDIs.COLORBARTITLE = 'dB';
    UDIs.SAVENAME = [saveDirectory 'coherence_matrix_4groups'];

    UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
    clear data*
    %%%
    
    %%% 2 Groups %%%
    data = cat(3,nanmean(aveCoherence(:,:,spikeGroupIdx{1}),3),nanmean(aveCoherence(:,:,spikeGroupIdx{2}),3));
    dataCell{1} = data(freqIdx,:,2) - data(freqIdx,:,1);
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.XAXIS = -1:2; %1:size(data,2);
    UDIs.YAXIS = freqs(freqIdx);
    UDIs.COLORMAP = 'jet';
    UDIs.XLABEL = {'Time (s)'};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = {[spikeGroups{2} ' - ' spikeGroups{1}]};
    UDIs.COLORBARTITLE = 'dB';
    UDIs.SAVENAME = [saveDirectory 'coherence_matrix_2groups'];

    UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
    clear data*
    %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Bar Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    
    % Colors for treatment groups
    groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
    
    thetaFreqIdx = freqs >= 6 & freqs <= 10;
    %alphaFreqIdx = freqs > 10 & freqs <= 14;
    
    % Bar group of difference of group averages
    figure('Position',[100 100 1000 500])
    data = squeeze(nanmean(dataDiff(thetaFreqIdx,:,:)));
    h = bar(data);
    set(gca,'XTickLabel',-1:2)
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    %legend(groups)
    xlabel('Time (s)')
    ylabel('Mean Coherence Difference')
    set(gca,'Position',[0.13 0.12 0.8 .4])
	set(h(1),'FaceColor',groupColor{2},'EdgeColor',groupColor{2})
	set(h(2),'FaceColor',groupColor{3},'EdgeColor',groupColor{3})
    set(h(3),'FaceColor',groupColor{4},'EdgeColor',groupColor{4})
    title('Difference from SHAM')
    save2pdf([saveDirectory 'coherence_bar_4groups'])
    clear data dataDiff
    %%%
    
     %%% 2 Groups %%%
    
    % Colors for treatment groups
    spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};
    
    thetaFreqIdx = freqs >= 6 & freqs <= 10;
    %alphaFreqIdx = freqs > 10 & freqs <= 14;
    
    % Bar group of difference of group averages
    figure('Position',[100 100 565 500])
    data = squeeze(nanmean(dataDiff(thetaFreqIdx,:,:)));
    h = bar(data);
    set(gca,'XTickLabel',-1:2)
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    %legend(groups)
    xlabel('Time (s)')
    ylabel('Mean Coherence Difference')
    set(gca,'Position',[0.14 0.12 0.708 .4])
	set(h(1),'FaceColor',spikeGroupColor{1},'EdgeColor',spikeGroupColor{1})
    title('Difference from Low Spike Rate')
    save2pdf([saveDirectory 'coherence_bar_2groups'])
    clear data dataDiff
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

soundAlarm