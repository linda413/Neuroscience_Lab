%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coherence Analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
directoryPath = '/shared/Katelynn/20495/Converted Files/Epileptogenesis/';

% Animal information excel file
codesFile = '/shared/Katelynn/20495/20495_Channels.xlsx';

% Task to analyze (enter only 1)
session = 'epileptogenesis3';

% Group to analyze
groupInput = {'Sham'};

% Treatment to analyze
treatmentInput = {'None'};

% Regions to analyze (enter only 2)
electrodes = {'vHPC2','dPFC'};

% Enter the frequency range for coherence anlaysis
freqRange = [2 30];

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [10 15];

% Path for folder where the results will be saved

saveDirectory = ['/shared/Katelynn/20495/Results/Coherence/' savedate '/' session '/' electrodes{1} '/'];

%% Identify files to analyze and calculate coherence

% Grabs the animals of interest
%fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{1},0,[]);
channels1 = readBadChannels(fileInfo,codesFile,2,session);
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{2},0,[]);
channels2 = readBadChannels(fileInfo,codesFile,2,session);
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
%groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

%% Loop through animals to calculate coherence

ctr = 1;

for iAnimal = 1:length(analysisAnimals)
    
    if ~isempty(channels1{iAnimal}) && ~isempty(channels2{iAnimal})
        % Load FILES
        fileName1 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels1{iAnimal} '.mat']);
        LFP1 = load(fileName1.name);
        
        fileName2 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels2{iAnimal} '.mat']);
        LFP2 = load(fileName2.name);
        
%         % Load FILES with if loops added 12/19 by KO
%         fileName1 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '*/' analysisAnimals{iAnimal} '_' session '*' channels1{iAnimal} '.mat']);
%         
%         if exist(fileName1.name,'file') == 2
%             LFP1 = load(fileName1.name); disp(fileName1.name);
%         else
%             continue
%         end
%         
%         fileName2 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '*/' analysisAnimals{iAnimal} '_' session '*' channels2{iAnimal} '.mat']);
%         if exist(fileName2.name,'file') == 2
%             LFP2 = load(fileName2.name); disp(fileName2.name);
%         else
%             continue
%         end
%         
        % Get time/freq info
        sRate = LFP1.sFreq;
        timeStamps = round(LFP1.timestamps - LFP1.timestamps(1),3);
        interestWindowIdx = find(timeStamps == windowOfInterest(1)):(find(timeStamps == windowOfInterest(2))-1);
        
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
        timeBins = timeBins + 2; % Step window size = 4 (see above), so the timing of the coherence estimate is shifted
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Frequencies of interest
        freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);
        freqThetaIdx = freqs >= frequenciesOfInterest(1) & freqs <= frequenciesOfInterest(2);
        
        % Find timeBins with artifact and remove
        artifactReject = zeros(1,length(interpolatedData1));
        artifactReject(unique([badIdx1 badIdx2])) = 1;
        artifactReject = logical(artifactReject(interestWindowIdx));
        newTime = 0:1/LFP1.sFreq:(windowOfInterest(2) - windowOfInterest(1))-(1/LFP1.sFreq);
        artifactTimes = round(newTime(artifactReject),3);
        for iReshape = 1:length(timeBins)
            timeBinsReshaped(iReshape,:) = timeBins(iReshape):(1/LFP1.sFreq):(timeBins(iReshape) + 1 - 1/LFP1.sFreq);
        end
        timeBinsReshaped = round(timeBinsReshaped,3);
        removeIdx = find(sum(ismember(timeBinsReshaped,artifactTimes),2));
        coherence(:,removeIdx) = NaN;
        
        % Save calculated information for all animals
        if ctr == 1
            coherenceAll = nan([size(coherence) length(analysisAnimals)]);
            behaviorMatAll = nan([sum(freqIdx) size(coherence,2) length(analysisAnimals)]);
            ctr = 2;
        end
        coherenceAll(:,:,iAnimal) = coherence;
        timeSqueezeCoherence = squeeze(nanmean(coherenceAll, 2));
        animalSqueezeCoherence = squeeze(nanmean(coherenceAll, 3));
        
        % Sanity Checks: Plot individual raw data and coherence
%         figure('Position',[100 100 1500 1000])
%         colormap(jet)
%         ax(1) = subplot(211);
%         plot(timeStamps(interestWindowIdx),[interpolatedData1(interestWindowIdx)' interpolatedData2(interestWindowIdx)'])
%         set(gca,'XTick',timeBins,'XTickLabel',timeBins)
%         legend(electrodes)
%         title(analysisAnimals{iAnimal})
%         ax(2) = subplot(212);
%         imagesc(timeBins+0.50,1:sum(freqIdx),coherence(freqIdx,:)) % Added the 0.50 so that the coherence for the first bin would be plotted at the tick mark and not have the tick mark in the middle of the bin
%         axis xy
%         set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx),'XTick',timeBins,'XTickLabel',newLabels)
%         linkaxes(ax,'x')
%         waitfor(gcf)
        
    else
        continue
    end

    clear fileName* LFP* interpolatedData* coherence timeBinsReshaped 
    
end

%%% Spectrograms %%%
% freqIdx = logical(freqs >= freqRange(1) & freqs <= freqRange(2));
% 
% clear dataCell
% UDIs.groups = groups([1 2 3]);
% UDIs.groupColor = {};
% UDIs.XAXIS = timeBins;
% UDIs.YAXIS = freqs(freqIdx);
% UDIs.COLORMAP = 'jet';
% UDIs.XLABEL = {'Time Bins', 'Time Bins'};
% UDIs.YLABEL = 'Frequency (Hz)';
% UDIs.TITLE = {'SHAM', 'PILO', 'STIM'};
% UDIs.COLORBARTITLE = 'Coherence';
% UDIs.SAVENAME = [saveDirectory 'coherence_' session '25s_' electrodes{1} '_' electrodes{2} '_sham-vs-TBI-NO'];
% dataCell{1}= squeeze(nanmean(coherenceAll(freqIdx, :, groupIdx{1}), 3));
% dataCell{2}= squeeze(nanmean(coherenceAll(freqIdx, :, groupIdx{2}), 3));
% dataCell{3}= squeeze(nanmean(coherenceAll(freqIdx, :, groupIdx{3}), 3));
% UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
% 
% 
