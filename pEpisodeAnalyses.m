%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P-Episode analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input into pEpisode)
%directoryPath = '/Volumes/Gurkoff Lab/Amber/Converted Files/NOR/';
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'novel1';

% Group to analyze
groupInput = {'SHAM';'PILO'};

% Treatment to analyze
treatmentInput = {'NO'; 'STIM';'BURST'};

% Region to analyze (enter only 1)
electrode = {'LV'};

% Enter a vector of frequencies (log-space) over which you would like to
% calculate pEpisode (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
backgroundFrequencies = logspace(log10(2),log10(30),30);

% Enter the frequency band of particular interest to be output to the .csv and figure files
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that will be used for the background spectrum
% Example: time window 0 min - 5 min -> backgroundWindow = [0 300];
backgroundWindow = [0 300];

% Enter the time window in seconds that you are curious to examine pEpisode
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrode{1} '/'];

% Calculate Pepidose for specific behaviors
addSpecificBehaviors = 1;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate Pepisode

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,1,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session);

% Index of range of frequencies of interest
[~,idxF1] = min(abs(backgroundFrequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(backgroundFrequencies - frequenciesOfInterest(2)));

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Allocations
PepisodeAll = nan(size(fileInfo,1),length(backgroundFrequencies));
freqTimeWinAll = nan(length(backgroundFrequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
powerAll = nan(length(backgroundFrequencies),length(backgroundWindow(1):backgroundWindow(2)-1)*1000,size(fileInfo,1));
behaviorPepisode = nan(size(fileInfo,1),1);
behaviorPepisodeAlign = nan(length(backgroundFrequencies),3000,size(fileInfo,1));
[numWindows,percentArtifactBackground,percentArtifactInterestWindow] = deal(nan(size(fileInfo,1),1));

%% Loop through animals and removes bad channels and computes p-episode
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels{iAnimal})
        
        % Load file
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = round(LFP.timestamps - LFP.timestamps(1),3); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps >= windowOfInterest(1) & timeStamps < windowOfInterest(2));
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP,'spike_intervals') % Check to see if this field exists
            LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% PEPISODE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%
        pEpisodeDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Pepisode/';
        pEpisodeFile = [pEpisodeDir 'Pepisode_' analysisAnimals{iAnimal} '_' session '_' channels{iAnimal} '.mat'];
        if exist(pEpisodeFile,'file') == 2
            load(pEpisodeFile)
        else
            % Filter raw signal for 60 Hz line noise
            d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
            LFP.values = filtfilt(d_60,LFP.values);
            %fvtool(d_60) %Look at filter reponse
            
            [Pepisode,freqTimeWin,power] = pEpisode_AS(LFP,backgroundFrequencies,frequenciesOfInterest,backgroundWindow,windowOfInterest,[]);
            
            if ~exist(pEpisodeDir,'dir')
                mkdir(pEpisodeDir)
            end
            save(pEpisodeFile(1:end-4),'backgroundFrequencies','frequenciesOfInterest','backgroundWindow','windowOfInterest','Pepisode','freqTimeWin','power')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Finds portion of signal that is artifact/spike
        artifactBackground = isnan(power(1,:));
        percentArtifactBackground(iAnimal) = mean(artifactBackground);
        percentArtifactInterestWindow(iAnimal) = mean(artifactBackground(interestWindowIdx));
        
        % Save calculated information for all animals
        PepisodeAll(iAnimal,:) = Pepisode;
        freqTimeWinAll(:,:,iAnimal) = freqTimeWin;
        powerAll(:,:,iAnimal) = power;
        
        %         % Sanity Checks: Plot individual raw data, power, & Pepisode
        %         figure('Position',[100 100 1500 1000])
        %         ax(1) = subplot(311);
        %         backWindowIdx = timeStamps >= backgroundWindow(1) & timeStamps < backgroundWindow(2);
        %         rawData = LFP.values(backWindowIdx); rawData(artifactBackground) = NaN;
        %         plot(timeStamps(interestWindowIdx),rawData(interestWindowIdx))
        %         title(saveName)
        %         ax(2) = subplot(312);
        %         imagesc(timeStamps(interestWindowIdx),1:length(backgroundFrequencies),squeeze(power(:,interestWindowIdx)))
        %         set(gca,'YTick',1:5:length(backgroundFrequencies),'YTickLabel',backgroundFrequencies(1:5:end))
        %         title('Power')
        %         ax(3) = subplot(313);
        %         imagesc(timeStamps(interestWindowIdx),1:length(backgroundFrequencies),squeeze(freqTimeWin))
        %         set(gca,'YTick',1:5:length(backgroundFrequencies),'YTickLabel',backgroundFrequencies(1:5:end))
        %         title('Pepisode')
        %         linkaxes(ax,'x')
        %         %saveName = [analysisAnimals{iAnimal} ' ' analysisChan ' Raw Power Pepisode'];
        %         %saveas(gcf,[saveDirectory saveName])
        
    else
        continue
    end
    
    % Calculates Pepisode for the time windows that correspond to the
    % behaviors of interest
    if addSpecificBehaviors
        
        [fileNames,timeWindow] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{session},behavior);
        
        % Checks to ensure there are not multiple files
        if isempty(fileNames) || isempty(timeWindow{1})
            continue
        else
            
            timeVec = zeros(1,length(interestWindowIdx));
            
            numWindows(iAnimal) = size(timeWindow{1},1);
            
            for iTime = 1:numWindows(iAnimal)
                
                timeVec(timeWindow{1}(iTime,1)*sRate - 999:timeWindow{1}(iTime,2)*sRate) = 1;
                timeVec = logical(timeVec);
                
                startTime = timeWindow{1}(iTime,1)*sRate - 999;
                endTime = timeWindow{1}(iTime,2)*sRate;
                freqTimeWinAlign(:,:,iTime) = freqTimeWin(:,startTime:endTime);
                
            end
            
            % Averages across all time points (multiple interacting windows)
            behaviorPepisode(iAnimal) = nanmean(nanmean(freqTimeWin(idxF1:idxF2,timeVec),2),1)*100;
            
            % Adds across all interacting windows
            behaviorPepisodeAlign(:,:,iAnimal) = nansum(freqTimeWinAlign,3);
            
        end
        
    end
    
    clear filename LFP pEpisodeFile fileNames timeWindow freqTimeWinAlign
    
end

soundAlarm