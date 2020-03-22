%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phase Analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
%directoryPath = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Converted Files/';
directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
task = 'novel1';

% Region to analyze (enter only 1)
electrode = {'LV'};

% Enter a vector of frequencies (log-space) over which you would like to
% calculate power (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
frequencies = logspace(log10(2),log10(30),50);

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/'];

% Calculate power for specific behaviors
addSpecificBehaviors = 1;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate power

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],{'SHAM';'PILO'},{'NO';'STIM';'BURST'},electrode,1,[]);
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
spikeTable = plotNORspikeTotals(codesFile,'/Volumes/SP PHD U3/AmberS/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
    if spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
        spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    elseif spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
        spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    end
end

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,task);

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(frequencies - frequenciesOfInterest(2)));

% Allocations
phaseAll = nan(length(frequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
numWindows = nan(size(fileInfo,1),1);
phaseAlign = cell(1,size(fileInfo,1));
%powerConvertAlign = nan(length(frequencies),3000,size(fileInfo,1));

%% Loop through animals and removes bad channels and computes power
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels{iAnimal})
        
        % Load file
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' task '_*/' analysisAnimals{iAnimal} '_' task '*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = LFP.timestamps - LFP.timestamps(1); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps > windowOfInterest(1) & timeStamps <= windowOfInterest(2)); % Index of window of interest
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP,'spike_intervals') % Check to see if this field exists
            LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
        end
        
        % Filter raw signal for 60 Hz line noise
        d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
        LFP.values = filtfilt(d_60,LFP.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%%%% POWER CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
        OUT(iAnimal) = peakFreqPower_AS(LFP,frequencies,windowOfInterest,[0 300],0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save calculated information for all animals
        phaseAll(:,:,iAnimal) = OUT(iAnimal).phase(:,interestWindowIdx);
        
    else
        continue
    end

    % Calculates power for the time windows that correspond to the
    % behaviors of interest
    if addSpecificBehaviors
        
        [fileNames,timeWindow] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{task},behavior);
        
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
                phaseAlign{iAnimal}(:,:,iTime) = phaseAll(:,startTime:endTime,iAnimal);
                
            end
            
            %phaseAlign(:,:,iAnimal) = nanmean(tempPower,3);
            
        end
    end
    
    clear fileName LFP fileNames timeWindow tempPower
    
end

%% Plot

ITPC = cellfun(@(x) abs(nanmean(exp(1i*x),3)),phaseAlign,'UniformOutput',false);
for iGroup = 1:4
    data = ITPC{groupIdx{iGroup}'};
    data = reshape(data,50,3000,[]);
    figure
    imagesc(nanmean(data,3))
end


for iWin = 1:3
    
    windows = [1 1000; 1001 2000; 2001 3000];
   
    data = phaseAlign{3}(idxF1:idxF2,windows(iWin,1):windows(iWin,2),:);
    figure
    polarhistogram(data(:),100)
    
    imagesc(ITPC{3})
    
end