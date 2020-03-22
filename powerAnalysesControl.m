%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Power Analyses for Gurkoff Lab
%
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
directoryPath = '/Volumes/Gurkoff Lab/Amber/Converted Files/Epileptogenesis/';
%directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'postpilostim04';

% Group to analyze
groupInput = {'SHAM'}; %;'PILO'

% Treatment to analyze
treatmentInput = {'NO'}; %; 'STIM';'BURST'

% Region to analyze (enter only 1)
electrode = {'LV'};

% Enter a vector of frequencies (log-space) over which you would like to
% calculate power (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
frequencies = logspace(log10(2),log10(30),30);

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Enter a "1" or "0" if you want to use baseline normalization
% If so, modify the baselineNormalization function to suit your needs.
useBaselineNormalization = 1;

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrode{1} '/'];

%% Identify files to analyze and calculate power

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
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session);

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(frequencies - frequenciesOfInterest(2)));

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Allocations
powerAll = nan(length(frequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
powerAllConvert = nan(length(frequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));

%% Loop through animals and removes bad channels and computes power
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels{iAnimal})
        
        % Load file
        %fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels{iAnimal} '.mat']);
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = LFP.timestamps - LFP.timestamps(1); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps >= windowOfInterest(1) & timeStamps < windowOfInterest(2)); % Index of window of interest
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP,'spike_intervals') % Check to see if this field exists
            LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
        end
        
        % Filter raw signal for 60 Hz line noise
        d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
        LFP.values = filtfilt(d_60,LFP.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%%%% POWER CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
        OUT(iAnimal) = power_AS(LFP,frequencies,windowOfInterest,0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save calculated information for all animals
        powerAll(:,:,iAnimal) = OUT(iAnimal).power(:,interestWindowIdx);

        % Conduct baseline normalization on raw power values
        % [powerDB,powerZ] = baselineNormalization(powerInfo,analysisAnimal,task)
        if useBaselineNormalization
            [powerDB,powerZ] = baselineNormalization(OUT(iAnimal).power,OUT(iAnimal).power(:,interestWindowIdx),analysisAnimals{iAnimal},session);
            powerAllConvert(:,:,iAnimal) = powerDB;
        end
        
    else
        continue
    end
    
    clear fileName LFP fileNames timeWindow tempPower
    
end

%% Insert NaNs and calculate the PSDs
PSD = squeeze(nanmean(powerAllConvert,2));
ctr = 2;

for iInsert = 50:50:600 % 60 = max number of spikes in two minutes for this study
    
    randIdx = sort(randi(size(powerAllConvert,2),iInsert,1));
    randIdx = [randIdx randIdx+350]; % 350 = size of spike window
    
    powerAllTemp = powerAllConvert;
    
    for iSpike = 1:size(randIdx,1)
        powerAllTemp(:,randIdx(iSpike,1):randIdx(iSpike,2),:) = NaN;
    end
    
    numberNAN(:,ctr-1) = sum(isnan(powerAllTemp(1,:,:)));
    PSD(:,:,ctr) = squeeze(nanmean(powerAllTemp,2));
    ctr = ctr + 1;
    
end

[h,p,~,~] = ttest2(PSD(:,:,1)',PSD(:,:,end)')

plot(squeeze(nanmean(PSD(:,:,[1 2]),2)),'LineWidth',2)

clearvars -except numberNAN