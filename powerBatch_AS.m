%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Batch Power for Gurkoff Lab 

% Amber Schedlbauer
% 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input into power)
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';

% Name of directory you would specifically like to analyze
% Example: dataDirectory = 'A4_novel1_070116';
dataDirectory = 'A4_novel1_070116';

% Enter a vector of frequencies (log-space) over which you would like to
% calculate power (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
backgroundFrequencies = logspace(log10(2),log10(54),30);

% Enter the time window in seconds that you are curious to examine power
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Enter the time window in seconds that will be used for baseline normalization
% Example: time window 0 min - 5 min -> backgroundWindow = [0 300];
backgroundWindow = [0 300];

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' dataDirectory '/'];

%% Power Batch 

% Grabs the files in the data directory
allFiles = rdir([directoryPath dataDirectory '/*.mat']); % FYI rdir is a non-Matlab function

for iFile = 1:length(allFiles)
    
    LFP = load(allFiles(iFile).name);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% INSERT SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert spike information into LFP structure as bad_intervals
    if isfield(LFP,'spike_intervals') % Check to see if this field exists
        LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
    else
        continue % If the field does not exist, then the channel was proabably determined to be artifact.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract names of files
    [~,name,~] = fileparts(allFiles(iFile).name);
    name = strrep(name,' ','_'); %Gets rid of any spaces and replaces with underscore
    colNames{iFile} = strrep(name,'-','_'); %Gets rid of any dashes and replaces with underscore
    
    % Calculate Power
    %OUT = power_AS(LFP,freqs,window,baselineWindow,PLOT)
    OUT(iFile) = power_AS(LFP,backgroundFrequencies,windowOfInterest,backgroundWindow,1);
    
    % Label figure name
    set(gcf,'Name',colNames{iFile})
   
    clear LFP name
    
end

% Save outputs
% User enter code here to save matrix of output or figures

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir') 
    mkdir(saveDirectory)
end



