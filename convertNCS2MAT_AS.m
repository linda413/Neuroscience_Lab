function convertNCS2MAT_AS(inputDir,outputDir,filesToAnalyze,sRate,suffix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convertNCS2MAT_AS(inputDir,outputDir,filesToAnalyze,trialNum,sRate,suffix)
%
% Function to convert Neuralynx files (.ncs) into Matlab (.mat) files
% If you want to modify any of the fields of the LFP structure, this would
% be the place to do it!
%
% Used code from Matt in UoA
%
% NOTE: Need to still address the multiple trials within a recording
% session.
%
% INPUTS:
%   inputDir: directory where the files you want to convert are located
%   outputDir: directory you want to save the converted files
%   filesToAnalyze: convert only a subset of the CSC files
%   sRate: sampling rate in Hz
%   suffix: string that is added on to the end of saved files
%
% Amber Schedlbauer - 2018 and Matthew 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get .ncs file names in input directory
cd(inputDir)
filesInDir = dir('CSC*.ncs');
filesInDir = {filesInDir.name};

% Extract only the .ncs files in the channels list
filesInDir = intersect(filesToAnalyze,filesInDir);   

% Checks to ensure there are not multiple files per channel
% if sum(contains(filesInDir,'_')) ~= 0
%     error('Uh oh! This directory has multiple files per channel.')
% end

% Loop through each .ncs file and convert to .mat file
for iFile = 1:length(filesInDir)
    
    % Get directory name without full path
    idx = strfind(inputDir,'/');
    outputLabel = inputDir(idx(end-1)+1:end-1);
    
    % Convert the data using the code provided by Neuralynx
    %[timeStamps,~,~,~,samples,headerInfo] = Nlx2MatCSC(filesInDir{iFile},[1 1 1 1 1],1,1,[]);
    [timeStamps,chanNum,sampleFrequency,numValSamples,samples,headerInfo] = Nlx2MatCSC(filesInDir{iFile},[1 1 1 1 1],1,1,[]);
    
    keyboard
    
    % Matt's code to handle header info
    version = (char(regexp(headerInfo{7},'(?<=\s)\S*','match')));
    
    switch version
        case '3.3.0'
            
            LFP.old_sFreq = str2double(char(regexp(headerInfo{14},'(?<=\s)\S*','match')));
            LFP.bits_to_volt_conversion = str2double(char(regexp(headerInfo{16},'(?<=\s)\S*','match')));
            
            trialStart = regexp(headerInfo{3},'(\d+):(\d+):(\d+)','tokens');
            LFP.Trial_Start = sum(str2double(trialStart{1,1}).*[3600,60,1]);
            
            trialEnd = regexp(headerInfo{4},'(\d+):(\d+):(\d+)','tokens');
            LFP.Trial_End = sum(str2double(trialEnd{1,1}).*[3600,60,1]);
            
        otherwise % If the version is 3.3 or older the timestamps move.
            
            LFP.old_sFreq = str2double(char(regexp(headerInfo{15},'(?<=\s)\S*','match')));
            LFP.bits_to_volt_conversion = str2double(char(regexp(headerInfo{17},'(?<=\s)\S*','match')));
            
            trialStart = regexp(headerInfo{8},'(\d+):(\d+):(\d+)','tokens');
            if isempty(trialStart)
                trialStart = regexp(headerInfo{3},'(\d+):(\d+):(\d+)','tokens');
            end
            LFP.Trial_Start = sum(str2double(trialStart{1,1}).*[3600,60,1]);
            
            trialEnd = regexp(headerInfo{9},'(\d+):(\d+):(\d+)','tokens');
            if isempty(trialEnd)
                trialEnd = regexp(headerInfo{4},'(\d+):(\d+):(\d+)','tokens');
            end
            LFP.Trial_End = sum(str2double(trialEnd{1,1}).*[3600,60,1]);
    end
    
    % Reshape the data into a single dimension vector, conver to mV
    values = reshape(samples,1,[]) * LFP.bits_to_volt_conversion * 1000;
    
    % The Timestamps are only recorded every 512 points, but if we want to
    % make analysis easier later, we need each point to be assigned a
    % timestamp.
    for tix = 2:length(timeStamps)
        holding = linspace(timeStamps(tix-1),timeStamps(tix),size(samples,1)+1)';
        TS(:,tix-1) = holding(1:end-1);
    end
  
    TS = reshape(TS,1,[]);
    holding = timeStamps(end):mean(diff(TS)):(timeStamps(end)+mean(diff(TS)) * size(samples,1));
    TS = [TS holding(1:end-1)];
    
    % Downsample the data -> resample produces ring 6/17/2016
    % Because decimate was also using a filter, it was causing the
    % timestamps to vary in delta-t. Just like spike2 we will now literally
    % only use one in 32 points, which is the usual ratio, for the final
    % data.
    decimationFactor = LFP.old_sFreq/sRate; % ratio of sampling frequencies
    LFP.timestamps = TS(1:decimationFactor:end)/10^6;
    LFP.values  = values(1:decimationFactor:end);
    LFP.units = {'mV','s'};
    LFP.sFreq = LFP.old_sFreq/decimationFactor;
    
    % Input for LFP structure fields
    LFP.Channel = str2double(filesInDir{iFile}(4:end-4));
    LFP.Trial = 1;
    LFP.Date = date; 
    LFP.Description = [outputLabel '_' filesInDir{iFile}(1:end-4)];

    if isempty(suffix)
        name = [outputLabel '_' date '_' filesInDir{iFile}(1:end-4)];
    else
        name = [outputLabel '_' date '_' filesInDir{iFile}(1:end-4) '_' suffix];
    end
    
    LFP.name = name; %store the filename of the channel
    
    LFP.bad_intervals = [1 2]; % Inelegant, but lets the artifact rejection run smoothly.
    
    display(['File ' num2str(iFile) '/' num2str(length(filesInDir)) ' Complete'])

    if ~exist([outputDir outputLabel '/'],'dir')
        mkdir([outputDir outputLabel '/'])
    end

    save([outputDir outputLabel '/' name],'-struct','LFP')
    
    clearvars -except inputDir outputDir filesToAnalyze sRate suffix filesInDir
    
end