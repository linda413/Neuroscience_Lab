                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Batch Count Interictal Spikes for Gurkoff Lab
%
% Amber Schedlbauer
% 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% UDIs

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/shared/Katelynn/20495/20495_Channels.xlsx';

% Enter the path of where all the data are located
mainDir = '/shared/Katelynn/20495/Converted Files/Behavior/NL/';

% Enter the recording session you would like to count
session = 'NL2';

% Group to analyze
groupInput = {'Pilo';'Sham'};

% Treatment to analyze
treatmentInput = {'None'; 'Stim'}; %

% Region to analyze (enter only 1)
electrode = {'MSN'};

% Enter the bad channel file/sheet number
badChanFile = '/shared/Katelynn/20495/20495_Channels.xlsx';
sheetNum = 2;

% Enter directory where you want to save your data
saveDir = '/shared/Katelynn/20495/Results/Spikes/Post-Stim/NL/';                  

%% Determine Animals and bad channels

% Identifies which animals to anlayze
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,0,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
animals = fileInfo(:,1);

% Identifies which channels for the identified animals
channels = readBadChannels(fileInfo,badChanFile,sheetNum,session);

%% Caculate spike totals and spike rates for all animals 
[spikeTotals,spikeRates] = countInterictalSpikes_AS(animals,mainDir,session,channels,saveDir);