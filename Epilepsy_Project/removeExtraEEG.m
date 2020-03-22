% These files contain ~1min of baseline that does not correspond to the NOR trial
% ONCE THIS IS RUN, DO NOT RUN AGAIN!

clc; clear; close all;

mainDir = '/Volumes/SP PHD U3/AmberS/Converted Files/';
animalNames = {'A4';'A5';'A6';'A7';'A8'};

% task = 'novel1';
% removeNovel1 = [155 65 94 65 67];

task = 'novel2';
removeNovel2 = [61 69 154 65 75];

for iAnimal = 1:length(animalNames)
    
    loadedDir = rdir([mainDir animalNames{iAnimal} '_' task '*/*.mat']);
    
    for iFile = 1:length(loadedDir)
        
        load(loadedDir(iFile).name,'timestamps','values','sFreq')
        
        timestamps(1:(removeNovel2(iAnimal)*sFreq)) = [];
        values(1:(removeNovel2(iAnimal)*sFreq)) = [];
        
        save(loadedDir(iFile).name,'timestamps','values','-append')
        clear timestamps values sFreq
        
    end
    
    clearvars -except mainDir animalNames task remove* iAnimal
    
end