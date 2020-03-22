function [performanceNOR,performanceBarnes]  = readPerformance(performanceFile,sheetNumNOR,colNumNOR,sheetNumBarnes,colNumBarnes,animals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Behavioral Performance

% Inputs
%   performanceFile  - name of excel file containing calculated performance
%   sheetNum - sheet number of excel file where info is located
%   animals - cell array of animal names for which you want their performance
%
% Output
%   performanceNOR - vector of NOR performance values
%   performanceBarnes - vector of Barnes performance values
%
% Amber Schedlbauer
% 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NOR
 [~,~,raw] = xlsread(performanceFile,sheetNumNOR);
 
 temp = cellfun(@(x) find(strcmp(x,raw(:,1))),animals,'UniformOutput',false);
 performanceNOR = cell2mat(raw(cell2mat(temp),colNumNOR));
 
 %% Barnes
 [~,~,raw] = xlsread(performanceFile,sheetNumBarnes);
 
 temp = cellfun(@(x) find(strcmp(x,raw(:,1))),animals,'UniformOutput',false);
 performanceBarnes = cell2mat(raw(cell2mat(temp),colNumBarnes));
 performanceBarnes = sum(performanceBarnes,2);
 
end