function statusInfo = readStatus(excelFile,sheetNum,animalsInput)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readNORstatus reads in the number of cycles of status

% INPUTS
%
% Excel file containing the animal informaiton
% Ex: excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Status.xlsx';
%
% sheetNum: number of excel sheet
% Ex: sheetNum = 2;
% 
% Choose which animals you would like to analyze
% You are allowed to leave empty: animalsInput = [];
% Ex: animalsInput = {'A24' 'A25' 'A26'};

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find relvant info like animal number (if not input)

% Read in excel file
[~,~,raw] = xlsread(excelFile,sheetNum);
raw = raw(2:end,:);

% Extract animal indices in the file
if isempty(animalsInput)
    % Looks for strings containing the prefix for animal naming scheme
    animalIdx = strfind(raw(:,1),raw{1,1}(isletter(raw{1,1}))); 
    animalIdx = find(~cellfun(@isempty,animalIdx))+1;
else
    animalIdx = find(ismember(raw(:,1),animalsInput));
end

statusInfo = raw(animalIdx,:);

end