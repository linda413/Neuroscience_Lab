%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code to convert files containing 3 or 4 animals to indiviudal files
% NOTE: Heavy modification is needed to indivualize for your desired file
% names!!! Once you are SURE you have appropriately named things, uncomment
% Line 32.

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctr = 1;

for i = 1:16 %17:32 33:48 49:64
    
    % Directory that contains the CSC files for multiple animals
    mainDir = '/shared/Greg/Converted Files 2/A58_57_54_day3/';
    
    % Grabs the files to change the names (i.e. CSC9 to CSC1)
    file2change = dir([mainDir '*_CSC' num2str(i) '.mat']);
    
    [~,f] = fileparts(file2change(1).name);
    
    % New name of folder for indivual animal
    newName1 = 'A54_day3';
    mkdir([mainDir newName1])
    newName2 = [f(15:30) num2str(ctr) '.mat'];
    
    disp(['Old file name: ' f])
    disp([newName1 newName2])
    
    %movefile([mainDir f '.mat'],[mainDir newName1 '/' newName1 newName2])
    
    ctr = ctr + 1;
    
end