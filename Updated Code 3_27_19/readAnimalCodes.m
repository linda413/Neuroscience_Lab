function fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,behaviorInput,saveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readNORcodes reads in the animal codes from an excel file

% Function Inputs
%
% Excel file containing the animal informaiton
% Ex: excelFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
% Ex: excelFile = '/home/amschedl/Desktop/NOR/Behavior/NOR Animal Codes.xlsx';
% 
% Choose which animals you would like to analyze
% You are allowed to leave empty: animalsInput = [];
% Ex: animalsInput = {'A24' 'A25' 'A26'};
% 
% Choose which group you would like to analyze
% Ex: groupInput = {'SHAM';'PILO'};
% 
% Choose which treatment you would like to analyze
% Ex: treatmentInput = {'NO'; 'STIM';'BURST'};
% 
% Choose which electrodes you would like to analyze
% Ex: electrodeInput = {'LPFC'; 'RPFC'; 'LMSN'; 'RMSN'; 'LV'; 'RV'; 'LD'; 'RD'}; 
% 
% Choose if you want analyze animals that completed the NOR task
% Ex: NORinput = 1;
% 
% Directory to save information
% Ex: saveDir = '/home/amschedl/Desktop/NOR/Code/Results/';

% Amber Schedlbauer
% 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find relvant info like animal number (if not input) and channel names of interest

% Read in excel file
[~,txt,~] = xlsread(excelFile,1);

% Extract animal indices in the file
if isempty(animalsInput)
    % Looks for strings containing the letter 'A'
    animalIdx = strfind(txt(2:end,1),txt{2,1}(isletter(txt{2,1}))); % 2nd argument of strfind finds the characters in the string
    animalIdx = find(~cellfun(@isempty,animalIdx))+1;
else
    animalIdx = find(ismember(txt(:,1),animalsInput));
end

% Extract group indices in the file
if isempty(groupInput)
    error('Please input the groups you would like to analyze.')
else
    groupIdx = find(ismember(txt(:,2),groupInput));
end

% Extract treatment indices in the file
if isempty(treatmentInput)
    error('Please input the treatments you would like to analyze.')
else
    treatmentIdx = find(ismember(txt(:,3),treatmentInput));
end

% Identify rows of the animals of interest in excel sheet
rows = intersect(animalIdx,groupIdx);
rows = intersect(rows,treatmentIdx);

% Extract electrode indices in the file
if isempty(electrodeInput)
    error('Please input the electrodes you would like to analyze.')
else
    electrodeIdx = find(ismember(txt(1,:),electrodeInput));
end

% Selected animals based upon the group, treatment, and electrodes
currentSelection = txt(rows,[1:3 electrodeIdx]);

% Removes animals that didn't complete NOR
if behaviorInput
    
    [~,txtNOR,~] = xlsread(excelFile,4);
    NORidx = strcmpi(txtNOR(:,4),'yes');
    animalNOR = txtNOR(NORidx,1);
    
    [~,~,keepIdx] = intersect(animalNOR,currentSelection(:,1),'stable');
    currentSelection = currentSelection(keepIdx,:);
    
end

% Identify stimulating electrodes and rename electrodes
justElectrodes = currentSelection(:,4:end);

for iReplace = 1:numel(justElectrodes)
    
    % Identifies if cell has a stimulating channel
    goodChan = strfind(justElectrodes{iReplace},'AD');
    if isempty(goodChan)
        justElectrodes{iReplace} = '';
    else
        elecNum = regexp(justElectrodes{iReplace},'\d*','Match');
        newElecNum = str2double(elecNum) + 1;
        justElectrodes{iReplace} = ['CSC' num2str(newElecNum)];
    end

end

justElectrodes = reshape(justElectrodes,size(currentSelection,1),[]);
currentSelection(:,4:end) = justElectrodes;

% Save information as a new variable
fileInfo = currentSelection;

% Find file names of files of interest and saves to variable
% Add outputFileName variable to input arguments (needs directory and task fields
% if exist(outputFileName,'var')
%     ctr = 1;
%     for iAnimal = 1:size(fileInfo,1)
%         for iChan = 1:size(fileInfo(:,4:end),2)
%             if ~isempty(fileInfo{iAnimal,iChan})
%                 temp = rdir([outputFileName.directory fileInfo{iAnimal,1} '_' outputFileName.task '_*/*' fileInfo{iAnimal,iChan} '.mat']);
%                 fileNames{ctr} = temp.name;
%                 ctr = ctr + 1;
%             end
%         end
%     end    
% end

% Save to specified directory
if ~isempty(saveDir)
    
    % Makes saveDirectory if it doesn't exist
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    save([saveDir 'files2analyze_' savedate],'fileInfo')
    
end
