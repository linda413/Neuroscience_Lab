function [fileNames,timeWindowOfInterest] = readNORbehavior(videoAnalysisFile,animals,task,behavior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readNORbehavior reads in the NOR Video Analysis file

%% UDIs

% Excel file of animal codes
% Ex: videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';

% A24 A25 A26 A28 A29 A30 A31 A33 A34 A35 A36 A37 A39 A40 A41 A42
% Ex: animals = {'A5'};

% Task: novel1 novel2
% Ex: task = {'nove'};

% Specific Behavior
% M (moving)
% NM (not moving)
% S (sitting)
% SNOT (sitting next to object top)
% SNOB (sitting next to object bottom)
% ITO (interacting top object)
% IBO (interacting bottom object)
% ATO (approaching top object)
% ABO (approaching bottom object)
% LTO (leaving top object)
% LBO (leaving bottom object)
% Ex: behavior = {'M';'IBO'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find the timepoints in seconds that correspond to behavior

% If behaviors are not specified, then use all time points in first 2 min.
if isempty(behavior)
    warning('You should probably enter in some specific behaviors.')
    timeWindowOfInterest = [1 120];
end

% Read in data from excel file
[~,txt,~] = xlsread(videoAnalysisFile,1,'C3:DS102');

% Extract animal and task listed in the file
[animalNames,remain] = strtok(txt(:,1),'_');
remain = cellfun(@(x) x(2:end),remain,'UniformOutput',false);
[taskNames,~] = strtok(remain,'_');

% Find the row of the animals of interest
if isempty(animals)
    animalIdx = 1:size(txt,1);
else
    animalIdx = find(ismember(animalNames,animals));
end

% Find the row of the task of interest
if isempty(task)
    taskIdx = 1:size(txt,1);
else
    taskIdx = find(ismember(taskNames,task));
end

% Find the row of both the animal and task of interest
rows = intersect(animalIdx,taskIdx);

% Find the file names
fileNames = txt(rows,1);

% Find the columns of the behaviors of interest
allBehaviorsAndWindows = 1;
findConsecutiveBehaviors = 0;

if isempty(rows)
    timeWindowOfInterest = [];
else
    
    for iRow = 1:length(rows)
        
        if allBehaviorsAndWindows
            
            tempIdx = [];

            for iBehavior = 1:length(behavior)
                
                if ismember(behavior{iBehavior},{'M';'S'})
                    temp = find(strcmp(txt(rows(iRow),2:end),behavior{iBehavior}));
                    tempIdx = [tempIdx temp];
                    clear temp
                else
                    temp = strfind(txt(rows(iRow),2:end),behavior{iBehavior});
                    tempIdx = [tempIdx find(~cellfun(@isempty,temp))];
                end
                
            end
            
            behaviorIdx{iRow} = sort(tempIdx);
            
            % Finds the start and end times of the windows of the behavior
            if isempty(behaviorIdx{iRow})
                startTime = [];
                endTime = [];
            else
                foo = behaviorIdx{iRow};
                startTime = foo([true (diff(foo) ~= 1)])';
                endTime = foo([(diff(foo) ~= 1) true])';
            end
            
            % All time windows
            timeWindows{iRow} = [startTime endTime];
            
            % All time windows at least 2 sec
            timeWindowsOfSize{iRow} = timeWindows{iRow}(endTime - startTime >= 1,:);
            
            % Identify onset of behavior and window around it
            onsets = startTime; % keep this flexible in case want to change to end of a behavior
            beforeWin = 1; afterWin = 1;
            
            if sum(onsets - beforeWin < 1)
                onsets(onsets - beforeWin < 1) = []; % Remove any onsets that are too close to the beginning to window
            end
            
            if sum(onsets + afterWin > 120)
                onsets(onsets + afterWin > 120) = []; % Remove any onsets that are too close to the end to window
            end
            
            timeWindowWindowed{iRow} = [(onsets - beforeWin) (onsets + afterWin)];
             
%             secVec = zeros(1,120); % Find overlapping windows and eliminate
%             for iSec = 1:size(timeWindowWindowed{iRow},1)
%                 secVec(timeWindowWindowed{iRow}(iSec,1):timeWindowWindowed{iRow}(iSec,2)) = 1;
%             end
%             connectComp = bwconncomp(secVec);
%             newOnsets = cellfun(@(x) x(1),connectComp.PixelIdxList);
%             timeWindowWindowed{iRow} = [newOnsets' (newOnsets' + beforeWin + afterWin)];
            
        elseif findConsecutiveBehaviors && length(behavior) == 2
            
            behav1 = find(ismember(txt(rows(iRow),2:end),behavior{1}));
            behav2 = find(ismember(txt(rows(iRow),2:end),behavior{2}));
            
            behav2 = behav2([true (diff( behav2) ~= 1)]);
            
            comboBehav = sort([behav1 behav2]);
            
            startTime = comboBehav(find(diff(sort([behav1 behav2])) == 1));
            endTime = behav2;
            
            timeWindowConsecutiveBehaviors{iRow} = [startTime' endTime'];
            
        end
        
        timeWindowOfInterest(iRow) = timeWindowWindowed;
        
        clear foo startTime endTime
        
    end
    
end

end