function channels = readBadChannels(fileInfo,badChanFile,sheetNum,task)

[~,txt,raw] = xlsread(badChanFile,sheetNum); % first sheet
channels = fileInfo(:,4);

if strfind(txt{1,2},'CSC') == 1
    
    for iAnimal = 1:size(fileInfo,1)
        
        animalIdx = strcmp(txt(:,1),fileInfo{iAnimal,1}); % row number of animal of interest in NOR Bad Channels.xlsx
        colIdx = ismember(txt(1,:),channels{iAnimal});
        if raw{animalIdx,colIdx} == 1
            channels(iAnimal) = cell(1);
        end
        
    end
    
else
    
    currentCol = strcmp(txt(1,:),task);
    
    animalIdx = ismember(txt(:,1),fileInfo(:,1));
    bad = raw(animalIdx,currentCol);
    badSession = zeros(length(bad),1);
    
    for iBad = 1:length(bad)
        
        if isempty(channels{iBad})
            badSession(iBad) = 1;
        else
            if ~isnan(bad{iBad})
                if strcmp(bad(iBad),'all')
                    badSession(iBad) = 1;
                else
                    temp1 = regexp(channels{iBad},'\d+','match'); % Number of channel
                    if isnumeric(bad{iBad})
                        bad{iBad} = num2str(bad{iBad});
                    end
                    temp2 = regexp(bad{iBad},'\d+','match'); % Numbers in excel sheet
                    badSession(iBad) = ismember(num2str(temp1{1}),temp2);
                end
            end
        end
    end
    channels(logical(badSession)) = cell(1,sum(badSession));
    
end

end