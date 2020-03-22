function[LFP] = DAVIS_MultiTrace_Manual_Artifact_Reject(LFP,window_size_sec,overwrite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function[LFP] = DAVIS_Manual_Artifact_Reject(LFP,window_size_sec,overwrite)
% LFP        - filename or LFP structure from DAVIS_Pre_Process. If left
% empty or blank it will pull up a gui to select the file.
% window_size - How big a window to look at each time. DEfauklt 30 seconds
% overwrite   - Whether to keep the bad intervals currently in the file(0)
% or completely replace them (1)
% Basic GUI stuff:
% Run with no arguments, this will let you open up a file and start using
% the Define Artifact button to select bad points. You can also click
% outside of the axis to select points not currently in view, like before
% the start of a recording, or in the next or preveious window. At the
% begining it detects the current filename and writes it into the LFP
% structure, for future refr ence. It'll ask you ifyou want to use the new
% filename if you've changed it before. Finally, at the end, you can
% overwrite the original file, or just end the program. The LFP structure
%  will still be output.
%  Mattenator 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LFP_fullpath = [];
if nargin < 1 | isempty(LFP);
    
    [LFP,LFP_fullpath] = uigetfile('*.mat','Select multiple LFP files with ctrl or shift','Multiselect','on');
    
end
if strcmp(LFP,'All')
    LFP = dir('*CSC*.mat');
    LFP = {LFP.name};
end
if ischar(LFP)
    LFP = {LFP}
end
if nargin < 2
    window_size_sec = 300;
end
if nargin < 3
    overwrite = 0;
end
for i_lfp = 1:length(LFP);
    filename = char(LFP{i_lfp});
    LFP{i_lfp} = load([LFP_fullpath filename]);
    try if ~strcmp(LFP{i_lfp}.name,filename);
            button = questdlg('Filename and LFP.name do not agree. Overwrite LFP.name?','Filename Overwrite');
            switch button;
                case 'Yes'
                    LFP{i_lfp}.name = filename;
                    disp('Setting name value of LFP to current file name.')
                case 'No'
                case 'Cancel'
                    return
            end
        end
    catch
        LFP{i_lfp}.name = filename;
        disp('No name in LFP.')
        disp('Setting name value of LFP to current file name.')
    end
    
    
    if any(diff(LFP{i_lfp}.timestamps)<0)
        disp('Timestamps are incorrect, were donwsampled at some point. Making up new ones.');
        LFP{i_lfp}.timestamps = linspace(LFP{i_lfp}.timestamps(1),LFP{i_lfp}.timestamps(end),length(LFP{i_lfp}.values));
    end
    TS{i_lfp} = LFP{i_lfp}.timestamps-LFP{i_lfp}.timestamps(1);
    values{i_lfp} = LFP{i_lfp}.values;
    badvals{i_lfp} = nan(size(LFP{i_lfp}.values));
    try LFP{i_lfp}.bad_intervals;
    catch
        LFP{i_lfp}.bad_intervals = [1 2];
    end
    if isempty(overwrite) || overwrite == 0;
        for i_interval = 1:size(LFP{i_lfp}.bad_intervals,1)
            values{i_lfp}(...
                floor(LFP{i_lfp}.bad_intervals(i_interval,1))...
                :ceil(LFP{i_lfp}.bad_intervals(i_interval,2)))...
                = NaN();
            badvals{i_lfp}(...
                floor(LFP{i_lfp}.bad_intervals(i_interval,1))...
                :ceil(LFP{i_lfp}.bad_intervals(i_interval,2)))...
                = LFP{i_lfp}.values(...
                floor(LFP{i_lfp}.bad_intervals(i_interval,1))...
                :ceil(LFP{i_lfp}.bad_intervals(i_interval,2)));
        end
    else
        LFP{i_lfp}.bad_intervals = [1 2];
    end
end
point_jump = round(window_size_sec*LFP{1}.sFreq);

ixis = [1:point_jump:length(LFP{1}.values) length(LFP{1}.values)];
interf = figure;
i = 2;
figure(interf)
for figix = 1:length(LFP)
    hero(figix) = subplot(length(LFP),1,figix)
    plot(TS{figix}(ixis(i-1):ixis(i)),values{figix}(ixis(i-1):ixis(i)),'b')
    hold on
    plot(TS{figix}(ixis(i-1):ixis(i)),badvals{figix}(ixis(i-1):ixis(i)),'r')
    hold on
    xlabel('Time (s)')
    ylabel({['Ch ' num2str(LFP{figix}.Channel)]; 'mV'})
end
movegui(interf,'west')
hb = uicontrol('style','pushbutton');
set(hb,'position',[1 1 120 20])
set(hb,'string','Define Artifact')
set(hb,'callback',{@define_artifact,interf})
nw = uicontrol('style','pushbutton');
set(nw,'position',[1 21 120 20])
set(nw,'string','Next Window')
set(nw,'callback',{@next_window,interf})
pw = uicontrol('style','pushbutton');
set(pw,'position',[121 21 120 20])
set(pw,'string','Preveious Window')
set(pw,'callback',{@prev_window,interf})
ga = uicontrol('style','pushbutton');
set(ga,'position',[121 1 120 20])
set(ga,'string','Global Artifact')
set(ga,'callback',{@define_global_artifact,interf})
c = uicontrol('style','pushbutton');
set(c,'position',[241 1 120 20])
set(c,'string','Cancel')
set(c,'callback',{@exit,interf})
m = uicontrol('style','pushbutton');
set(m,'position',[241 21 120 20])
set(m,'string','Mirror Axes')
set(m,'callback',{@mirror,interf})
t = timer;
t.TimerFcn = @active_color
t.ExecutionMode = 'fixedRate';
t.Period = .25;
t.ErrorFcn = @(~,~) beep
start(t)
linked = 0;
    function mirror(m,~,~)
        if linked == 0;
            linkaxes(hero);
            linked = 1;
            set(m,'string','Separate Axes')
        else
            linkaxes(hero,'off')
            linked = 0;
            set(m,'string','Mirror Axes')
        end
    end

    function active_color(~,~)
        old = hero == gca;
        set(hero(old),'Color',[0.8 1 0.8]);
        set(hero(~old),'Color',[1 1 1]);
        
    end
    function define_artifact(hb, ~, interf)
        figure(interf)
        sub_ix = hero == gca;
        hold on
        title('Click Artifact Start')
        [xf(1),~] = ginput(1);
        %x(1) = xf(1) +ixis(i-1);
        [~, x(1)] = min(abs(TS{sub_ix}-xf(1)));
        line([xf(1) xf(1)],ylim,'Color',[0 1 0])
        title('Click Artifact End')
        [xf(2),~] = ginput(1);
        %x(2) = xf(2) +ixis(i-1);
        [~, x(2)] = min(abs(TS{sub_ix}-xf(2)));
        line([xf(2) xf(2)],ylim,'Color',[1 0 0])
        title('Right Click if OK')
        [~, ~, button] = ginput(1);
        set(hero(sub_ix),'Color',[1 1 1])
        if button == 3;
            LFP{sub_ix}.bad_intervals = [LFP{sub_ix}.bad_intervals; x;];
            LFP{sub_ix}.bad_intervals(LFP{sub_ix}.bad_intervals < 1) = 1;
            LFP{sub_ix}.bad_intervals(LFP{sub_ix}.bad_intervals > length(LFP{sub_ix}.values)) = length(LFP{sub_ix}.values);
        else
        end
        for i_interval = 1:size(LFP{sub_ix}.bad_intervals,1)
            values{sub_ix}(...
                floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)))...
                = NaN();
            badvals{sub_ix}(...
                floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)))...
                = LFP{sub_ix}.values(...
                floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)));
        end
        hold off
        plot(TS{sub_ix}(ixis(i-1):ixis(i)),values{sub_ix}(ixis(i-1):ixis(i)),'b')
        hold on
        plot(TS{sub_ix}(ixis(i-1):ixis(i)),badvals{sub_ix}(ixis(i-1):ixis(i)),'r')
        hold on
        xlabel('Time (s)')
        ylabel({['Ch ' num2str(LFP{sub_ix}.Channel)]; 'mV'})
    end
    function define_global_artifact(hb, ~, interf)
        figure(interf)
        old = hero == gca;
        hold on
        axes(hero(1));
        title('Click Artifact Start')
        axes(hero(old));
        [xf(1),~] = ginput(1);
        %x(1) = xf(1) +ixis(i-1);
        for sub_ix = 1:length(LFP);
            axes(hero(sub_ix));
            [~, x(1)] = min(abs(TS{sub_ix}-xf(1)));
            line([xf(1) xf(1)],ylim,'Color',[0 1 0])
        end
        axes(hero(1));
        title('Click Artifact End')
        axes(hero(old));
        [xf(2),~] = ginput(1);
        for sub_ix = 1:length(LFP);
            axes(hero(sub_ix));%x(2) = xf(2) +ixis(i-1);
            [~, x(2)] = min(abs(TS{sub_ix}-xf(2)));
            line([xf(2) xf(2)],ylim,'Color',[1 0 0])
        end
        axes(hero(1));
        title('Right Click if OK')
        axes(hero(old));
        [~, ~, button] = ginput(1);
        if button == 3;
            for sub_ix = 1:length(LFP);
                LFP{sub_ix}.bad_intervals = [LFP{sub_ix}.bad_intervals; x;];
                LFP{sub_ix}.bad_intervals(LFP{sub_ix}.bad_intervals < 1) = 1;
                LFP{sub_ix}.bad_intervals(LFP{sub_ix}.bad_intervals > length(LFP{sub_ix}.values)) = length(LFP{sub_ix}.values);
            end
        else
            
        end
        for sub_ix = 1:length(LFP);
            for i_interval = 1:size(LFP{sub_ix}.bad_intervals,1)
                values{sub_ix}(...
                    floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                    :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)))...
                    = NaN();
                badvals{sub_ix}(...
                    floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                    :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)))...
                    = LFP{sub_ix}.values(...
                    floor(LFP{sub_ix}.bad_intervals(i_interval,1))...
                    :ceil(LFP{sub_ix}.bad_intervals(i_interval,2)));
            end
            axes(hero(sub_ix));
            hold off
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),values{sub_ix}(ixis(i-1):ixis(i)),'b')
            hold on
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),badvals{sub_ix}(ixis(i-1):ixis(i)),'r')
            hold on
            xlabel('Time (s)')
     
        ylabel({['Ch ' num2str(LFP{sub_ix}.Channel)]; 'mV'})
        end
        axes(hero(old));
    end

    function next_window(nw, ~, interf)
        old = hero == gca;
        figure(interf)
        
        i = i+1;
        for sub_ix = 1:length(LFP);
            axes(hero(sub_ix));
            hold off
            try
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),values{sub_ix}(ixis(i-1):ixis(i)),'b')
                hold on
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),badvals{sub_ix}(ixis(i-1):ixis(i)),'r')
                hold on
                      
                xlabel('Time (s)')
             
        ylabel({['Ch ' num2str(LFP{sub_ix}.Channel)]; 'mV'})
            catch
                axes(hero(1));
                title('End of Recording. Click Cancel')
                axes(hero(sub_ix));
            end
        end
        axes(hero(old));
    end
    function prev_window(nw, ~, interf)
        figure(interf)
        old = hero == gca;
        hold off
        i = i-1;
        for sub_ix = 1:length(LFP);
            axes(hero(sub_ix));
            hold off
            try
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),values{sub_ix}(ixis(i-1):ixis(i)),'b')
                hold on
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),badvals{sub_ix}(ixis(i-1):ixis(i)),'r')
                hold on
                xlabel('Time (s)')
             
        ylabel({['Ch ' num2str(LFP{sub_ix}.Channel)]; 'mV'})
            catch
                axes(hero(1));
                title('No Previous Window.')
                axes(hero(old));
                i = i+1;
                break
            end
        end
        axes(hero(old));
    end
    function exit(c, ~, interf)
        button = questdlg('Save these bad intervals and changes to file?','LFP Overwrite');
        switch button;
            case 'Yes'
                for sub_ix = 1:length(LFP);
                    saveme = LFP{sub_ix};
                    save([LFP_fullpath LFP{sub_ix}.name],'-struct','saveme')
                end
            case 'No'
                disp('The LFP structure is still output into your workspace.')
            case 'Cancel'
               
                return
        end
        stop(t)
        close all
        return
    end
end