function[LFP] = DAVIS_MultiTrace_Manual_Artifact_Reject(LFP,LFP_fullpath,window_size_sec,overwrite)
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
%------------------------------------------------------------------------
% Using file split function
% select start and end point of the area WHERE THE GRAPH FLATLINES
% like you would while selecting artifacts
% Click split file to generate the files. WAIT till the done message pops
% up.
% If incorrect break points have been set, click reset breakpoints and save
% it does not show the changes immediately until it is reloaded.
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
    window_size_sec = 3000;
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
    brkpts{i_lfp} = nan(size(LFP{i_lfp}.values));
    try LFP{i_lfp}.bad_intervals;
    catch
        LFP{i_lfp}.bad_intervals = [1 2];
    end
    try LFP{i_lfp}.break_points;
    catch
        LFP{i_lfp}.break_points = [1 2];
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
        for i_interval = 1:size(LFP{i_lfp}.break_points,1)
            brkpts{i_lfp}(...
                floor(LFP{i_lfp}.break_points(i_interval,1))...
                :ceil(LFP{i_lfp}.break_points(i_interval,2)))...
                = LFP{i_lfp}.values(...
                floor(LFP{i_lfp}.break_points(i_interval,1))...
                :ceil(LFP{i_lfp}.break_points(i_interval,2)));        
        end
    else
        LFP{i_lfp}.bad_intervals = [1 2];
        LFP{i_lfp}.break_points = [1 2];
    end
end
point_jump = round(window_size_sec*LFP{1}.sFreq);

%ixis = [1:point_jump:length(LFP{1}.values) length(LFP{1}.values)];
interf = figure;
i = 2;
figure(interf)
for figix = 1:length(LFP)
    ixis = [1:point_jump:length(LFP{figix}.values) length(LFP{figix}.values)];
    hero(figix) = subplot(length(LFP),1,figix)
    plot(TS{figix}(ixis(i-1):ixis(i)),values{figix}(ixis(i-1):ixis(i)),'b')
    hold on
    plot(TS{figix}(ixis(i-1):ixis(i)),badvals{figix}(ixis(i-1):ixis(i)),'r')
    hold on
    plot(TS{figix}(ixis(i-1):ixis(i)),brkpts{figix}(ixis(i-1):ixis(i)),'g')
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
bp = uicontrol('style','pushbutton');
set(bp,'position',[362 21 120 20])
set(bp,'string','Set Breaks')
set(bp,'callback',{@set_breaks,interf})
bf = uicontrol('style','pushbutton');
set(bf,'position',[362 1 120 20])
set(bf,'string','Split File')
set(bf,'callback',{@split_file,interf})
rbp = uicontrol('style','pushbutton');
set(rbp,'position',[483 21 120 20])
set(rbp,'string','Reset Breakpoints')
set(rbp,'callback',{@reset_bp,interf})
t = timer;
t.TimerFcn = @active_color
t.ExecutionMode = 'fixedRate';
t.Period = .25;
t.ErrorFcn = @(~,~) beep
start(t)
linked = 0;
    
    function reset_bp(rbp,~,interf)
        for sub_ix = 1:length(LFP)
            LFP{sub_ix}.break_points = [1,2];
        end
    end
    
    function split_file(bf,~,interf)
        for i_lfp = 1:length(LFP)
        t_val = LFP{i_lfp}.values;
        t_timestamps = LFP{i_lfp}.timestamps;
        t_break_points = LFP{i_lfp}.break_points;
        t_bad_intervals = LFP{i_lfp}.bad_intervals;
        starting = 1;
        ending = 0;
        n_bad_intervals = [1,2];
        for idx = 2:length(LFP{i_lfp}.break_points)+1
            if idx == length(LFP{i_lfp}.break_points)+1
                ending = length(LFP{i_lfp}.timestamps);
            else
                ending = LFP{i_lfp}.break_points(idx,1);
            end
            for indx = 1:size(LFP{i_lfp}.bad_intervals,1)
                if LFP{i_lfp}.bad_intervals(indx,1) >= starting && LFP{i_lfp}.bad_intervals(indx,2) <= ending
                    n_bad_intervals = [n_bad_intervals;LFP{i_lfp}.bad_intervals(indx,1)-starting+1 LFP{i_lfp}.bad_intervals(indx,2)-starting+1];
                elseif LFP{i_lfp}.bad_intervals(indx,1) < starting && LFP{i_lfp}.bad_intervals(indx,2) > ending
                    n_bad_intervals = [n_bad_intervals;1 ending-starting+1];
                elseif (LFP{i_lfp}.bad_intervals(indx,1) >= starting && LFP{i_lfp}.bad_intervals(indx,1) <=ending) && LFP{i_lfp}.bad_intervals(indx,2) > ending
                    n_bad_intervals = [n_bad_intervals;LFP{i_lfp}.bad_intervals(indx,1)-starting+1 ending-starting+1];
                elseif LFP{i_lfp}.bad_intervals(indx,1) < starting && (LFP{i_lfp}.bad_intervals(indx,2) <= ending && LFP{i_lfp}.bad_intervals(indx,2) >= starting)
                    n_bad_intervals = [n_bad_intervals;1 LFP{i_lfp}.bad_intervals(indx,2)-starting+1];
                else
                    ;
                end
               
     
            end
            LFP{i_lfp}.values = LFP{i_lfp}.values(starting:ending);
            LFP{i_lfp}.timestamps = LFP{i_lfp}.timestamps(starting:ending); 
            LFP{i_lfp}.bad_intervals = n_bad_intervals;
            LFP{i_lfp}.break_points = [];
            n_bad_intervals = [];
                    saveme = LFP{i_lfp};
                    if idx-1 == 1
                        new_name = [LFP{i_lfp}.name(1:end-4) '.mat'];
                        save([LFP_fullpath new_name],'-struct','saveme')
                    else
                        new_name = [LFP{i_lfp}.name(1:end-4) '_000' num2str(idx-2) '.mat'];
                        save([LFP_fullpath new_name],'-struct','saveme')
                    end         
            LFP{i_lfp}.values = t_val;
            LFP{i_lfp}.timestamps = t_timestamps;
            LFP{i_lfp}.bad_intervals = t_bad_intervals;
            LFP{i_lfp}.break_points = t_break_points;
            if idx <= length(LFP{i_lfp}.break_points)
               starting = LFP{i_lfp}.break_points(idx,2);
            end
        end
       end
        msgbox('The file has been split.')
    end

    function set_breaks(bp,~,interf)
        figure(interf)
        old = hero == gca;
        hold on
        axes(hero(1));
        title('Click first breakpoint')
        axes(hero(old));
        [xf(1),~] = ginput(1);
        %x(1) = xf(1) +ixis(i-1);
        for sub_ix = 1:length(LFP);
            axes(hero(sub_ix));
            [~, x(1)] = min(abs(TS{sub_ix}-xf(1)));
            line([xf(1) xf(1)],ylim,'Color',[0 1 0])
        end
        axes(hero(1));
        title('Click second breakpoint')
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
                LFP{sub_ix}.break_points = [LFP{sub_ix}.break_points; x;];
                LFP{sub_ix}.break_points(LFP{sub_ix}.break_points < 1) = 1;
                LFP{sub_ix}.break_points(LFP{sub_ix}.break_points > length(LFP{sub_ix}.values)) = length(LFP{sub_ix}.values);
            end
        else
            
        end
        for sub_ix = 1:length(LFP);
            for i_interval = 1:size(LFP{sub_ix}.break_points,1)
                values{sub_ix}(...
                    floor(LFP{sub_ix}.break_points(i_interval,1))...
                    :ceil(LFP{sub_ix}.break_points(i_interval,2)))...
                    = NaN();
                brkpts{sub_ix}(...
                    floor(LFP{sub_ix}.break_points(i_interval,1))...
                    :ceil(LFP{sub_ix}.break_points(i_interval,2)))...
                    = LFP{sub_ix}.values(...
                    floor(LFP{sub_ix}.break_points(i_interval,1))...
                    :ceil(LFP{sub_ix}.break_points(i_interval,2)));
            end
            axes(hero(sub_ix));
            hold off
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),values{sub_ix}(ixis(i-1):ixis(i)),'b')
            hold on
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),badvals{sub_ix}(ixis(i-1):ixis(i)),'r')
            hold on
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),brkpts{sub_ix}(ixis(i-1):ixis(i)),'g')
            hold on
            xlabel('Time (s)')
     
        ylabel({['Ch ' num2str(LFP{sub_ix}.Channel)]; 'mV'})
        end
        axes(hero(old));
    end
        
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
        plot(TS{sub_ix}(ixis(i-1):ixis(i)),brkpts{sub_ix}(ixis(i-1):ixis(i)),'g')
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
            plot(TS{sub_ix}(ixis(i-1):ixis(i)),brkpts{sub_ix}(ixis(i-1):ixis(i)),'g')
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
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),brkpts{sub_ix}(ixis(i-1):ixis(i)),'g')
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
                plot(TS{sub_ix}(ixis(i-1):ixis(i)),brkpts{sub_ix}(ixis(i-1):ixis(i)),'g')
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