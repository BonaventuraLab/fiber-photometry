%% analyze external events
delay = 855*1e-3;
Tvid = RAWDATA(2).Data.Data + delay;
GREENvid = interp1(T, cGREEN, Tvid); 


% Set up the Import Options and import the event data
opts = spreadsheetImportOptions("NumVariables", 4);
% Specify sheet and range
opts.Sheet = "video";
opts.DataRange = "A2:D41";
% Specify column names and types
opts.VariableNames = ["Frame", "Time", "event", "type"];
opts.VariableTypes = ["double", "double", "char", "char"];
% Specify variable properties
opts = setvaropts(opts, ["event", "type"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["event", "event"], "EmptyFieldRule", "auto");
% Import the data
events = readtable(uigetfile('*.xls*'), opts, "UseExcel", false);
clear opts

events.Time = Tvid(events.Frame);

F_SIstart = events.Frame (strcmp(events.event, 'START')& strcmp(events.type, 'SI'));
F_SIstop = events.Frame (strcmp(events.event, 'STOP')& strcmp(events.type, 'SI'));
F_ECstart = events.Frame (strcmp(events.event, 'START')& strcmp(events.type, 'EC'));
F_ECstop = events.Frame (strcmp(events.event, 'STOP')& strcmp(events.type, 'EC'));
%%
figure
plot (Tvid, GREENvid, 'g', 'LineWidth', 2);hold on
plot ([Tvid(F_SIstart) Tvid(F_SIstop)], [0 0], 'k-', 'LineWidth', 2)
plot ([Tvid(F_ECstart) Tvid(F_ECstop)], [0 0], 'r-', 'LineWidth', 2)

%% plot average
Fs = round(length(RAWDATA(2).Data.Data)/RAWDATA(2).Data.Data(end));
data = GREENvid;
pre = 5; %time (s) before event
pe = pre; pre = pre * Fs;
post = 5; %time (s) after event
po = post; post = post * Fs;

STIMLOC = F_SIstart;
clearvars event
for ii = 1:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        event(1:interval,ii)= data(in:fi);
    else
        col = size(event); r = col(1,2);
        r = r + 1;
        event(1:interval,r)= data(in:fi);
    end
end

event = event(1:col(1,1)-1,:);
clearvars ii col r in fi

S= size(event);
Enorm = ones (S);
    
for jj = 1:S(1,2)
    
    med = median(event([1:pre],jj));sd = std(event([1:pre],jj));
        for ii = 1:S(1,1)
        Enorm(ii,jj) = ((event(ii,jj)-med)./med); %dFF
        %Enorm(ii,jj) = ((events(ii,jj)-med)); %dF
        zEnorm(ii,jj) = ((event(ii,jj)-med)./sd); %zScore
    end
end 
clearvars ii jj 


figure
subplot (2,1,1)
Eavg = mean(zEnorm', 'OmitNaN')';
Eerr1 = mean(zEnorm', 'OmitNaN')' + std(zEnorm', 'OmitNaN')';
Eerr2 = mean(zEnorm', 'OmitNaN')' - std(zEnorm', 'OmitNaN')';
ax = gca;

fill ([1:length(Eerr1') fliplr(1:length(Eerr2'))], [Eerr1' fliplr(Eerr2')], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on
plot(Eavg, '-b', 'LineWidth', 2);
set(ax,  'XLim', [0 interval+1],'XTick', 1:pre:interval+1);
set(ax, 'XTickLabel', -pe:pe:po);
ax.YLabel.String = 'dFF';
title('Social Interaction')


%empty cage
pre = 5; %time (s) before event
pe = pre; pre = pre * Fs;
post = 5; %time (s) after event
po = post; post = post * Fs;

STIMLOC = F_ECstart;
clearvars event
for ii = 1:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        event(1:interval,ii)= data(in:fi);
    else
        col = size(event); r = col(1,2);
        r = r + 1;
        event(1:interval,r)= data(in:fi);
    end
end

event = event(1:col(1,1)-1,:);
clearvars ii col r in fi

S= size(event);
Enorm = ones (S);
    
for jj = 1:S(1,2)
    
    med = median(event([1:pre],jj));sd = std(event([1:pre],jj));
        for ii = 1:S(1,1)
        Enorm(ii,jj) = ((event(ii,jj)-med)./med); %dFF
        %Enorm(ii,jj) = ((events(ii,jj)-med)); %dF
        zEnorm(ii,jj) = ((event(ii,jj)-med)./sd); %zScore
    end
end 
clearvars ii jj 


subplot(2,1,2)
Eavg = mean(zEnorm', 'OmitNaN')';
Eerr1 = mean(zEnorm', 'OmitNaN')' + std(zEnorm', 'OmitNaN')';
Eerr2 = mean(zEnorm', 'OmitNaN')' - std(zEnorm', 'OmitNaN')';
ax = gca;

fill ([1:length(Eerr1') fliplr(1:length(Eerr2'))], [Eerr1' fliplr(Eerr2')], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on
plot(Eavg, '-b', 'LineWidth', 2);
set(ax,  'XLim', [0 interval+1],'XTick', 1:pre:interval+1);
set(ax, 'XTickLabel', -pe:pe:po);
ax.YLabel.String = 'dFF';
title ('Empty Cage')
%% filter Gpeaktbl
T_SIstart = events.Time (strcmp(events.event, 'START')& strcmp(events.type, 'SI'));
T_SIstop = events.Time (strcmp(events.event, 'STOP')& strcmp(events.type, 'SI'));
T_ECstart = events.Time (strcmp(events.event, 'START')& strcmp(events.type, 'EC'));
T_ECstop = events.Time (strcmp(events.event, 'STOP')& strcmp(events.type, 'EC'));

filtered_SI = table(); 
for i = 1:length(T_SIstart)
    onset = T_SIstart(i);
    offset = T_SIstop(i);
    
     idx = Gpeaktbl.Gpk_time >= onset & Gpeaktbl.Gpk_time <= offset; 
    if any(idx) % Check if any events match the condition
        filtered_SI = [filtered_SI; Gpeaktbl(idx, :)]; 
    end 
end
clear idx

filtered_EC = table(); 
for i = 1:length(T_ECstart)
    onset = T_ECstart(i);
    offset = T_ECstop(i);

   idx = Gpeaktbl.Gpk_time >= onset & Gpeaktbl.Gpk_time <= offset; 
    if any(idx) % Check if any events match the condition
        filtered_EC = [filtered_EC; Gpeaktbl(idx, :)]; 
    end  
end
clear idx

%other
T_ALLstart = events.Time (strcmp(events.event, 'START'));
T_ALLstop = events.Time (strcmp(events.event, 'STOP'));
filtered_other = Gpeaktbl; % Initialize with the full table

for i = 1:length(T_ALLstart)
    onset = T_ALLstart(i);
    offset = T_ALLstop(i);
    
    % Find indices of events *not* within the condition
    idx = ~(filtered_other.Gpk_time >= onset & filtered_other.Gpk_time <= offset); 
    
    % Keep only the events that do not meet the current condition
    filtered_other = filtered_other(idx, :); 
end

clear idx; 

%% calculate time
timeSI = sum(T_SIstop-T_SIstart);
timeEC = sum(T_ECstop-T_ECstart);

disp(['Time in Social Interaction: ', num2str(timeSI), ' s'])
disp(['Time in Empty Cage: ', num2str(timeEC), ' s'])
disp(['Total time: ', num2str(Tvid(end)), ' s'])