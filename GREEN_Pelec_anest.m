%% IMPORT AND PROCESS DORIC DATA FROM DORIC FILE
%  JB 20/7/2023 (jbonaventura@ub.edu)

clear
close all

filename = uigetfile({'*.doric'}, 'Pick the FP file'); %select a file (open split down menu and select all file to show the CSV files)

name = input("Enter a name: ", 's');
if isempty(name)
    name = strsplit(filename,".");
    name = string(name(1));
end
RAWDATA = ExtractDataAcquisition(filename);

%attr = h5readatt(filename, '/DataAcquisition/FPConsole', 'DifferenceMasterStartToFirstData')

GREEN = RAWDATA(4).Data(1).Data;
ISOSBESTIC = RAWDATA(5).Data(1).Data;
%RED = RAWDATA(3).Data(2).Data  ;
Times = RAWDATA(4).Data(2).Data  ;
sr = round(length(Times)/Times(end)); % update to your sampling rate
%sr= 1205;

Fig1 = figure; 
subplot (2,1,1)
plot (Times, GREEN, 'g--');hold on;
%plot( Times, RED, 'r--');hold on;
plot (Times, ISOSBESTIC, 'm--')

% smooth the control trace using the running average
%window_size = round(0.2*sr); % in (C*sr), C is the window size in seconds
%smooth_ISOS = conv(ISOSBESTIC, ones(1, window_size) / window_size, 'same');  
smooth_ISOS = smoothdata(ISOSBESTIC, 'sgolay', 100);

%smooth the GREEN and RED
GREEN_raw = GREEN;
GREEN = smoothdata(GREEN, 'sgolay', 100);

%RED_raw = RED;
%RED = smoothdata(RED, 'sgolay', 100);

%% ACTIVATE to IMPORT DIGITAL OR ANALOG SIGNALS if NECESSARY
% data sources will need to be defined

%DIin1pre = RAWDATA(5).Data(1).Data;
%DIin4pre = RAWDATA(4).Data(2).Data;
%AOut1pre = RAWDATA(3).Data(1).Data;
AOut4pre = RAWDATA(2).Data(3).Data;

% fix the timing/align with locked-in data
time2 = RAWDATA(2).Data(4).Data;

AOut4 = interp1(time2, AOut4pre, Times, 'nearest', 'extrap');


%% TRIM DATA TO REMOVE ARTIFACTS

trimON = 1; %seconds to remove at the start
trimOFF = 1; %seconds to remove at the end

GREEN = GREEN(trimON*sr:(length(Times)-trimOFF*sr));
%RED = RED(trimON*sr:(length(Times)-trimOFF*sr));
CON = smooth_ISOS(trimON*sr:(length(Times)-trimOFF*sr));
T = Times(trimON*sr:(length(Times)-trimOFF*sr));

% if you have imported signals, activate/rename:
%Pelec = AOut4(trimON*sr:(length(AOut4)-trimOFF*sr));
%Popto = AOut1(trimON*sr:(length(AOut4)-trimOFF*sr));
%IN1 = DIin1(trimON*sr:(length(DIin1)-trimOFF*sr));
%IN2 = DIin2(trimON*sr:(length(DIin2)-trimOFF*sr));


subplot(2,1,1)
plot (T, GREEN, 'g','LineWidth', 2); hold on;
plot(T, CON, 'm','LineWidth', 2)
%plot(T, RED, 'r','LineWidth', 2)
%plot (T, Pelec*0.1)


%% correct using isosbestic

% SCALE ISOSBESTIC CONTROL TO BIOSENSOR SIGNAL LEVELS
% it uses the minima of the signal values
 [~, gminloc] = findpeaks (-GREEN, 'MinPeakDistance', 2*sr); %find minima locs
 p = polyfit (CON(gminloc), GREEN(gminloc), 1); %fits a a 1st order equation y = a*x + b
 %use polyfit to find the parameters

% baseline has no peaks, ingore the minima, take only baseline period (~20s)
%base = 20*sr; %the first 20s of data
%p = polyfit (CON(1:base), GREEN(1:base), 1); %fits a a 1st order equation y = a*x + b


a = p(1); b = p(2);

%scale the control
greenCON = CON*a + b;
%clear p a b
cGREEN = (GREEN - greenCON)./greenCON;

%if isosbestic correction causes noise, ignore it, use cGREEN2
%cGREEN2 = ((GREEN -mean(GREEN(1:base)))./mean(GREEN(1:base)));
% same for red channel
%[~, rminloc] = findpeaks (-RED, 'MinPeakDistance', 2*sr); %find minima locs
%p = polyfit (CON(rminloc), RED(rminloc), 1); %fits a a 1st order equation y = a*x + b
%a = p(1); b = p(2);

%scale the control
%redCON = CON*a + b;
%clear p a b
%cRED = (RED -redCON)./redCON;

subplot (2,1,1)
plot (T, greenCON, 'm','LineWidth', 2); 

subplot (2,1,2)
%plot(T, cRED, 'r','LineWidth', 2); hold on;
plot (T, cGREEN, 'g','LineWidth', 2); hold on
%plot (T, cGREEN2, 'b','LineWidth', 1); hold on

%plot (T, Pelec*0.02-0.01, 'k');
legend('Corrected w/ iso', 'Corrected with baseline', 'Stim')

%% trim around stimulus

pre = 5; % seconds before and after
post = 20;

[~, STIMLOC] = findpeaks (double(Pelec));

event = cGREEN(STIMLOC(1)-pre*sr:STIMLOC(1)+post*sr);
p_event = Pelec(STIMLOC(1)-pre*sr:STIMLOC(1)+post*sr);
t_event = (-pre:1/sr:post);
plot (t_event+T(STIMLOC(1)), event, 'k--', 'LineWidth', 1)

% downsample (to Fs) to export 
Fs = 30; %new sampling rate
tds_event = resample (event, Fs, sr)';

figure
plot (t_event, event);hold on;
plot (t_event, p_event*0.01-0.01)

%% quantify peaks

[ampl, loc, width, prom] = findpeaks (event, ...
    'MinPeakProminence',0.05,'MinPeakWidth', 0,'MinPeakDistance',0);

%if peaks are too close, will need to calculate slope manually
%activate to troubleshoot
%figure
%findpeaks (event, ...
%    'MinPeakProminence',0.008,'MinPeakWidth', 0,'MinPeakDistance',0, ...
%    'Annotate', 'extent');


peak_time = t_event(loc)*1000; peak_time = peak_time'; %in ms
width = width/sr; % in s

% if using single pulses
%[~, pulse_loc] = findpeaks(p_event);

% if using trains
[~, pulse_loc] = findpeaks(p_event);
ppt = 10; %pulses per train
pulse_loc = pulse_loc(1:ppt:end);

rise_time = (loc-pulse_loc)/sr*1000; % in ms
amp_increase = ampl-event(pulse_loc);
slope = amp_increase ./ rise_time*1000; % in dFF/s




peak_tbl = table (peak_time, ampl, amp_increase, prom, width, rise_time, slope);
disp (peak_tbl)
%% save figures and selected data
% activate/deatctivate your option

%save(name, 'p_event', 'event', 't_event', 'tds_event'); 
save(name, 'p_event', 'event', 't_event', 'tds_event', 'peak_tbl'); 

%savefig (Fig1, strcat(name,"_Fig1"))
%writetable (Gpeaktbl,strcat("summaryPeaksGREEN_",name,".xls"))
%writetable (Rpeaktbl,strcat("summaryPeaksRED_",name,".xls"))



