%% IMPORT AND PROCESS DORIC DATA and VIDEO and DOWNSAMPLE
% for frame matching
%  JB 30/11/2023 (jbonaventura@ub.edu)

clear
close all


% Import the data

filename = uigetfile({'*.doric'}, 'Pick the FP file'); %select a file (open split down menu and select all file to show the CSV files)

name = input("Enter a name: ", 's');
if isempty(name)
    name = strsplit(filename,".");
    name = string(name(1));
end
RAWDATA = ExtractDataAcquisition(filename);

GREEN = RAWDATA(2).Data(2).Data;
ISOSBESTIC = RAWDATA(3).Data(2).Data;
%RED = RAWDATA(2).Data(2).Data  ;
Times = RAWDATA(2).Data(1).Data  ;
sr = 1205;

FrameTimes = RAWDATA(1).Data.Data;

Fig1 = figure; 
subplot (2,1,1)
plot (Times, GREEN, 'g--');
hold on;
%plot (Times, RED, 'r--')
plot (Times, ISOSBESTIC, 'm--')

% smooth the control trace using the running average
%window_size = round(0.5*sr); % in (C*sr), C is the window size in seconds
%smooth_ISOS = conv(ISOSBESTIC, ones(1, window_size) / window_size, 'same');  

smooth_ISOS = smoothdata(ISOSBESTIC, 'sgolay', 100);
%smooth the GREEN and RED
GREEN_raw = GREEN;
GREEN = smoothdata(GREEN, 'sgolay', 100);

%RED_raw = RED;
%RED = smoothdata(RED, 'sgolay', 300);

% ACTIVATE to IMPORT SIGNALS if NECESSARY
% data sources will need to be defined

%DIin1pre = RAWDATA(5).Data(1).Data;
%AOut4pre = RAWDATA(4).Data(3).Data;
%timeshit = round(length(DIin1pre)-length(Times));
%AOut4 =  (AOut4pre(1+timeshit/2:length(AOut4pre)-timeshit/2));
%DIin1 =  (DIin1pre(1+timeshit/2:length(DIin1pre)-timeshit/2));



%% TRIM DATA TO REMOVE ARTIFACTS

trimON = 1 ; %seconds to remove at the start
trimOFF = 1 ; %seconds to remove at the end

GREEN = GREEN(trimON*sr:(length(Times)-trimOFF*sr));
%RED = RED(trimON*sr:(length(Times)-trimOFF*sr));
CON = smooth_ISOS(trimON*sr:(length(Times)-trimOFF*sr));
%CON = AIn1IsosbesticControl(trimON*sr:(length(Times)-trimOFF*sr));
T = Times(trimON*sr:(length(Times)-trimOFF*sr));

subplot(2,1,1)
plot (T, GREEN, 'g','LineWidth', 2); hold on;
plot(T, CON, 'm','LineWidth', 2)
%plot(T, RED, 'r','LineWidth', 2)


%% correct using isosbestic

% SCALE ISOSBESTIC CONTROL TO BIOSENSOR SIGNAL LEVELS
% it uses the minima of the signal values
[~, gminloc] = findpeaks (-GREEN, 'MinPeakDistance', 2*sr); %find minima locs
% use polyfit to find the parameters
p = polyfit (CON(gminloc), GREEN(gminloc), 1); %fits a a 1st order equation y = a*x + b
a = p(1); b = p(2);

%scale the control
greenCON = CON*a + b;
plot (T, greenCON, 'c', 'LineWidth', 2)
clear p a b
cGREEN = (GREEN - greenCON)./greenCON;

% same for red channel
%[~, rminloc] = findpeaks (-RED, 'MinPeakDistance', 2*sr); %find minima locs
%p = polyfit (CON(rminloc), RED(rminloc), 1); %fits a a 1st order equation y = a*x + b
%a = p(1); b = p(2);

%scale the control
%redCON = CON*a + b;
%clear p a b
%cRED = (RED -redCON)./redCON;

subplot (2,1,2)
%plot(T, cRED, 'r','LineWidth', 2); hold on;

plot (T, cGREEN, 'g','LineWidth', 2);hold on

%% downsample FP to fit video frames
% takes time
%for i = 1:length(FrameTimes)
%    [~, idx] = min(abs(T - FrameTimes(i)));
%    FrameTimeLoc (i) = idx;
%end
FrameTimeLoc = interp1(T, 1:numel(T), FrameTimes, 'nearest', 'extrap');

cGREEN_vid = cGREEN (FrameTimeLoc);
GREEN_vid = GREEN (FrameTimeLoc);
greenCON_vid = greenCON (FrameTimeLoc);
plot (FrameTimes, cGREEN_vid, 'k.','LineWidth', 1)
%% save
save (strcat("DS_vid__",name),"cGREEN_vid", "GREEN_vid","sr", "name", ...
    "greenCON_vid", "FrameTimes")
