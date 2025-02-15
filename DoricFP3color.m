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

GREEN = RAWDATA(1).Data(2).Data;
ISOSBESTIC = RAWDATA(2).Data(2).Data;
RED = RAWDATA(3).Data(2).Data  ;
Times = RAWDATA(1).Data(1).Data  ;
sr = 1205; % update to your sampling rate

Fig1 = figure; 
subplot (2,1,1)
plot (Times, GREEN, 'g--');hold on;
plot( Times, RED, 'r--');hold on;
plot (Times, ISOSBESTIC, 'm--')

% smooth the control trace using the running average
window_size = round(0.5*sr); % in (C*sr), C is the window size in seconds
smooth_ISOS = conv(ISOSBESTIC, ones(1, window_size) / window_size, 'same');  

%smooth the GREEN and RED
GREEN_raw = GREEN;
GREEN = smoothdata(GREEN, 'sgolay', 100);

RED_raw = RED;
RED = smoothdata(RED, 'sgolay', 100);

%% ACTIVATE to IMPORT DIGITAL OR ANALOG SIGNALS if NECESSARY
% data sources will need to be defined

%DIin1pre = RAWDATA(5).Data(1).Data;
%DIin4pre = RAWDATA(4).Data(2).Data;
%AOut1pre = RAWDATA(3).Data(1).Data;
%AOut4pre = RAWDATA(4).Data(3).Data;

%timeshit = round(length(AOut4pre)-length(Times));

%AOut1 =  (AOut1pre(1+timeshit/2:length(AOut1pre)-timeshit/2));
%AOut4 =  (AOut4pre(1+timeshit/2:length(AOut4pre)-timeshit/2));
%DIin1 =  (DIin1pre(1+timeshit/2:length(DIin1pre)-timeshit/2));
%DIin4 =  (DIin4pre(1+timeshit/2:length(DIin4pre)-timeshit/2));


%% TRIM DATA TO REMOVE ARTIFACTS

trimON = 2; %seconds to remove at the start
trimOFF = 2; %seconds to remove at the end

GREEN = GREEN(trimON*sr:(length(Times)-trimOFF*sr));
RED = RED(trimON*sr:(length(Times)-trimOFF*sr));
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
plot(T, RED, 'r','LineWidth', 2)
%plot (T, Pelec*0.1)


%% correct using isosbestic

% SCALE ISOSBESTIC CONTROL TO BIOSENSOR SIGNAL LEVELS
% it uses the minima of the signal values
[~, gminloc] = findpeaks (-GREEN, 'MinPeakDistance', 2*sr); %find minima locs
% use polyfit to find the parameters
p = polyfit (CON(gminloc), GREEN(gminloc), 1); %fits a a 1st order equation y = a*x + b
a = p(1); b = p(2);

%scale the control
greenCON = CON*a + b;
%clear p a b
cGREEN = (GREEN - greenCON)./greenCON;

% same for red channel
[~, rminloc] = findpeaks (-RED, 'MinPeakDistance', 2*sr); %find minima locs
p = polyfit (CON(rminloc), RED(rminloc), 1); %fits a a 1st order equation y = a*x + b
a = p(1); b = p(2);

%scale the control
redCON = CON*a + b;
clear p a b
cRED = (RED -redCON)./redCON;

subplot (2,1,1)
plot (T, greenCON, 'm','LineWidth', 2); 

subplot (2,1,2)
plot(T, cRED, 'r','LineWidth', 2); hold on;
plot (T, cGREEN, 'g','LineWidth', 2); hold on

%% CORRECT DECAY WITHOUT USING ISOSBESTIC

%[~, gminloc] = findpeaks (-GREEN, 'MinPeakDistance', 2*sr); %in case you didnt run the previous section
%decay = fit (T(gminloc), GREEN(gminloc), "exp2"); %check other options besides exp2
%coefs2 = coeffvalues (decay); 

%a=coefs2(1);b=coefs2(2);c=coefs2(3);d=coefs2(4);

%GREEN_f0 = a*exp(b*T) + c*exp(d*T);
%clearvars a b c d coefs2

%cniGREEN = (GREEN -GREEN_f0)./GREEN_f0;


%same with red channel 
%[~, rminloc] = findpeaks (-RED, 'MinPeakDistance', 2*sr); %in case you didnt run the previous section
%decay = fit (T(rminloc), RED(rminloc), "exp2"); %check other options besides exp2
%coefs2 = coeffvalues (decay); 
%a=coefs2(1);b=coefs2(2);c=coefs2(3);d=coefs2(4);

%RED_f0 = a*exp(b*T) + c*exp(d*T);
%clearvars a b c d coefs2

%cniRED = (RED -RED_f0)./RED_f0;
%subplot (2,1,2)
%figure
%plot(T, cniRED, 'r','LineWidth', 1); hold on;
%plot (T, cniGREEN, 'c','LineWidth', 1)

%% DETECT PEAKS

prom_threshold = 0.05; % adjust prominence threshold for the green Ch
%change if you want to use another corrected signal
[Gpk_dFFcon, Gloc, Gwid, Gprom] = findpeaks (cGREEN,  'MinPeakProminence',prom_threshold);
Gpk_width = Gwid/sr;
Gpk_time = T(Gloc);
Gpeaktbl = table(Gloc, Gpk_time, Gpk_dFFcon, Gpk_width, Gprom);
plot (Gpeaktbl.Gpk_time, Gpk_dFFcon, 'ob');hold on;

Gpkpm = length (Gpk_time)/(length(cGREEN)/(sr*60));
disp(strcat ('GREEN Peaks per minute:  ', num2str(Gpkpm)))



%prom_threshold = 210; % adjust prominence threshold for the red Ch

%[Rpk_dFFcon, Rloc, Rwid, Rprom] = findpeaks (cRED,  'MinPeakProminence',prom_threshold);
%Rpk_width = Rwid/sr;
%Rpk_time = T(Rloc);
%Rpeaktbl = table(Rloc, Rpk_time, Rpk_dFFcon, Rpk_width, Rprom);
%plot (Rpeaktbl.Rpk_time, Rpk_dFFcon, 'ok')

%Rpkpm = length (Rpk_time)/(length(cRED)/(sr*60));
%disp(strcat ('RED Peaks per minute:  ', num2str(Rpkpm)))

%% save figures and selected data
% activate/deatctivate your option

save(name)   
%savefig (Fig1, strcat(name,"_Fig1"))
%writetable (Gpeaktbl,strcat("summaryPeaksGREEN_",name,".xls"))
%writetable (Rpeaktbl,strcat("summaryPeaksRED_",name,".xls"))



