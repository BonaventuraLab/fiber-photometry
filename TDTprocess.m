%% IMPORT AND PROCESS TDT DATA
%  JB 17/12/2024 (jbonaventura@ub.edu)

clear
close all


% Import the data

filename = uigetdir( 'Pick the FP file'); %select a file (open split down menu and select all file to show the CSV files)

name = input("Enter a name: ", 's');
if isempty(name)
    name = strsplit(filename,"\");
    name = string(name(end));
end
RAWDATA = TDTbin2mat(filename);

GREEN = RAWDATA.streams.x65C.data;
ISOSBESTIC = RAWDATA.streams.x05C.data;
sr = round(RAWDATA.streams.x05C.fs);
Times = 1:length(GREEN);
Times = Times/sr+RAWDATA.streams.x05C.startTime;
figure;
subplot (2,1,1)
plot (Times, GREEN, 'g--');hold on;
plot (Times, ISOSBESTIC, 'm--')

% smooth the control trace using the running average

smooth_ISOS = smoothdata(ISOSBESTIC, 'sgolay', 100);

%smooth the GREEN
GREEN_raw = GREEN;
GREEN = smoothdata(GREEN, 'sgolay', 100);

%import relevant TTL data
PC0 = RAWDATA.epocs.PC0_.onset; %end of injection
PC1 = RAWDATA.epocs.PC1_.onset; %animal pick up

plot([PC0 PC0], [0 max(GREEN)], 'k--')
plot([PC1(1) PC1(1)], [0 max(GREEN)], 'k--')
%% TRIM DATA TO REMOVE ARTIFACTS

trimON = 10; %seconds to remove at the start
trimOFF = 2; %seconds to remove at the end

GREEN = GREEN(trimON*sr:(length(Times)-trimOFF*sr));
CON = smooth_ISOS(trimON*sr:(length(Times)-trimOFF*sr));
T = Times(trimON*sr:(length(Times)-trimOFF*sr));

[~, InjEnd] = min(abs(T - PC0));
[~, BlineEnd] = min(abs(T - PC1(1)));

subplot(2,1,1)
plot (T, GREEN, 'g','LineWidth', 2); hold on;
plot(T, CON, 'm','LineWidth', 2)


%% correct using isosbestic

% SCALE ISOSBESTIC CONTROL TO BIOSENSOR SIGNAL LEVELS
% it uses the minima of the signal values
% it only takes into account the baseline period

[~, gminloc] = findpeaks (-GREEN(1:BlineEnd), 'MinPeakDistance', 2*sr); %find minima locs
% use polyfit to find the parameters
p = polyfit (CON(gminloc), GREEN(gminloc), 1); %fits a a 1st order equation y = a*x + b
a = p(1); b = p(2);

%scale the control
greenCON = CON*a + b;
%clear p a b
cGREEN = (GREEN - greenCON)./greenCON;

subplot (2,1,1)
plot (T, greenCON, 'c','LineWidth', 1); 

subplot (2,1,2)
plot (T, cGREEN, 'g','LineWidth', 2); hold on
plot([PC0 PC0], [min(cGREEN) max(cGREEN)], 'k--')
plot([PC1(1) PC1(1)], [min(cGREEN) max(cGREEN)], 'k--')


%% CORRECT DECAY WITHOUT USING ISOSBESTIC
% ignore for now

%[~, gminloc] = findpeaks (-GREEN(1:BlineEnd), 'MinPeakDistance', 2*sr); %in case you didnt run the previous section
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

%% EXTRACT AND ALIGN BLINE AND DRUG PERIODS

%readjust times if necessary
bline=cGREEN(BlineEnd-600*sr:BlineEnd); %10 min of bline since picking up the animal
drug=cGREEN(InjEnd:InjEnd+1800*sr); %30 min after injection ends
alignedGREEN= horzcat(bline, drug); %reattach data

%downsample to usable sr (Fs)
Fs = 25; %25Hz
alignedGREEN = resample(alignedGREEN, Fs,sr);
alignedTime = [-10:1/Fs/60:30];
figure
plot (alignedTime, alignedGREEN, 'g', 'LineWidth', 2); hold on;
plot([0 0], [min(alignedGREEN) max(alignedGREEN)]);

%% DETECT PEAKS

prom_threshold = 0.1; % adjust prominence threshold for the green Ch
%change if you want to use another corrected signal
[Gpk_dFFcon, Gloc, Gwid, Gprom] = findpeaks (alignedGREEN,  'MinPeakProminence',prom_threshold);
Gpk_width = Gwid/Fs;
Gpk_time = alignedTime(Gloc);
peaktbl = table(Gloc', Gpk_time', Gpk_width', Gprom');
peaktbl.Properties.VariableNames = {'Location', 'Time', 'Width', 'Prominence'};
plot (Gpk_time, Gpk_dFFcon, 'ob');hold on;

Gpkpm = length (Gpk_time)/(length(cGREEN)/(Fs*60));
disp(strcat ('Global GREEN Peaks per minute:  ', num2str(Gpkpm)))

%% BIN peak DATA
timebin = 5; % 5 minute bins
for i = 1:8 % adjust if longer recording
 greenwid = Gpk_width (Gpk_time > timebin*(i-1)-10 & Gpk_time < timebin*i-10);
 greenamp = Gprom (Gpk_time > timebin*(i-1)-10 & Gpk_time < timebin*i-10);
 gbins (i,1)= length (greenwid)/(timebin);
 gbins (i,2)= mean (greenwid);
 gbins (i,3)= mean (greenamp);
end
disp ('GREEN')
disp ('Pks per min / Width / Prom ...in 5 min bins')
disp (gbins)



%% save figures and selected data

save(name)   
%savefig (Fig1, strcat(name,"_Fig1"))
%writetable (Gpeaktbl,strcat("summaryPeaksGREEN_",name,".xls"))
%writetable (Rpeaktbl,strcat("summaryPeaksRED_",name,".xls"))




