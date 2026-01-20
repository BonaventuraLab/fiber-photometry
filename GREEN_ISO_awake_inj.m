%% IMPORT AND PROCESS DORIC DATA FROM DORIC FILE
%  (jbonaventura@ub.edu)

% update: uses ALS calculated baseline to fit isosbestic
% update 20/01/2026: bin peak data

clear
close all

filename = uigetfile({'*.doric'}, 'Pick the FP file'); %select a file (open split down menu and select all file to show the CSV files)

name = input("Enter a name: ", 's');
if isempty(name)

    name = strsplit(filename,".");
    name = string(name(1));
end
RAWDATA = ExtractDataAcquisition(filename);

GREEN = RAWDATA(3).Data(2).Data;
ISOSBESTIC = RAWDATA(4).Data(2).Data;
% RED = RAWDATA(5).Data(1).Data  ;
Times = RAWDATA(3).Data(1).Data  ;
sr = round(length(Times)/Times(end)); % update to your sampling rate
%sr= 1205;

Fig1 = figure; 
subplot (2,1,1)
plot (Times, GREEN, 'g--');hold on;
%plot( Times, RED, 'r--');hold on;
plot (Times, ISOSBESTIC, 'm--')

% smooth the control trace using the running average
window_size = round(0.5*sr); % in (C*sr), C is the window size in seconds
smooth_ISOS = conv(ISOSBESTIC, ones(1, window_size) / window_size, 'same');  

%smooth the GREEN and RED
GREEN_raw = GREEN;
GREEN = smoothdata(GREEN, 'sgolay', 0.1*sr);

%RED_raw = RED;
%RED = smoothdata(RED, 'sgolay', 100);

%% ACTIVATE to IMPORT DIGITAL OR ANALOG SIGNALS if NECESSARY
% data sources will need to be defined

%DIin1pre = RAWDATA(5).Data(1).Data;
%DIin4pre = RAWDATA(4).Data(2).Data;
%AOut1pre = RAWDATA(3).Data(1).Data;
%AOut4pre = RAWDATA(2).Data(3).Data;

% fix the timing/align with locked-in data
%time2 = RAWDATA(2).Data(4).Data;

%AOut4 = interp1(time2, AOut4pre, Times, 'nearest', 'extrap');

%attr = h5readatt(filename, '/DataAcquisition/FPConsole', 'DifferenceMasterStartToFirstData')

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


%OLD METHOD 
    %[~, gminloc] = findpeaks (-GREEN(1:100*sr), 'MinPeakDistance', 3*sr); %find minima locs
    %p = polyfit (CON(gminloc), GREEN(gminloc), 1); %fits a a 1st order equation y = a*x + b

% if baseline has no peaks, ingore the minima, take only baseline period (~20s)
    %base = 20*60*sr; %the first 20min of data
    %p = polyfit (CON(1:base), GREEN(1:base), 1); %fits a a 1st order equation y = a*x + b

% NEW METHOD

% Define parameters for ALS
lambda = sr*1e4; % Smoothing parameter (higher value = smoother baseline)
p = 0.01; % Asymmetry parameter (adjust for how peaks affect baseline)


baseline_G = als_baseline(GREEN, lambda, p);

flatzone = 60*sr;  %define the zone in which the signal should be flat
                    % i.e. before the injection/treatment

p = polyfit (CON(1:flatzone), baseline_G(1:flatzone), 1);

a = p(1); b = p(2);

%scale the control
greenCON = CON*a + b;
%clear p a b
cGREEN = (GREEN - greenCON)./greenCON;
%cGREEN2 = (GREEN - baseline_G)./baselineG; % high pass filter, removes trends

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
plot (T, baseline_G, 'y','LineWidth', 1); 

subplot (2,1,2)
%plot(T, cRED, 'r','LineWidth', 2); hold on;
plot (T, cGREEN, 'g','LineWidth', 2, 'DisplayName','Corrected w/ iso'); hold on
%plot (T, cGREEN2, 'b','LineWidth', 1, 'DisplayName','Corrected w/ baseline'); hold on

%plot (T, Pelec*0.02-0.01, 'k');
legend('Corrected w/ iso')

%% trim around injection
% keep 15 pre injection and 30 min post
 before = 15 * 60 *sr;
 after = 30 * 60 *sr;
% 
 INJstart = 20*60 + 5;
 INJend = 20*60 + 15;
% 
 [~, idxSTART] = min(abs(T - INJstart));
 [~, idxEND] = min(abs(T - INJend));
% 
 cGREENbefore = cGREEN(idxSTART-before:idxSTART);
 cGREENafter = cGREEN(idxEND:idxEND+after);
 cGREENcut = vertcat(cGREENbefore(2:end), cGREENafter);
% 
% %cGREENcut(1357727:1359845)=NaN;
% %cGREENcut(1361776:1366921)=NaN;
% 
figure
newT = (-10:1/sr/60:30)';
plot (newT, cGREENcut);hold on;
plot ([0 0], [-0.1 0.2], 'k--')
ylim([-0.1 0.2])


%% quantify peaks

[ampl, loc, width, prom] = findpeaks (cGREENcut, ...
    'MinPeakProminence',0.01,'MinPeakWidth', 0,'MinPeakDistance',0);

width = width/sr; % in s
pk_time = newT(loc);
peak_tbl = table (loc, ampl, prom, width);
plot (pk_time, ampl, 'or')
title('dFF corrected data')
%% zScore
    % and its peaks

zScore = (cGREENcut-median(cGREENbefore))./std(cGREENbefore);

[zampl, zloc, zwidth, zprom] = findpeaks (zScore, ...
    'MinPeakProminence',1.6,'MinPeakWidth', 0,'MinPeakDistance',0);

zwidth = zwidth/sr; % in s
zpk_time = newT(zloc);
zscore_peak_tbl = table (zloc, zampl, zprom, zwidth);

figure 
plot (newT, zScore); hold on;
plot (zpk_time, zampl, 'or')
title('z-score corrected data')

% bin peak data
timebin = 5;        % bin size in minutes
t0 = -10;            % STARTING TIME in minutes
nbins = 8;           % number of bins
gbins = nan(nbins,4);
for i = 1:nbins
    t_start = t0 + timebin*(i-1);
    t_end   = t0 + timebin*i;
    idx = zpk_time >= t_start & zpk_time < t_end;
    greenwid = zwidth(idx);
    greenamp = zprom(idx);
    gbins(i,1)= t_end;
    gbins(i,2) = numel(greenwid) / (timebin);  % peaks per min
    gbins(i,3) = mean(greenwid, 'omitnan');
    gbins(i,4) = mean(greenamp, 'omitnan');
end

disp(['Time bin / Pks per min / Width / Prom ...in: ', num2str(timebin), ' minute bins'])
disp(gbins)




%% save figures and selected data
% activate/deatctivate your option

save (name)
%save(name, "T", 'newT', 'cGREENcut', "cGREEN", "peak_tbl","greenCON",'GREEN', 'sr'); 

%savefig (Fig1, strcat(name,"_Fig1"))
%writetable (Gpeaktbl,strcat("summaryPeaksGREEN_",name,".xls"))
%writetable (Rpeaktbl,strcat("summaryPeaksRED_",name,".xls"))



%%
% ALS function
function baseline = als_baseline(y, lambda, p)
    L = length(y);
    D = diff(speye(L), 2);
    w = ones(L, 1);
    for i = 1:10
        W = spdiags(w, 0, L, L);
        Z = W + lambda * (D' * D);
        baseline = Z \ (w .* y(:));
        w = p * (y > baseline) + (1 - p) * (y <= baseline);
    end
end


