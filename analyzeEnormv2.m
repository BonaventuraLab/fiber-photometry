%% Analyze Enorm

% please load events file
%pEavg = mean(pEnorm')';
%data = pEavg;

data = Enorm;
%data = pEnorm;
%data= filtered_signal;
before = pre/sr;
after = post/sr;
%data = [baseline ketamine];
siz = size(data);
traces = siz (2);
fs = sr;

% smooth the data with the running average

%window_size = round(1*fs);
%for ii = 1:traces
%    smooth_data (:,ii) = conv(data(:,ii), ones(1, window_size) / window_size, 'same');  
%end

%data = smooth_data;


event_size = after + before;
%fs = length(data)/event_size;



%create time vector (not used anymore)
times = [-before:1/fs:after]';
times = times(1:length(data));

%find peaks


tbl = table (NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
         'VariableNames',{'loc' 'pks' 'w' 'p' 'tau', 'AUC', 'AUC_adjusted'});  
blank_row = table (NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
         'VariableNames',{'loc' 'pks' 'w' 'p' 'tau', 'AUC', 'AUC_adjusted'});     
Fig5 = figure;
for ii = 1:traces   
   
    %for not normalized events (remove baseline)
    %data(:,ii)=data(:,ii)-mean(data(1:before,ii));


     [pks, loc, w, p] = findpeaks (data(:,ii), fs, ...
         'MinPeakProminence',0.05,'MinPeakWidth', 1,'MinPeakDistance',0);
     
     %warning, use this for stimulated peaks
     %start of the first peak is the start of the stimulation
     pk_start = before*fs;
     
     %start of the first peak is the start of the trial
   %pk_end pk_start = 1 %this will lead to errors

     plot (data(:,ii)); hold on
     
     for j = 1:numel(loc)

     tau(j) = find(data(loc(j)*fs:end,ii) <= pks(j)/2, 1)/fs;
     
     pk_end = int16(find(data(loc(j)*fs:end,ii) <= 0.05*pks(j), 1)+loc(j)*fs);
    
     baseline = linspace(data(pk_start, ii), data(pk_end,ii), pk_end - pk_start + 1);
     plot (pk_start:pk_start+(numel(baseline)-1), baseline)
     plot (loc(j)*fs+before,pks(j), 'vk')
     %use this for normalized signal using 0 as a baseline
     AUC(j) = sum(data(pk_start:pk_end,ii))/fs;
        
     %use this for baseline defined as the start/end of the peak
     %might lead to errors
       
     peak_adjusted = data(pk_start:pk_end, ii) - baseline'; 
     AUC_adjusted(j) = sum(peak_adjusted)/fs;
    

     pk_start = pk_end; %update start of the next peak
     end


    xticks([1:fs:length(data)])
    xticklabels([-before:after])

     tau = tau';
     AUC = AUC';
     AUC_adjusted = AUC_adjusted';
     loc = loc - before;
     tbl_temp = table(loc, pks, w, p, tau, AUC, AUC_adjusted);
     tbl = vertcat(tbl,tbl_temp)
     %tbl = vertcat(tbl,blank_row)
     clearvars pks loc w p tau AUC p AUC_adjusted
    
end

%% 


writetable (tbl,strcat("peaksDA2_",name,".xls"))
%writetable (tbl,strcat("peaksCalcium_",name,".xls"))

