% the DLC analysis file should be open and the area of interest selected
% make sure the frame rate (fps) is correct, if not.. this is useless
% when it prompts, select the excel file obtained from DoricFP3color.m

inzoneTTL= zeros(length(X),1);
inzoneTTL(inZone)=1;
time = (0:1/fps:(length(X)-1)/fps)';
%plot(time, inzoneTTL)
edges = diff(inzoneTTL); edges(length(inzoneTTL)+1) = 0; 

% a table with start and end times of every period it enters/exits the Zone
enters = time (edges == 1);
exits = time (edges == -1);

if size (enters) == size(exits)
    if enters(1) > exits(1)
        inzone_periods(:,1) = vertcat (0, enters(1:end-1));
        inzone_periods(:,2) = exits;
    else
        inzone_periods(:,1) = enters;
        inzone_periods(:,2) = exits;
    end
end

if size (enters) > size(exits)
    inzone_periods(:,1) = enters(1:end-1);
    inzone_periods(:,2) = exits;
end

if size (enters) < size(exits)
    inzone_periods(:,1) = vertcat (0, enters);
        inzone_periods(:,2) = exits;
end
peak_table = readtable(uigetfile({'*.xls'}, 'Pick the Excel file with peaks'))

event_times = peak_table.pk_time; %times of the events
events_within_zone = false(size(event_times));


for i = 1:numel(event_times)
    % Check if the event time falls within any of the state periods
    for j = 1:size(inzone_periods, 1)
        if event_times(i) >= inzone_periods(j, 1) && event_times(i) <= inzone_periods(j, 2)
            % The event occurs within the current state period
            events_within_zone(i) = true;
            break; % Exit the inner loop since the event has been classified
        end
    end
end

% Filter the events based on whether they occur within the state periods or not
peaks_inZone = peak_table(events_within_zone,:)
peaks_outZone = peak_table(events_within_zone==false,:);

writetable (peaks_inZone,strcat("peaksINzone_",name,".xls"))
writetable (peaks_outZone,strcat("peaksOUTzone_",name,".xls"))
