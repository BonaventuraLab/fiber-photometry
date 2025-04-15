%% Sort signals from DIin1 depending on their length
% requires data from DoricFPprocess
% signals from FED3 are stored in the variable IN
% this version does not plot graphs but puts all trials together
% and generates predictor vectors

data = (cRED);
IN = IN1;

edges = diff(IN); edges(length(edges)+1) = 0; 
rise = find (edges == 1);
fall = find (edges == -1);
len = (fall-rise);

pellet = rise(len >= 10 & len <= 14);
rightpoke = rise (len <= 3);
cue = rise (len >=5 & len <= 7);
leftpoke = rise(len >= 23);
randomized = randi (length(data), 15, 1);



%% JB (6/23/2017)- look up events
% uptaded in May 2023
% uses data generated from doricFPprocess.m
%stripped down version to use with FED

n = 1; %legacy from LUevents.m

sr = round(sr); %sampling rate of the recording
pre = 5; %time (s) before event
pe = pre; pre = pre * sr;
post = 5; %time (s) after event
po = post; post = post * sr;

%% trials aligned to pellet retrieval
STIMLOC = pellet;

clearvars events
for ii = 1:n:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        events(1:interval,ii)= data(in:fi);
    else
        col = size(events); r = col(1,2);
        r = r + 1;
        events(1:interval,r)= data(in:fi);
    end
end

events = events(1:col(1,1)-1,:);
clearvars  col r in fi

% normalization to baseline
Enorm = events;
pellet_Enorm = Enorm -mean(Enorm(1:pre,:),1); %subtract baseline (time before onset)


%% trials aligned to cue onset
STIMLOC = cue;

clearvars events
for ii = 1:n:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        events(1:interval,ii)= data(in:fi);
    else
        col = size(events); r = col(1,2);
        r = r + 1;
        events(1:interval,r)= data(in:fi);
    end
end

events = events(1:col(1,1)-1,:);
clearvars  col r in fi

% normalization to baseline
Enorm = events;
cue_Enorm = Enorm -mean(Enorm(1:pre,:),1); %subtract baseline (time before onset)


%% trials aligned to inactive poke
STIMLOC = rightpoke;

clearvars events
for ii = 1:n:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        events(1:interval,ii)= data(in:fi);
    else
        col = size(events); r = col(1,2);
        r = r + 1;
        events(1:interval,r)= data(in:fi);
    end
end

events = events(1:col(1,1)-1,:);
clearvars  col r in fi

% normalization to baseline
Enorm = events;
inactive_Enorm = Enorm -mean(Enorm(1:pre,:),1); %subtract baseline (time before onset)

%% random trials
STIMLOC = randomized;

clearvars events
for ii = 1:n:length(STIMLOC) %change if events are in another DI
    in = round (STIMLOC(ii)-pre);
    fi = round (STIMLOC(ii)+post);
    interval = fi - in + 1;
    if ii == 1       
        events(1:interval,ii)= data(in:fi);
    else
        col = size(events); r = col(1,2);
        r = r + 1;
        events(1:interval,r)= data(in:fi);
    end
end

events = events(1:col(1,1)-1,:);
clearvars  col r in fi

% normalization to baseline
Enorm = events;
random_Enorm = Enorm -mean(Enorm(1:pre,:),1); %subtract baseline (time before onset)

%% put all events together in 'trials' and make predictor vectors

trials = [pellet_Enorm cue_Enorm inactive_Enorm random_Enorm];

npellet = size (pellet_Enorm, 2);
ncue = size (cue_Enorm, 2);
ninactive = size (inactive_Enorm, 2);
nrandom = size (random_Enorm, 2);

X_pellet = zeros(size(trials, 2), 1);
X_pellet(1:npellet) = 1;
current_position = npellet+1;

X_cue = zeros(size(trials, 2), 1);
X_cue(current_position:current_position+ncue-1) = 1;
current_position = current_position+ncue+1;

X_inactive = zeros(size(trials, 2), 1);
X_inactive(current_position:current_position+ninactive-1) = 1;
current_position = current_position+ninactive+1;

X_random = zeros(size(trials,2), 1);
X_random(current_position:current_position+nrandom-1) = 1; 

%edit manually
X_drug = ones(size(trials, 2), 1);

%% save

clearvars ax p pe po med ii jj interval n S

save (strcat("all_events_DA_",name))


