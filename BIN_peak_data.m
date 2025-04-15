
timebin = 60; %in seconds
for i = 1:18
 greenwid = Gpk_width (Gpk_time > timebin*(i-1) & Gpk_time < timebin*i);
 greenamp = Gprom (Gpk_time > timebin*(i-1) & Gpk_time < timebin*i);
 gbins (i,1)= length (greenwid)/(timebin/60);
 gbins (i,2)= mean (greenwid);
 gbins (i,3)= mean (greenamp);
end
disp ('GREEN')
disp ('Pks per min / Width / Prom ...in X min bins')
disp (gbins)
%%
for i = 1:6
 redwid = Rpk_width (Rpk_time > timebin*(i-1) & Rpk_time < timebin*i);
 redamp = Rprom (Rpk_time > timebin*(i-1) & Rpk_time < timebin*i);
 rbins (i,1)= length (redwid)/(timebin/60);
 rbins (i,2)= mean (redwid);
 rbins (i,3)= mean (redamp);
end
disp ('RED')
disp ('Pks per min / Width / Prom ...in 5 min bins')
disp (rbins)