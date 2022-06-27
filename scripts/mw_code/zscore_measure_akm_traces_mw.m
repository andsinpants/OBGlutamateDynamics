preodor=1:500;
odoron=301:825;
numtrials = 3960;
t1=350:425;   %use for calculating t2-t1/tmax? use max of 500 msec post-odor onset.
t2=750:825;  %use max of last 500 msec before odor offset
for i = 1:numtrials;
    trace=results_mat(i,6:1507);
    t2t1tmax=(max(trace(t2))-max(trace(t1)))/results_mat(i,3);
    results_mat(i,1508)=t2t1tmax;
end;
zerotrials=isnan(results_mat(:,5));
results_mat(zerotrials,:)=[];
