close all
clear expdecay rmse peakvalues fitresults fitgoodness;
expdecay = zeros(size(rawchosentraces,1),1);
rmse = zeros(size(rawchosentraces,1),1);
peakvalues = zeros(size(rawchosentraces,1),1);
fitresults = cell(size(rawchosentraces,1),1);
fitgoodness = cell(size(rawchosentraces,1),1);
clear peakind r tmpind
for r = 1:size(rawchosentraces,1)
    times = 0:1/150:(length(rawchosentraces(r,:))-1)/150;
    [peakvalues(r),peakind] = max(rawchosentraces(r,:));
    [~,tmpind] = find(rawchosentraces(r,peakind:end)<= 0.9*peakvalues(r),1,'first');
    tmpind = peakind+tmpind-1;
%     %[fitresults{r},fitgoodness{r}] = fit(times(tmpind:end)',rawchosentraces(r,tmpind:end)','exp2');     
    %[fitresults{r},fitgoodness{r}] = fit(times(tmpind:tmpind+150)',rawchosentraces(r,tmpind:tmpind+150)','exp1'); %     use for a184v,s,s72a, STA
    [fitresults{r},fitgoodness{r}] = fit(times(tmpind:tmpind+50)',rawchosentraces(r,tmpind:tmpind+50)','exp1'); %use for original iglusnfr
    figure; plot(times,rawchosentraces(r,:)); tmpax = axis;
    hold on; plot(fitresults{r}); axis(tmpax); xlabel('Time (s)'), ylabel('{\Delta}F/F')
    text(1.5,rawchosentraces(r,tmpind),sprintf('decay constant = %f',fitresults{r}.b));
    expdecay(r) = fitresults{r}.b;
    reexpdecay(r) = -(1/expdecay(r));
    rmse(r) = fitgoodness{r}.rmse;
    % akm start here
    time90=find(rawchosentraces(r,:)>=(.90*peakvalues(r)));%
    time90ind(r) = time90(1);
    time10=find(rawchosentraces(r,:)>=(.10*peakvalues(r)));
    time10ind(r) = time10(1);
    peak90times(r)=times(1,time90ind(r));%
    peak10times(r)=times(1,time10ind(r));
    risetimes(r) = peak90times(r)-peak10times(r);%
    onsetlat(r) = peak10times(r)-1;
    pause;
end
reexpdecay = reexpdecay';
risetimes = risetimes';
onsetlat = onsetlat';
clear peakind r tmpind times peak10times peak90times time10 time90 time90ind time10ind;
%clear fitresults fitgoodness;

%% try normalized traces now
% close all
% clear expdecay rmse peakvalues fitresults fitgoodness;
% expdecay = zeros(size(normchosentraces,1),1);
% rmse = zeros(size(normchosentraces,1),1);
% peakvalues = zeros(size(normchosentraces,1),1);
% fitresults = cell(size(normchosentraces,1),1);
% fitgoodness = cell(size(normchosentraces,1),1);
% clear peakind r tmpind
% for r = 1:size(normchosentraces,1)
%     times = 0:1/150:(length(normchosentraces(r,:))-1)/150;
%     [peakvalues(r),peakind] = max(normchosentraces(r,:));
%     [~,tmpind] = find(normchosentraces(r,peakind:end)<= 0.9*peakvalues(r),1,'first');
%     tmpind = peakind+tmpind-1;
%     [fitresults{r},fitgoodness{r}] = fit(times(tmpind:end)',normchosentraces(r,tmpind:end)','exp2');
%     figure; plot(times,normchosentraces(r,:)); tmpax = axis;
%     hold on; plot(fitresults{r}); axis(tmpax); xlabel('Time (s)'), ylabel('{\Delta}F/F')
%     text(1.5,normchosentraces(r,tmpind),sprintf('decay constant = %f',fitresults{r}.b));
%     expdecay(r) = fitresults{r}.b;
%     reexpdecay(r) = -(1/expdecay(r));
%     rmse(r) = fitgoodness{r}.rmse;
%     % akm start here
%     time90=find(normchosentraces(r,:)>=(.90*peakvalues(r)));%
%     time90ind(r) = time90(1);
%     time10=find(normchosentraces(r,:)>=(.10*peakvalues(r)));
%     time10ind(r) = time10(1);
%     peak90times(r)=times(1,time90ind(r));%
%     peak10times(r)=times(1,time10ind(r));
%     risetimes(r) = peak90times(r)-peak10times(r);%
%     onsetlat(r) = peak10times(r)-1;
%     pause;
% end
% reexpdecay = reexpdecay';
% risetimes = risetimes';
% onsetlat = onsetlat';
% clear peakind r tmpind times peak10times peak90times time10 time90 time90ind time10ind;
% clear fitresults fitgoodness;