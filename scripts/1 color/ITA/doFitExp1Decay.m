close all
clear expdecay rmse peakvalues fitresults fitgoodness;
expdecay = zeros(size(rawchosentraces,1),1);
rmse = zeros(size(rawchosentraces,1),1);
peakvalues = zeros(size(rawchosentraces,1),1);
fitresults = cell(size(rawchosentraces,1),1);
fitgoodness = cell(size(rawchosentraces,1),1);
clear peakind r tmpind
for r = 1:size(rawchosentraces,1)
    times = 0:1/150:(length(rawchosentraces(r,2:end))-1)/150;%return to (r,:) if not using first column to index rois
    [peakvalues(r),peakind] = max(rawchosentraces(r,2:end));%return to (r,:) if not using first column to index rois
    [~,tmpind] = find(rawchosentraces(r,peakind:end)<= 0.9*peakvalues(r),1,'first');
    tmpind = peakind+tmpind-1;
    %[fitresults{r},fitgoodness{r}] = fit(times(tmpind:end)',rawchosentraces(r,tmpind:end)','exp1');
    [fitresults{r},fitgoodness{r}] = fit(times(tmpind:tmpind+150)',rawchosentraces(r,tmpind:tmpind+150)','exp1');
    figure; plot(times,rawchosentraces(r,2:end)); tmpax = axis; %return to (r,:) if not using first column to index rois
    hold on; plot(fitresults{r}); axis(tmpax); xlabel('Time (s)'), ylabel('{\Delta}F/F')
    text(1.5,rawchosentraces(r,tmpind),sprintf('decay constant = %f',fitresults{r}.b));
    text(1.5,rawchosentraces(r,tmpind+2),sprintf('ROI number = %f',rawchosentraces(r,1))); %added in so that the ROI can print with the decay constant
    expdecay(r) = fitresults{r}.b;
    reexpdecay(r) = -(1/expdecay(r));
    rmse(r) = fitgoodness{r}.rmse;
    %pause;
end
% akm start here
 todoron = 152; %1 sec; use 151 for non-roi indexed rois
 todoroff = 452; %3 sec; use 451 for non-roi indexed rois (use for 602 frames)
 %todoroff = 352; %3 sec; use 451 for non-roi indexed rois (use for 436 frames)
for r=1:size(rawchosentraces(:,1))
    timemat = 0:1/150:(length(rawchosentraces(r,2:end))-1)/150;
    [R,LT,UT] = risetime(rawchosentraces(r,2:end)); %change back to (r,2:end) if using unfiltered data
    if numel(R)>1 
        R=R(1);
    if isempty(R)
        R=0;
    end
    end
    if numel(LT)>1
        LT=LT(1);
    if isempty(LT)
        LT=0;
    end
    end
    if numel(UT)>1 
        UT=UT(1);
    if isempty(UT)
        UT=0;
    end
    end
    if R==[] & UT==[] & LT==[]
        continue
    end
    risetimes(r) = R;
    pct10times(r) = LT;
    pct90times(r) = UT;
    onsetlat(r) = pct10times(r)-rawchosentraces(r,todoron);
    onsetlatrd(r) = round(onsetlat(r));
    onsetlatspan = timemat(1:onsetlatrd(r));
    onsetlattimes(r) = onsetlatspan(end);
    pct10timesrd(r) = round(pct10times(r));
    pct90timesrd(r) = round(pct90times(r));
    pct10dfoverf(r) = rawchosentraces(r,todoron+pct10timesrd(r));
    pct90dfoverf(r) = rawchosentraces(r,todoron+pct90timesrd(r));
    peakvalues(r) = max(rawchosentraces(r,todoron:todoroff))*.90;
    rt_frames(r) = pct90times(r)-pct10times(r);
    rt_times(r) = rt_frames(r)/150;
    ttpeak_frames = find(rawchosentraces(r,todoron:end)>peakvalues(r));
    ttpeak_frame(r) = ttpeak_frames(1);
    ttpeak_times(r) = ttpeak_frame(r)/150;
    dur_fwhm(r) = fwhm(timemat,rawchosentraces(r,2:end));
end

% for r=1:size(rawchosentraces(:,1))
% % %     timemat = 0:1/150:(length(rawchosentraces(r,:))-1)/150;
% %     otimes = timemat(todoron:todoroff);
% %     rtrng = otimes(1:risetimes(r));
% % %     rtrng = timemat(1:risetimes(r));
% %     rt_times(r) = rtrng(end);
% %     pct90round(r) = round(pct90times(r));
% %     time90(r) = otimes(pct90round(r));
% %     pct10round(r) = round(pct10times(r));
% %     time10(r) = otimes(pct10round(r));
%     rt_frames = pct90times-pct10times;
%     rt_times = rt_frames/150;
% %     onlat_frames = 1+pct10times-todoron;
% %     onlat_times = onlat_frames/150;
%     ttpeak_frames = peakvalues-todoron;
%     ttpeak_times = ttpeak_frames/150;
%     dur_fwhm(r) = fwhm(timemat,rawchosentraces(r,2:end));
% end


%close all

%  for i=1:(size(rawchosentraces,1))%% now plot onsetlat todoron and top10% of values
x1=zeros(size(rawchosentraces,1),2);
x2=zeros(size(rawchosentraces,1),2);
x3=zeros(size(rawchosentraces,1),2);

 for i =1:size(rawchosentraces,1)
    x1(i,:) = [onsetlatrd(i) onsetlatrd(i)];
 end
 for i =1:size(rawchosentraces,1)
    x2(i,:) = [todoron todoron];
 end
 for i=1:size(rawchosentraces,1)
    x3(i,:) = [pct90timesrd(i),pct90timesrd(i)];
 end
 
 x1time = timemat(x1);
 x2time = timemat(x2);
 x3time = timemat(x3);
for i =1:size(rawchosentraces,1)
%     hold on
    figure(i),plot(timemat,rawchosentraces(i,2:end)); tmpax = axis;
line(x1time(i,:),ylim,'Color','b');
line(x2time(i,:),ylim,'Color','r');
line(x3time(i,:),ylim,'Color','k');
text(1.5,rawchosentraces(i,tmpind),sprintf('risetime = %f',rt_times(i)))
text(1.5,rawchosentraces(i,tmpind+1),sprintf('onset latency = %f',onsetlattimes(i)))
text(1.5,rawchosentraces(i,tmpind+2),sprintf('fwhm = %f',dur_fwhm(i)))
text(1.5,rawchosentraces(i,tmpind+3),sprintf('time2peak = %f',ttpeak_times(i)))
legend('Risetime (t10:t90)', 'Onset latency (todoron:t10)','Odor and Sniff (todoron)')
% xticks = [];
end

%% to use for zoom
% for i =1:size(rawchosentraces,1)
% %     hold on
%     figure(i),plot(timemat,rawchosentraces(i,2:end)); tmpax = axis;
% line(x1time(i,:)+1,ylim,'Color','b');
% line(x2time(i,:),ylim,'Color','r');
% line(x3time(i,:)+1,ylim,'Color','k');
% legend('Risetime (t10:t90)', 'Onset latency (todoron:t10)','Odor and Sniff (todoron)')
% xlim([0.5 2]);
% clear xlim
% end


% for
% f1(i) = x1(1,i);
% f2(i) = x2(1,i);
% f3(i) = x3(1,i);
% % hold on
% text(x1(i,:),[0 0],sprintf('pct90timesrd is %f',f1(i)));
% text(x2(i,:),[10 10],sprintf('odoron is 1 (151)','r'));
% text(x3(i,:),[15 15], sprintf('onsetlat is %f',f3(i)));
% legend('Onset Latency', 'Odor and Sniff', 'Top 10% of values')
% % hold off
% end

for i=1
    tf = size(rawchosentraces,1)<=27;
        if tf==1
        autoArrangeFigures();
        else 
        break
    end
end

for i=1
    tf = ttpeak_times<onsetlattimes;
    if tf==1
        disp('Time to peak is less than onset latency!')
    else
        break
    end
end
    
rt_times = rt_times';
dur_fwhm = dur_fwhm';
reexpdecay = reexpdecay';
onsetlattimes = onsetlattimes';
pct90dfoverf = pct90dfoverf';
ttpeak_times = ttpeak_times';


clearvars -except rt_times pct90dfoverf reexpdecay dur_fwhm onsetlattimes rawchosentraces ttpeak_times fitresults fitgoodness
%%clear fitresults fitgoodness;

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
%     [fitresults{r},fitgoodness{r}] = fit(times(tmpind:end)',normchosentraces(r,tmpind:end)','exp1');
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