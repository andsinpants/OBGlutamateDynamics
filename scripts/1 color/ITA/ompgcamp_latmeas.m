%% import OMP-GCaMP files (only 1x odor?)
clc
clearvars
close all

%% load odor file
%for g=1:7 %loop through each avgfile from exp.
disp('Select Odor Gun File')
uiopen('load')

for i = 1:length(myplotdata.avgfile.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.avgfile.roi(i);
end

prompt = ('Is this a 3 to 7s odor presentation?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    myplot_times = [myplot_r(1).odor(1).avgtrial.time]; %%4.5sec to 8.5sec have no 'odor' field. WHY?
            for rnum = 1:length(myplot_r) %rnum=roi number
            allplot_odor(:,rnum) = myplot_r(rnum).odor.avgtrial.series';
            end
end
allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

figure;imagesc(allplot_odor);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%clearvars -except allplot_odor timemat_odor myplot_onum onum

figure;plot(allplot_odor');
%% check for significance
%% find times of stimulus onset and offset
for i=1
    if tf==1
        odoron = 3;
        odoroff = 7;
    end
    if tf==0
         prompt = 'Odor on?';
         odoron = input(prompt);
         prompt = 'Odor off?';
         odoroff = input(prompt);
    end
end
% prompt = 'Odor on?';
% odoron = input(prompt);
% prompt = 'Odor off?';
% odoroff = input(prompt);
prompt = 'Odorant (use concentration and initials)?';
odornum = input(prompt,'s');
timemat_odoron = find(timemat_odor(1,:)>=odoron);
timemat_odoroff = find(timemat_odor(1,:)<=odoroff);
todoron = timemat_odoron(1);
todoroff = timemat_odoroff(end);
tpreodor = 1:timemat_odoron(1);

%% find mean of odor responses (ODOR)
for g = 1:size(allplot_odor)
    meanvals_odor(g) = mean(allplot_odor(g,todoron:todoroff),2);
    meanvals_preodor(g) = mean(allplot_odor(g,tpreodor)); 
end
meanvals_odor = meanvals_odor';
meanvals_preodor = meanvals_preodor';

%% find 95_max of odor responses (ODOR)
for g = 1:size(allplot_odor)
    maxvals_odor(g) = max(allplot_odor(g,todoron:todoroff));
end
maxvals_odor = maxvals_odor';
max95_odor = .95*(maxvals_odor);

%% find 5pct of odor responses (ODOR)
for g = 1:size(allplot_odor)
qprctile_odor(g,:) = qprctile(allplot_odor(g,todoron:todoroff), [5.0 95.0]);
end
pct5_odor = qprctile_odor(:,1);

%% find std of odor responses (ODOR)
for g = 1:size(allplot_odor);
std_preodor(g,:) = std(allplot_odor(g,1:todoron));
end
%close all%remove after writing script
%% zscore data
for g= 1:size(allplot_odor)
    zdata_odor(g,:) = ((allplot_odor(g,:)-meanvals_preodor(g))./std_preodor(g));
end
 %% smooth responses display%%
%   allplot_smth = [];
% for i=1:length(allplot_odor(:,1))
%     ys=allplot_odor(i,:);
%     smoothb=fastsmooth(ys,15,3,1);
%     allplot_smth=[allplot_smth;smoothb];
% end
allplot_smth = allplot_odor; %no need to smooth epi signals
figure;plot(allplot_smth');

%% zscore data
for g= 1:size(allplot_smth,1)
    zdata_odorsmth(g,:) = ((allplot_smth(g,:)-meanvals_preodor(g))./std_preodor(g));
end
figure;plot(zdata_odorsmth')

%% get std descriptive stats
for g=1:size(zdata_odorsmth,1)
b_odor(g,:) = zdata_odorsmth(g,todoron:todoroff);
b_odormax(g) = max(b_odor(g,:));
b_odormin(g) = min(b_odor(g,:));
b_odormean(g) = mean(b_odor(g,:));
b_odorpct5(g,:) = qprctile(zdata_odorsmth(g,todoron:todoroff), [5.0 95.0]);
%b_odorpct5(g) = b_odorpct5(:,1)
end
b_odormax = b_odormax';
b_odormin = b_odormin';
b_odorpct5 = b_odorpct5(:,1);
b_odormean = b_odormean';
b_odormax95 = b_odormax(:,:)*.95;

%% find significant zscored data by responses above 4SD 
aboveposThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmth(:,1))
    posThreshold = 4; %%change for threshold
    aboveposThreshold(i,:) = (zdata_odorsmth(i,todoron:todoroff))>= posThreshold; 
end

for i=1:length(zdata_odorsmth(:,1))
    spanLengthpos = regionprops(aboveposThreshold(i,:));
    if isempty(spanLengthpos)
     spanLengthpos_mat(i) = 0;
    else
    spanLengthpos = spanLengthpos.Area;
    spanLengthpos_mat(i) = spanLengthpos;
    end
end

goodspans_exc = find(spanLengthpos_mat>=30); %%change for span length
excsigrois = allplot_smth(goodspans_exc,:);
%end %uncomment when finished

%% now for sig. raw odor traces and max/5th percetile
for l=1:size(goodspans_exc,2)
g = goodspans_exc(l);
a_odor(g,:) = allplot_smth(g,todoron:todoroff);
a_odormax(g) = max(a_odor(g,:));
a_odormin(g) = min(a_odor(g,:));
a_odormean(g) = mean(a_odor(g,:));
a_odorpct5(g,:) = qprctile(allplot_smth(g,todoron:todoroff), [5.0 95.0]);
%b_odorpct5(g) = b_odorpct5(:,1)
end
a_odormax = a_odormax';
a_odormin = a_odormin';
a_odorpct5 = a_odorpct5(:,1);
a_odormean = a_odormean';
a_odormax95 = a_odormax(:,:)*.95;

%% extract risetimes, onset latency  
for i=1
    if tf==0
    for r=1:size(excsigrois(:,1))
        timemat = 0:1/150:(length(excsigrois(r,:))-1)/150;
        risetime(excsigrois(r,todoron:todoron+450));
        [R,LT,UT,LL,UL] = risetime(excsigrois(r,todoron:todoron+450));
        risetimes(r) = R;
        pct10times(r) = LT;
        pct90times(r) = UT;
        onsetlat(r) = pct10times(r)-excsigrois(r,todoron);
        onsetlatrd(r) = round(onsetlat(r));
        onsetlatspan = timemat(1:onsetlatrd(r));
        onsetlattimes(r) = onsetlatspan(end);
        pct10timesrd(r) = round(pct10times(r));
        pct90timesrd(r) = round(pct90times(r));
        pct10dfoverf(r) = excsigrois(r,pct10timesrd(r));
        pct90dfoverf(r) = excsigrois(r,pct90timesrd(r));
        peakvalues(r) = max(excsigrois(r,todoron:todoron+450));
    end
    end
    if tf==1
    for r=1:size(excsigrois(:,1))
        timemat = 0:1/150:(length(excsigrois(r,:))-1)/150;
        risetime(excsigrois(r,todoron:todoroff));
        [R,LT,UT,LL,UL] = risetime(excsigrois(r,todoron:todoroff));
        risetimes(r) = R;
        pct10times(r) = LT;
        pct90times(r) = UT;
        onsetlat(r) = pct10times(r)-excsigrois(r,todoron);
        onsetlatrd(r) = round(onsetlat(r));
        onsetlatspan = timemat(1:onsetlatrd(r));
        onsetlattimes(r) = onsetlatspan(end);
        pct10timesrd(r) = round(pct10times(r));
        pct90timesrd(r) = round(pct90times(r));
        pct10dfoverf(r) = excsigrois(r,pct10timesrd(r));
        pct90dfoverf(r) = excsigrois(r,pct90timesrd(r));
        peakvalues(r) = max(excsigrois(r,todoron:todoroff));
    end
    end
end

for r=1:size(excsigrois(:,1))
    timemat = 0:1/150:(length(excsigrois(r,:))-1)/150;
    otimes = timemat(todoron:todoroff);
    rtrng = timemat(1:risetimes(r));
    rt_times(r) = rtrng(end);
    pct90round(r) = round(pct90times(r));
    time90(r) = otimes(pct90round(r));
    pct10round(r) = round(pct10times(r));
    time10(r) = otimes(pct10round(r));
end

rt_times = rt_times';
onsetlattimes = onsetlattimes';
peakvalues = peakvalues';

%close all
%clearvars -except rt_times onsetlattimes peakvalues excsigrois goodspans_exc 
%% chose ROIs and save
% answer = inputdlg('Enter ROI numbers',...
%              'What excitatory ROIs do you want to keep for ', [1 50]);
% roiinput = str2num(answer{1});
% 
% for k=1:(length(roiinput))
%     n = roiinput(k);
%     ind = excsigrois(n,:);
%     allplot_smth(k,:) = ind;
% end

prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');

for i=1
    if tf==1
    save ('ompgcamp_lat_peak.mat', 'rt_times', 'onsetlattimes', 'peakvalues', 'excsigrois', 'goodspans_exc')
    if tf==0
        break
    end
    end
    clearvars -except rt_times onsetlattimes peakvalues excsigrois goodspans_exc
end

