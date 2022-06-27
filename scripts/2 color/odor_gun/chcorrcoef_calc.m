clc
clear all
close all

%% mutiselect files to convert
disp('Select combined dataset with Ch1 and Ch2')
uiopen('load')

%% loop through each and plot correlations across each timeseries, 95% max, 5th percentile, and only during odor onset
% start with complete timeseries
for i =1:size(x1)
tic
chccoef = corrcoef(x1(i,:),x2(i,:)); %change to numbers to when odorant is on?
chccvals = akmccoef(akmccoef~=1);
chtestvals(i,1) = akmccvals(1);
toc
end

chmean = mean(chtestvals);
chmovmean = movmean(chtestvals,5);

% correlation during odor presentation
for i =1:size(x1)
tic
chccoef_odor = corrcoef(x1(i,XXX:XXX),x2(i,XXX:XXX)); %change to numbers to when odorant is on?
chccvals_odor = akmccoef(akmccoef~=1);
chtestvals_odor(i,1) = akmccvals(1);
toc
end

chmean_odor = mean(chtestvals_odor);
chmovmean_odor = movmean(chtestvals_odor,5);

prompt = ('Do you want to plot corrcoef vals? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i =1
    if tf=1
        figure;plot(chmean),ylim(0:1),hold on,xlabel('G-O Pairs'),ylabel('Correlation value across entire concentration')
        plot(chmovemean,'k');
        figure;plot(chmean_odor),ylim(0:1),hold on,xlabel('G-O Pairs'),ylabel('Correlation values across odor duration')
        plot(chmovemean_odor,'k');
    else
    end
end

% now find peaks during odorant presentation
% now find troughs during odorant presentation


prompt = ('Do you want to save vars? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');