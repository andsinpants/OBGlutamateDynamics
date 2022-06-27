%% clear console and vars
clc
clearvars
%% load odor file
disp('Select Lowest Conc. Odor Plotdata File')
uiopen('load')

for i = 1:length(myplotdata.file.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.file.roi(i);
end

for j = 1:length(myplot_r) 
allplot_odor1x(j,:) = myplot_r(j).odor.avgtrial.series';
end

myplot_times = myplot_r(1).odor.avgtrial.time;
% allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

figure;imagesc(allplot_odor1x);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

figure;plot(allplot_odor1x');
%% find times of stimulus onset and offset
prompt = 'Odor on?';
odoron = input(prompt);
prompt = 'Odor off?';
odoroff = input(prompt);
prompt = 'Odorant (use concentration and initials)?';
odornum = input(prompt,'s');
timemat_odoron = find(timemat_odor(1,:)>=odoron);
timemat_odoroff = find(timemat_odor(1,:)<=odoroff);
todoron = timemat_odoron(1);
todoroff = timemat_odoroff(end);
tpreodor = 1:timemat_odoron(1);

%% find mean of odor responses (ODOR)
for g = 1:size(allplot_odor1x)
    meanvals_odor(g) = mean(allplot_odor1x(g,todoron:todoroff),2);
    meanvals_preodor(g) = mean(allplot_odor1x(g,tpreodor)); 
end
meanvals_odor = meanvals_odor';
meanvals_preodor = meanvals_preodor';

%% find 95_max of odor responses (ODOR)
for g = 1:size(allplot_odor1x)
    maxvals_odor(g) = max(allplot_odor1x(g,todoron:todoroff));
end
maxvals_odor = maxvals_odor';
max95_odor = .95*(maxvals_odor);

%% find 5pct of odor responses (ODOR)
for g = 1:size(allplot_odor1x)
qprctile_odor(g,:) = qprctile(allplot_odor1x(g,todoron:todoroff), [5.0 95.0]);
end
pct5_odor = qprctile_odor(:,1);

%% find std of odor responses (ODOR)
for g = 1:size(allplot_odor1x)
std_preodor(g,:) = std(allplot_odor1x(g,1:todoron));
end
%close all%remove after writing script
%% zscore data
for g= 1:size(allplot_odor1x)
    zdata_odor(g,:) = ((allplot_odor1x(g,:)-meanvals_preodor(g))./std_preodor(g));
end
 %% smooth responses display%%
  allplot_smth = [];
for i=1:length(allplot_odor1x(:,1))
    ys=allplot_odor1x(i,:);
    smoothb=fastsmooth(ys,15,3,1);
    allplot_smth=[allplot_smth;smoothb];
end

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

goodspans_exc1x = find(spanLengthpos_mat>=30); %%change for span length
excsigrois1x = allplot_smth(goodspans_exc1x,:);
%excroisz = zdata_odorsmth(goodspans,:)
[~,r] = sortrows(excsigrois1x,'descend');

%clearvars -except allplot_odor1x timemat_odor

%% sort excitatory ROIs by onset time after odor
for i = 1:length(excsigrois1x(:,1))
    if isempty(excsigrois1x)
    break
    end
    for i=1:length(excsigrois1x(:,1))
    sorted_excplotsmth1x = excsigrois1x(r,:);
    end
end

for j=1
    if isempty(excsigrois1x)
        break
    end
    figure;imagesc(sorted_excplotsmth1x),colorbar;title('Significantly Excited/Sorted ROIs');colormap(gca,nawhimar_auto);
    xlabel('Frames')
    ylabel('ROI Number')%df/f allplot of all ROIs
end
%% find significant zscored data by responses below 4SD 
belownegThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmth(:,1))
    negThreshold = -4; %%change for threshold
    belownegThreshold(i,:) = (zdata_odorsmth(i,todoron:todoroff))<= negThreshold; 
end

for i=1:length(zdata_odorsmth(:,1))
    spanLengthneg = regionprops(belownegThreshold(i,:));
    if isempty(spanLengthneg)
     spanLengthneg_mat(i) = 0;
    else
    spanLengthneg = spanLengthneg.Area;
    spanLengthneg_mat(i) = spanLengthneg;
    end
end

goodspans_sup1x = find(spanLengthneg_mat>=30); %%change for span length
supsigrois1x = allplot_smth(goodspans_sup1x,:);
%excroisz = zdata_odorsmth(goodspans,:)
[~,r] = sortrows(supsigrois1x,'ascend');
%figure;imagesc(excrois),colorbar;%%excitatory ROIs

%% sort suppressed ROIs by onset time after odor presentation
for i = 1:length(supsigrois1x(:,1))
    if isempty(supsigrois1x)
    break
    end
    for i=1:length(supsigrois1x(:,1))
    sorted_supplotsmth1x = supsigrois1x(r,:);
    end
end

for j=1
    if isempty(supsigrois1x)
        break
    end
    figure;imagesc(sorted_supplotsmth1x),colorbar;title('Significantly Suppressed/Sorted ROIs');colormap(gca,nawhimar_auto);
    xlabel('Frames')
    ylabel('ROI Number')%df/f allplot of all ROIs
end

%% check if goodspans_exc1x and goodspans_sup1x are overlapping
for i=1
    overlapspans = intersect(goodspans_exc1x,goodspans_sup1x);
        if isempty(overlapspans)
            break
        end
        disp('Both excitatory and suppressive spans are overlapping!')
end

%% vertically cocantenate arrays
sigsortapodor_smth1x = [];
for j=1
    if isempty(sorted_excplotsmth1x)
    if isempty(sorted_supplotsmth1x)
        break
    end
sigsortapodor_smth1x=vertcat(sorted_excplotsmth1x,sorted_supplotsmth1x);
    end
end

for j=1
    if isempty(sorted_excplotsmth1x)
    if isempty(sorted_supplotsmth1x)
        break
    end
    figure;imagesc(sigsortapodor_smth1x),colorbar;title('Significantly Excited/Suppressed Sorted ROIs (todoron:todoroff)');colormap(gca,nawhimar_auto);
    end
end

%% normalize each exc and sup significantly sorted ROI matrix separately
normexc1x=[];
for i=1:length(sorted_excplotsmth1x(:,1))
    if ~exist('sorted_excplotsmth1x','var') || isempty(sorted_excplotsmth1x)
        break
    end
for j=1:size(sorted_excplotsmth1x,1) %% in progress
    normexc1x(j,:) = normalised_diff(sorted_excplotsmth1x(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normsup1x = [];
for i=1
    if ~exist('sorted_supplotsmth1x','var') || isempty(sorted_supplotsmth1x)
        break
    end
for j=1:size(sorted_supplotsmth1x,1) %% in progress
    normsup1x(j,:) = normalised_diff(sorted_supplotsmth1x(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

allplot_sigodorexc1x = [];
allplot_sigodorsup1x = [];
for i=1:length(goodspans_exc1x)
allplot_sigodorexc1x(i,:) = allplot_odor1x(goodspans_exc1x(i),:);
end
for i=1:length(goodspans_sup1x)
allplot_sigodorsup1x(i,:) = allplot_odor1x(goodspans_sup1x(i),:);
end
allplot_sigodor1x = vertcat(allplot_sigodorexc1x,allplot_sigodorsup1x);
%% now normalize each excitatory and suppressive matrix together
normall1x = vertcat(normexc1x,normsup1x);
figure;plot(normexc1x','r');hold on;plot(normsup1x','b'),title(['Significantly Excited and Suppressed ROIs Normalized and Sorted (Traces) {\color{red}Odor ID:}', num2str(odornum)])
%figure;imagesc(normall1x);colorbar;title(['Significantly Excited and Suppressed ROIs Normalized and Sorted {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
%c = colorbar;
%c.Label.String = 'Norm. \DeltaF/F';
%gca;xlabel({'Frames'});ylabel({'ROIs'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% create imagesc fig
figure;imagesc(normexc1x);colorbar;title(['Significantly Excited ROIs Normalized and Sorted (1x) {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
figure;imagesc(normsup1x);colorbar;title(['Significantly Suppressed ROIs Normalized and Sorted (1x) {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
figure;imagesc(normall1x);colorbar;title(['Significantly Excited and Suppressed ROIs Normalized and Sorted (1x) {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
%% change variable names
excrois_1xsmth = excsigrois1x;
suprois_1xsmth = supsigrois1x;
%% save variables
prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if ~exist('supsigrois1x','var') || isempty(supsigrois1x)
        save ('excandsuproiunsorted_1x.mat', 'excrois_1xsmth','allplot_sigodor1x')
     elseif tf==1 && ~isempty(supsigrois1x)
        save ('excandsuproisunsorted_1x.mat', 'excrois_1xsmth', 'suprois_1xsmth','allplot_sigodor1x')
    end
end

for i=1
    if ~exist('supsigrois1x','var') || isempty(supsigrois1x)
        save ('excandsuprois_sorted1x.mat', 'sorted_excplotsmth1x')
    elseif tf==1 && ~isempty(supsigrois1x)
        save ('excandsuprois_sorted1x.mat', 'sorted_excplotsmth1x', 'sorted_supplotsmth1x')
    end
end
      
for i=1
    if ~exist('supsigrois1x','var') || isempty(supsigrois1x)
        save ('excandsuprois_sortandnorm1x.mat', 'normexc1x')
    elseif tf==1 && ~isempty(supsigrois1x)
        save ('excandsuprois_sortandnorm1x.mat', 'normexc1x','normsup1x')
    end
end

for i=1
    if ~exist('supsigrois1x','var') || isempty(supsigrois1x)
        save ('excandsuprois_1x.mat', 'goodspans_exc1x','allplot_odor1x')
    elseif tf==1 && ~isempty(supsigrois1x)
        save ('excandsuprois_1x.mat', 'goodspans_exc1x', 'goodspans_sup1x')
    end
end
        
%% close all figures
prompt = ('Do you want to close all figures? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf == 1
    close all
else
end
%clearvars -except excrois1x suprois1x 