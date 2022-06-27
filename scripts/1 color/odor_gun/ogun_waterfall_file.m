%% clear console and vars
clc
clearvars
%% load odor file
disp('Select Odor Gun File')
uiopen('load')

for i = 1:length(myplotdata.file.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.file.roi(i);
end

myplot_onum = [myplot_r(1).odor.number];
myplot_times = [myplot_r(1).odor(1).avgtrial.time];

disp(myplot_onum)
disp(1:12)
prompt = 'Odor Number?'
onum = input(prompt);
odornum = myplot_onum(onum)

for rnum = 1:length(myplot_r) %rnum=roi number
allplot_odor(:,rnum) = myplot_r(rnum).odor(onum).avgtrial.series';
end
allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

figure;imagesc(allplot_odor);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%clearvars -except allplot_odor timemat_odor myplot_onum onum

figure;plot(allplot_odor');
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
% 
% figure;plot(allplot_smth');

%% zscore data
for g= 1:size(allplot_odor,1) %change 'allplot_odor' to 'allplot_smth' for smooth responses
    zdata_odorsmth(g,:) = ((allplot_odor(g,:)-meanvals_preodor(g))./std_preodor(g));  %change 'allplot_odor' to 'allplot_smth' for smooth responses 
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
excsigrois = allplot_odor(goodspans_exc,:); %change 'allplot_odor' to 'allplot_smth' for smooth responses
%excroisz = zdata_odorsmth(goodspans,:)
[~,r] = sortrows(excsigrois,'descend');

%clearvars -except allplot_odor timemat_odor

%% sort excitatory ROIs by onset time after odor
sorted_excplotsmth = [];
for i = 1:length(excsigrois(:,1))
    if isempty(excsigrois)
    break
    end
    for i=1:length(excsigrois(:,1))
    sorted_excplotsmth = excsigrois(r,:);
    end
end

for j=1
    if isempty(excsigrois)
        break
    end
    figure;imagesc(sorted_excplotsmth),colorbar;title('Significantly Excited/Sorted ROIs');colormap(gca,nawhimar_auto);
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

goodspans_sup = find(spanLengthneg_mat>=30); %%change for span length
supsigrois = allplot_odor(goodspans_sup,:); %change 'allplot_odor' to 'allplot_smth' for smooth responses
%excroisz = zdata_odorsmth(goodspans,:)
[~,r] = sortrows(supsigrois,'ascend');
%figure;imagesc(excrois),colorbar;%%excitatory ROIs

%% sort suppressed ROIs by onset time after odor presentation
sorted_supplotsmth = [];
for i = 1:length(supsigrois(:,1))
    if isempty(supsigrois)
    break
    end
    for i=1:length(supsigrois(:,1))
    sorted_supplotsmth = supsigrois(r,:);
    end
end

for j=1
    if isempty(supsigrois)
        break
    end
    figure;imagesc(sorted_supplotsmth),colorbar;title('Significantly Suppressed/Sorted ROIs');colormap(gca,nawhimar_auto);
    xlabel('Frames')
    ylabel('ROI Number')%df/f allplot of all ROIs
end

%% check if goodspans_exc and goodspans_sup are overlapping
for i=1
    overlapspans = intersect(goodspans_exc,goodspans_sup);
        if isempty(overlapspans)
            break
        end
        disp('Both excitatory and suppressive spans are overlapping!')
end

%% vertically cocantenate arrays
% sorted_excplotsmth = [];
% sorted_supplotsmth = [];
for j=1
    if ~exist('sorted_excplotsmth','var') 
    if ~exist('sorted_supplotsmth','var')
        break
    else 
sigsortapodor_smth=vertcat(sorted_excplotsmth,sorted_supplotsmth);
    end
    end
end

for j=1
    if ~exist('sorted_excplotsmth','var')
    if ~exist('sorted_supplotsmth','var')
        break
    else
    figure;imagesc(sigsortapodor_smth),colorbar;title('Significantly Excited/Suppressed Sorted ROIs (todoron:todoroff)');colormap(gca,nawhimar_auto);
    end
    end
end

%% normalize each exc and sup significantly sorted ROI matrix separately
normexc=[];
for i=1
    if ~exist('sorted_excplotsmth','var') || isempty(sorted_excplotsmth)
        break
    end
for j=1:size(sorted_excplotsmth,1) %% in progress
    normexc(j,:) = normalised_diff(sorted_excplotsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normsup = [];
for i=1
    if ~exist('sorted_supplotsmth','var') || isempty(sorted_supplotsmth)
        break
    end
for j=1:size(sorted_supplotsmth,1) %% in progress
    normsup(j,:) = normalised_diff(sorted_supplotsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

%% now normalize each excitatory and suppressive matrix together
normall = vertcat(normexc,normsup);
figure;plot(normexc','r');hold on;plot(normsup','b'),title(['Significantly Excited and Suppressed ROIs Normalized and Sorted (Traces) {\color{red}Odor ID:}', num2str(odornum)])
%figure;imagesc(normall);colorbar;title(['Significantly Excited and Suppressed ROIs Normalized and Sorted {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
%c = colorbar;
%c.Label.String = 'Norm. \DeltaF/F';
%gca;xlabel({'Frames'});ylabel({'ROIs'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% create imagesc fig
figure;imagesc(normexc);colorbar;title(['Significantly Excited ROIs Normalized and Sorted{\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
figure;imagesc(normsup);colorbar;title(['Significantly Suppressed ROIs Normalized and Sorted{\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
figure;imagesc(normall);colorbar;title(['Significantly Excited and Suppressed ROIs Normalized and Sorted{\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});

%% save variables
prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');

for i=1
    if ~exist('supsigrois','var') || isempty(supsigrois)
        save ('excandsuproisunsorted.mat', 'excsigrois')
     elseif tf==1 && ~isempty(supsigrois)
        save ('excandsuproisunsorted.mat', 'excsigrois', 'supsigrois')
    end
end

for i=1
    if ~exist('supsigrois','var') || isempty(supsigrois)
        save ('excandsuprois_sorted.mat', 'sorted_excplotsmth')
    elseif tf==1 && ~isempty(supsigrois)
        save ('excandsuprois_sorted.mat', 'sorted_excplotsmth', 'sorted_supplotsmth')
    end
end
      
for i=1
    if ~exist('supsigrois','var') || isempty(supsigrois)
        save ('excandsuprois_sortandnorm.mat', 'normexc')
    elseif tf==1 && ~isempty(supsigrois)
        save ('excandsuprois_sortandnorm.mat', 'normexc','normsup')
    end
end

for i=1
    if ~exist('supsigrois','var') || isempty(supsigrois)
        save ('excandsuprois.mat', 'goodspans_exc')
    elseif tf==1 && ~isempty(supsigrois)
        save ('excandsuprois.mat', 'goodspans_exc', 'goodspans_sup')
    end
end
        
%% close all figures
prompt = ('Do you want to close all figures? Y/N?')
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf == 1
    close all
else
end

%clearvars -except excrois suprois 