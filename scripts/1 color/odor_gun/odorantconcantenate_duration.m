clc
clearvars
close all
%% cocantenate across odorants
prompt = ('Insert Odorant Letters and Name:');
odorname = input(prompt,'s')


prompt = ('Does tbt51 odorant file exist?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
disp('Load FOV Cocantenation across a Single Odorant from tbt51')
cd('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR\tbt51_indiv_files\tbt51_struct')
uiopen('*.mat');
tbt51_allfovrois = allFOVodortraces;
tbt51_excfovrois = excFOVodortraces;
tbt51_supfovrois = supFOVodortraces;
    else
tbt51_allfovrois = [];
tbt51_excfovrois = [];
tbt51_supfovrois = [];  
    end
end

prompt = ('Does tbt57 odorant file exist?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
disp('Load FOV Cocantenation across a Single Odorant from tbt57')
cd('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR\tbt57_indiv_files\tbt57_struct')
uiopen('*.mat');
tbt57_allfovrois = allFOVodortraces;
tbt57_excfovrois = excFOVodortraces;
tbt57_supfovrois = supFOVodortraces;
    else
tbt57_allfovrois = [];
tbt57_excfovrois = [];
tbt57_supfovrois = [];  
    end
end

prompt = ('Does tbt61 odorant file exist?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
disp('Load FOV Cocantenation across a Single Odorant from tbt61')
cd('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR\tbt61_indiv_files\tbt61 struct')
uiopen('*.mat');
tbt61_allfovrois = allFOVodortraces;
tbt61_excfovrois = excFOVodortraces;
tbt61_supfovrois = supFOVodortraces;
    else
tbt61_allfovrois = [];
tbt61_excfovrois = [];
tbt61_supfovrois = [];  
    end
end

excglusnfrrois = vertcat(tbt51_excfovrois(:,(2:end)),tbt57_excfovrois(:,(2:end)),tbt61_excfovrois(:,(2:end)));
supglusnfrrois = vertcat(tbt51_supfovrois(:,(2:end)),tbt57_supfovrois(:,(2:end)),tbt61_supfovrois(:,(2:end)));
totglusnfrrois = vertcat(excglusnfrrois,supglusnfrrois);

excroinums = vertcat(tbt51_excfovrois,tbt57_excfovrois,tbt61_excfovrois);
excroinums = excroinums(:,1);
suproinums = vertcat(tbt51_supfovrois,tbt57_supfovrois,tbt61_supfovrois);
for i=1
    if isempty(supglusnfrrois)
        break
    else
suproinums = suproinums(:,1);
    end
end
totroinums = vertcat(excroinums,suproinums);

%% plot all glusnfr rois across experiments with labeled ROIs
figure;imagesc(totglusnfrrois);title(['Significant Excitatory/Suppressive Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
for j=1:size(totglusnfrrois)
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YTick = (1:j);
ax.YTickLabel = {totroinums(:,1)};
end
%% plot normalized all glusnfr rois across experiments with labeled ROIs
todoron=450;%%times when odor is on and off (all matricies are from 3:5 seconds CHANGE IF NEEDED
todoroff=750;

normexcglusnfrrois=[];
for j=1:size(excglusnfrrois,1) %% in progress
    normexcglusnfrrois(j,:) = normalised_diff(excglusnfrrois(j,1:todoron:todoroff:end)); %normalize to odor only
end

normsupglusnfrrois = [];
for i=1
    if ~exist('supglusnfrrois','var') || isempty(supglusnfrrois)
        break
    end
for j=1:size(supglusnfrrois,1) %% in progress
    normsupglusnfrrois(j,:) = normalised_diff(supglusnfrrois(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normallglusnfrrois = vertcat(normexcglusnfrrois,normsupglusnfrrois);
figure;imagesc(normallglusnfrrois);title(['Significant Excitatory/Suppressive Normalized Glomerular ROIs {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
for j=1:size(totglusnfrrois)
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YTick = (1:j);
ax.YTickLabel = {totroinums(:,1)};
end
%% plot exc glusnfr rois normalized & sorted to peak  %%include ROI number?
% for k = 1:size(excglusnfrrois,1) %%%%UNCOMMENT FOR SORTED BUT NOT NORMALIZED RESPONSES!
%    xpeaks(k) = max(excglusnfrrois(k,todoron:todoroff));
% end
% for k=1:size(excglusnfrrois,1)
%     peaktimes(k)=find(excglusnfrrois(k,:)==xpeaks(k));
% end
% peaktimes = peaktimes';
% [r,~] = sortrows(excglusnfrrois,peaktimes,'ascend');
% 
% figure;imagesc(r);title(['Significant Excitatory/Suppressive Sorted Glomerular ROIs {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
% c = colorbar;
% c.Label.String = '\DeltaF/F';
% gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% patch([450 450 750 750],[1 1 1 1], gray);
xfwhm=[];
fwmind=[];
x=1:1502;
times=(x)./150;

close all


for k = 1:size(normexcglusnfrrois,1)
   xfwhm(k) = fwhm(times,normexcglusnfrrois(k,:));
end

% for k = 1:size(normexcglusnfrrois,1)
%    xfwhm_odor(k) = fwhm(times,normexcglusnfrrois(k,450:750));
% end

[sorteddurtimes,fwhmind] = sort(xfwhm,'ascend');
normsorte = normexcglusnfrrois(fwhmind,:);

% [sorteddurtimes_odor,fwhmind_odor] = sort(xfwhm_odor,'ascend');
% normsorte_odor = normexcglusnfrrois(fwhmind_odor,:);

figure;imagesc(normsorte);title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted across time series {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
clims = [-1 1];
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
% figure;bar(xfwhm(fwhmind),1:length(xfwhm),'color','k')

% figure;imagesc(normsorte_odor);title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted by duration during odor{\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
% c = colorbar;
% c.Label.String = '\DeltaF/F';
% clims = [-1 1];
% gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% patch([450 450 750 750],[1 1 1 1], gray);
% % figure;bar(xfwhm(fwhmind),1:length(xfwhm),'color','k')
%% plot sup glusnfr rois normalized & sorted to peak
xfwhm=[];
fwmind=[];
for i=1
    if isempty(normsupglusnfrrois)
        normsorts = [];
    else 
        for k = 1:size(normsupglusnfrrois,1)
            xfwhm(k) = fwhm(times,normsupglusnfrrois(k,:));
        end
    end
[sorteddurtimes,fwhmind] = sort(xfwhm,'descend');
normsorts = normsupglusnfrrois(fwhmind,:);

figure;imagesc(normsorts);title(['Significant Suppressive Norm & Sorted Glomerular ROIs {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
%figure;bar(peaktimes(peakind),1:length(peaktimes),'color','k')
end
%% plot all norm and sort to peak
normsortes = vertcat(normsorte,normsorts);
figure;imagesc(normsortes,clims);title(['Significant Excitatory/Suppressive Norm & Sorted Glomerular ROIs {\color{red}Odor ID:}', num2str(odorname)]);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
hax = gca;
times=(0:1501)./150;
for i = [1 3 5 7 9]
xtimes(i)=find(times==i);
end
hax.XTick = [151,450,750,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

%%save unsorted responses for PCA analysis
prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
        save('allodormat', 'excroinums','excglusnfrrois','suproinums','supglusnfrrois')
    end
end