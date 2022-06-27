%%plotconcwaterfalls_clut2b2.m
%close all
clearvars
clc
%% choose folder and load in variables
disp('Select Folder')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\pcd iGluSnFR 3x odors');
cd(selpath);
prompt = 'Odorant (use initials)?';
odornum = input(prompt,'s');
%% now manually load in variables
disp('Load in variables')
vars = uigetfile('*.mat*','MultiSelect','on');

for i=1:length(vars)
    load(vars{i})
end

%load('excandsuproisunsorted_10x.mat')
%% check for missing sup files
for i=1
    if ~exist('normexc1x','var')
        normexc1x=[];
    end
end

for i=1
    if ~exist('normsup1x','var')
        normsup1x=[];
    end
end

for i=1
    if ~exist('normexc3x','var')
        normexc3x=[];
    end
end

for i=1
    if ~exist('normsup3x','var')
        normsup3x=[];
    end
end

for i=1
    if ~exist('normexc10x','var')
        normexc10x=[];
    end
end

for i=1
    if ~exist('normsup10x','var')
        normsup10x=[];
    end
end

for i=1
    if ~exist('excsigrois10x','var')
        excsigrois10x=[];
    end
end

for i=1
    if ~exist('supsigrois10x','var')
        supsigrois10x=[];
    end
end

for i=1
    if ~exist('excsigrois3x','var')
        excsigrois3x=[];
    end
end

for i=1
    if ~exist('supsigrois3x','var')
        supsigrois3x=[];
    end
end

for i=1
    if ~exist('goodspans_exc1x','var')
        goodspans_exc1x=[];
    end
end

for i=1
    if ~exist('goodspans_sup1x','var')
        goodspans_sup1x=[];
    end
end

for i=1
    if ~exist('excrois_1xsmth','var')
        excrois_1xsmth=[];
    end
end

for i=1
    if ~exist('suprois_1xsmth','var')
        suprois_1xsmth=[];
    end
end

for i=1
    if ~exist('excrois_3xsmth','var')
        excrois_3xsmth=[];
    end
end

for i=1
    if ~exist('suprois_3xsmth','var')
        suprois_3xsmth=[];
    end
end

for i=1
    if ~exist('excrois_10xsmth','var')
        excrois_10xsmth=[];
    end
end

for i=1
    if ~exist('suprois_10xsmth','var')
        suprois_10xsmth=[];
    end
end

for i=1
    if ~exist('rawall1x','var')
        rawall1x=[];
    end
end

for i=1
    if ~exist('rawall3x','var')
        rawall3x=[];
    end
end

for i=1
    if ~exist('rawall10x','var')
        rawall10x=[];
    end
end

rawall1x = vertcat(excrois_1xsmth,suprois_1xsmth);
rawall3x = vertcat(excrois_3xsmth,suprois_3xsmth);
rawall10x = vertcat(excrois_10xsmth,suprois_10xsmth);
allspans = vertcat(goodspans_exc1x',goodspans_sup1x');
allspans1x = horzcat(allspans,rawall1x);
allspans3x = horzcat(allspans,rawall3x);
allspans10x = horzcat(allspans,rawall10x);


% excrois_1xsmth = excsigrois1x;
% suprois_1xsmth = supsigrois1x;
%whitespc = zeros(50,length(normexc10x));  
load('timemat_odor_conc.mat') %% load in time series data for normalization
for i = 1:length(myplotdata.file.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.file.roi(i);
end

for j = 1:length(myplot_r) 
allplot_odor(j,:) = myplot_r(j).odor.avgtrial.series';
end

myplot_times = myplot_r(1).odor.avgtrial.time;
% allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented
%prompt = 'Odor on?';
odoron = 3;
%prompt = 'Odor off?';
odoroff = 7;
timemat_odoron = find(timemat_odor(1,:)>=odoron);
timemat_odoroff = find(timemat_odor(1,:)<=odoroff);
todoron = timemat_odoron(1);
todoroff = timemat_odoroff(end);
tpreodor = 1:timemat_odoron(1);

% renormalize for heatmaps
normall1x=[];
for i=1:length(rawall1x(:,1))
    if ~exist('rawall1x','var') || isempty(rawall1x)
        break
    end
 for j=1:size(rawall1x,1) %% in progress
    normall1x(j,:) = normalised_diff(rawall1x(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normall3x=[];
for i=1:length(rawall3x(:,1))
    if ~exist('rawall3x','var') || isempty(rawall3x)
        break
    end
for j=1:size(rawall3x,1) %% in progress
    normall3x(j,:) = normalised_diff(rawall3x(j,1:todoron:todoroff:end)); %normalize to odor only
end
end


normall10x=[];
for i=1:length(rawall10x(:,1))
    if ~exist('rawall10x','var') || isempty(rawall10x)
        break
    end
for j=1:size(rawall10x,1) %% in progress
    normall10x(j,:) = normalised_diff(rawall10x(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

clearvars allplot_odor
% % for i=1
% %     if isempty(excrois_10xsmth)% || isempty(excsigrois10x)
% %         break
% %     else
% % for j=1:size(excrois_10xsmth,1) %% in progress
% %     normexc10x(j,:) = normalised_diff(excrois_10xsmth(j,1:end)); %normalize to odor only
% % end
% %     end
% % end
% % 
% % for i=1
% %     if isempty(suprois_10xsmth)%~exist('supsigrois10x','var') || isempty(supsigrois10x)
% %         break
% %     else
% % for j=1:size(suprois_10xsmth,1) %% in progress
% %     normsup10x(j,:) = normalised_diff(suprois_10xsmth(j,1:end)); %normalize to odor only
% % end
% %     end
% % end
 
prompt = ('Ready to run? Y/N?');
runinput = input(prompt,'s');
tf=strcmpi(runinput,'Y');

for i=1
if tf==0
    break
else
% normall1x = vertcat(normexc1x,normsup1x);%%load in 1x Odor
clims = [-.2 1];
figure; subplot(1,3,1),imagesc(normall1x,clims);title(['1x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,clut2b2),yticklabels('manual') ;
gca;xlabel({'Time (s)'});ylabel({'ROIs'});
gray = [.9 .9 .9];
patch([451 451 1051 1051],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
end
end

% goodspans_exc1x=[];
% goodspans_sup1x=[];
load('excandsuprois_1x.mat') % create tick marks for significant ROIs
for i=1
    if ~exist('goodspans_sup1x', 'var')
    break
else
allspans = vertcat(goodspans_exc1x',goodspans_sup1x');
for j=1:length(allspans)
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YTick = (1:j);
ax.YTickLabel = {allspans(1:j)};
end
    end
end


for i=1
if tf==0
    break
else
% normall3x = vertcat(normexc3x,normsup3x);%%load in 1x Odor 1
clims = [-.2 1];
hold on; subplot(1,3,2),imagesc(normall3x,clims);title(['3x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,clut2b2);
gca;xlabel({'Time (s)'});ylabel({'ROIs'});
gray = [.9 .9 .9];
patch([451 451 1051 1051],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
end
end

for i=1
    if ~exist('goodspans_sup1x', 'var')
    break
else
allspans = vertcat(goodspans_exc1x',goodspans_sup1x');
for j=1:length(allspans)
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YTick = (1:j);
ax.YTickLabel = {allspans(1:j)};
end
    end
end


for i=1
if tf==0
    break
else
% normall10x = vertcat(normexc10x,normsup10x);%%load in 1x Odor 1
clims = [-.2 1];
hold on; subplot(1,3,3),imagesc(normall10x,clims);title(['10x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,clut2b2);%colorbar
% c = colorbar;
% c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Time (s)'});ylabel({'ROIs'});
gray = [.9 .9 .9];
patch([451 451 1051 1051],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
end
end

for i=1
    if ~exist('goodspans_sup1x', 'var')
    break
else
allspans = vertcat(goodspans_exc1x',goodspans_sup1x');
for j=1:length(allspans)
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YTick = (1:j);
ax.YTickLabel = {allspans(1:j)};
end
    end
end




% excrois_1xsmth = [];
% excrois_3xsmth = [];
% excsigrois10x = [];
% suprois_1xsmth = [];
% suprois_3xsmth = [];
% supsigrois10x = [];
%% choose ROI timeseries for selected or interesting odors 
prompt = ('Do you want to look at particular ROIs? Y/N?'); %%%%%%%SOMETHING HERE IS adding an extra row on
runinput = input(prompt,'s');
tf_roi=strcmpi(runinput,'Y');

% for i=1
%     if tf_roi==1
%         load('excandsuproisunsorted_1x.mat')
%     end
% end
% 
% for i=1
%     if tf_roi==1
%         load('excandsuproisunsorted_3x.mat')
%     end
% end
% 
% for i=1
%     if tf_roi==1
%         load('excandsuproisunsorted_10x.mat')
%     end
% end
% excrois_1xsmth = [];
% excrois_3xsmth = [];
% excsigrois10x = [];
% suprois_1xsmth = [];
% suprois_3xsmth = [];
% supsigrois10x = [];


answer = inputdlg('Enter ROI numbers',...
             'What ROIs do you want to look at?', [1 50]);
roiinput = str2num(answer{1}); %%%%%%%%%%%%%%%%%%%%%%NEED TO PLACE IN ROI CHOICE NUMBER BASED ON YTICKLABEL!

for k=1:(length(roiinput))
    n = roiinput(k);
    ind = find(n==allspans);
if tf_roi ~= 1
    break
else
    gray = [.9 .9 .9]; %odor onset and offset signal overlay 
    patch([451 451 1050 1051],[0 0 0 0], gray);
    figure;subplot(1,3,1);plot(rawall1x(ind,:));title(['1x {\color{red}Odor ID:}', num2str(odornum),'{\color{blue}ROI:}',num2str(roiinput(k))]);hold on
    gray = [.9 .9 .9]; %odor onset and offset signal overlay 
    patch([451 451 1050 1051],[0 0 0 0], gray);
    subplot(1,3,2);plot(rawall3x(ind,:));title(['3x {\color{red}Odor ID:}', num2str(odornum),'{\color{blue}ROI:}',num2str(roiinput(k))]);hold on
    gray = [.9 .9 .9]; %odor onset and offset signal overlay 
    patch([451 451 1050 1051],[0 0 0 0], gray);
    subplot(1,3,3);plot(rawall10x(ind,:));title(['10x {\color{red}Odor ID:}', num2str(odornum),'{\color{blue}ROI:}',num2str(roiinput(k))]);hold on
    gray = [.9 .9 .9]; %odor onset and offset signal overlay 
    patch([451 451 1050 1051],[0 0 0 0], gray);
end
end
autoArrangeFigures()
%% close all figures
prompt = ('Do you want to close all figures? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf == 1
    close all
else
end
    