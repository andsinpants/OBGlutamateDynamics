clearvars
clc
clear all
%close all
%% now load odornumID matrices
load('allodorantsinbank_catmatricesOdorNumID.mat') %%make sure to be in folder with experimental matrix
totglusnfr_bothbanks = vertcat(allexcglurois_bothbanks,allsupglurois_bothbanks);
sorttotglusnfr_bothbanks = sortrows(totglusnfr_bothbanks,1);
expname = 'tbt57';
%% plot all glusnfr rois across experiments with labeled ROIs
figure;imagesc(sorttotglusnfr_bothbanks);title(['Significant Excitatory/Suppressive Non-normalized/Unsorted Glomerular ROIs {\color{red}Exp ID:}', num2str(expname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);

%% plot normalized all glusnfr rois across experiments with labeled ROIs
todoron=451;%%times when odor is on and off (all matricies are from 3:5 seconds CHANGE IF NEEDED)
todoroff=751;

normallglurois_bothbanks=[];
for j=1:size(sorttotglusnfr_bothbanks,1) %% in progress
    normallglurois_bothbanks(j,:) = normalised_diff(sorttotglusnfr_bothbanks(j,3:todoron:todoroff:end)); %normalize to odor only
end

% normallexcglurois_bothbanks = normallexcglurois_bothbanks(:,3:end);
% normallsupglurois_bothbanks = normallsupglurois_bothbanks(:,3:end);

figure;imagesc(normallglurois_bothbanks);title(['Significant Excitatory/Suppressive Normalized/Sorted Glomerular ROIs {\color{red}Exp ID:}', num2str(expname)]);ylabel('Glom-Odor Pairs based on Odor ID');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm .\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
%% set XTick Labels
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
%% set YTick Labels
hax = gca;
nyticks = [1:size(sorttotglusnfr_bothbanks)];
nyticks = nyticks';
nytickslabel = sorttotglusnfr_bothbanks(:,1);
yticklabs = unique(nytickslabel);

for k=1:length(yticklabs)
    spanLength = ismember(nytickslabel,yticklabs(k));
    spanLength_range(k) = sum(spanLength>0);
end

spanLength_range = spanLength_range';
IDandspanLength = horzcat(yticklabs,spanLength_range);

tickstarts = cumsum(spanLength_range);

for g=1:length(tickstarts)
tickstarts(g) = tickstarts(g)+1;
end

tickstarts(end)=[];
tickstarts = vertcat(1,tickstarts);

hax.YTick= [tickstarts];
hax.YTickLabel = [yticklabs];