cd('L:\DropboxInstall\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SFiGluSnFR_jRGECO1a\tbt_jRGECO_SFiGluSnFR ogun only\tbt_jRGECO_SFiGluSnFR_ogun analysis');
clearvars
close all
load('2color_tbt132g_137f_sigtraces.mat');

figure;imagesc(sortedgtraces);title(['{\color{green}Ch1} Significant Excitatory Norm & Sorted Glomerular ROIs by sum iodor']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
cmin = -1;
cmax = 1;
caxis([cmin cmax]);
c.Label.String = 'norm. z-scored \DeltaF/F';
% gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651,1951,2251];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11';'13';'15'};

figure;imagesc(sortedrtraces);title(['{\color{red}Ch2} Significant Excitatory Norm & Sorted Glomerular ROIs by sum iodor']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
cmin = -1;
cmax = 1;
caxis([cmin cmax]);
c.Label.String = 'norm. z-scored \DeltaF/F';
% gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651,1951,2251];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11';'13';'15'};

