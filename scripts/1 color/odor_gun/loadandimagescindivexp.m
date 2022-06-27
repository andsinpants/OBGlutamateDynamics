clc
clearvars
close all
%% load files
disp('Select all Odorant .mat Files from Odor Bank')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR');
cd(selpath)

file = selpath; %this is just getting the name automatically from the folder
[filepath,name,ext] = fileparts(file);

f = dir('*.mat');
for i =1:length(f)
   fnames=f.name;
end

load(fnames);

figure;imagesc(excFOVodortraces(:,2:end));title(['Significant Excitatory, Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(name)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);


figure;imagesc(supFOVodortraces(:,2:end));title(['Significant Suppressive, Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(name)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);


figure;imagesc(allFOVodortraces(:,2:end));title(['Significant Excitatory & Suppressive Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(name)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);

prompt = ('Do you want to save the last figure? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
pathshift = 'C:\Users\matt\Desktop';
for i=1
    if tf==1
        cd(pathshift) 
        saveas(gcf,[num2str(name)],'pdf')
    end
end