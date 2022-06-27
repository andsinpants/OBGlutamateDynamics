clc
clearvars
close all
%% load files
disp('Select all Odorant .mat Files from Odor Bank')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR');
cd(selpath)

file = selpath; %this is just getting the name automatically from the folder
[filepath,name,ext] = fileparts(file);

f = dir2(); %check if a directory
for i=1:8 % change to number of files, since f is including non-mat files
% for i =1:length(f)
    tf(i) = isdir(f(i).name);
    if tf(i) == 0
    f(i) = [];
    end
end

% for i = 1:length(f)
for i=1:8
   fnames = f(i).name;
   cd(fnames)
   fnamesload(i) = load('acrossFOVodortraces_unsorted.mat');
   cd(selpath);
end

allexcglurois = vertcat(fnamesload.excFOVodortraces);
allsupglurois = vertcat(fnamesload.supFOVodortraces);


prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
        save('allodorantsinbank_unsortedwithOdorNumID.mat', 'allexcglurois','allsupglurois')
    end
end