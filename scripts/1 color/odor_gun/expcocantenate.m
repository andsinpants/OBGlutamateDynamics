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
   fnames=f(i).name;
   fnamesload(i) = load(fnames);
end

allexcglurois = vertcat(fnamesload.excglusnfrrois);
allexcgluroinums = vertcat(fnamesload.excroinums);
allsupglurois = vertcat(fnamesload.supglusnfrrois);
allsupgluroinums = vertcat(fnamesload.suproinums);

