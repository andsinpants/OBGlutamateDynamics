clc
close all
clear all

disp('Select folder with Ch1 Odorant .mat files from Odor Bank')
selpathch1=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');

disp('Select folder with Ch2 Odorant .mat files from Odor Bank')
selpathch2=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');

fpathch1 = dir2(selpathch1);
fpathch2 = dir2(selpathch2);%this is just getting the name automatically from the folder

numnames_ch1 = length(fpathch1);
numnames_ch2 = length(fpathch2);
if numnames_ch1 > numnamesch1
   


for i=1:length(
A = [[fpathch1(:).name];[fpathch2(:).name]]';
B = [[b(:).x];[fpathch2(:).name]]';
[C, ia, ib] = intersect( A, B, 'rows' );
c = a(ia);

for i = 1
cd(fpathch1(i).name);
akmtest_ch1(i) = load('acrossFOVodortraces_unsorted_ch1.mat');
cd(fpathch2(i).name);
akmtest_ch2(i) = load('acrossFOVodortraces_unsorted_ch2.mat');
end


[C,ia,ib] = intersect(akmtest_ch1(:,2),akmtest_ch2(:,2));

%%unfinished
