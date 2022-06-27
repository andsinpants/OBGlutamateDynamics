disp('Select odorant folder with Ch1 Odorant .mat files from Odor Bank')
selpathch1=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');
cd(selpathch1)

disp(selpathch1)

prompt = 'What is the Odor Number?';
OdorNumID = input(prompt);
save('OdorNumID.mat', 'OdorNumID')