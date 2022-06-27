disp('Select odorant folder with Ch2 Odorant .mat files from Odor Bank')
selpathch2=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');
cd(selpathch2)

disp(selpathch2)

prompt = 'What is the Odor Number?';
OdorNumID = input(prompt);
save('OdorNumID.mat', 'OdorNumID')