clc
clearvars
close all
%% chose folder with ch1 odor gun files
disp('Select folder with Ch1 Odorant .mat files from Odor Bank')
selpathch1=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');
cd(selpathch1)

fpathch1 = dir2(selpathch1); %this is just getting the name automatically from the folder
for i = 1:length(fpathch1)
cd(fpathch1(i).name);


for k=1
    load('excandsuprois_ch1.mat');
    tk = exist('goodspans_sup');
        if tk ~= 0
        load('excandsuproisunsorted_ch1.mat')
        else
        load('excandsuproiunsorted_ch1.mat');
        end
    load('OdorNumID.mat');
    goodspans_exc=goodspans_exc';
    OdorNumIDrep = repelem(OdorNumID,size(excsigrois,1))';
    excFOVodortraces_ch1 = horzcat(OdorNumIDrep,goodspans_exc,excsigrois);
end

for k=1
    ch = exist('supsigrois');
    if ch ~= 0
        goodspans_sup=goodspans_sup';
        OdorNumIDrep = repelem(OdorNumID,size(supsigrois,1))';
        supFOVodortraces_ch1 = horzcat(OdorNumIDrep,goodspans_sup,supsigrois);
    else
        supFOVodortraces_ch1 = [];
    end
end

allFOVodortraces_ch1 = vertcat(excFOVodortraces_ch1,supFOVodortraces_ch1);%vertically concatenate both

prompt = (['Do you want to save ', num2str(fpathch1(i).name) ' matrix?  Y/N?']);
    savequestion = input(prompt,'s');
    tf=strcmpi(savequestion,'Y');
        for i=1
            if tf==1
                save('acrossFOVodortraces_unsorted_ch1.mat', 'allFOVodortraces_ch1','excFOVodortraces_ch1','supFOVodortraces_ch1','OdorNumID')
            end
        end
     clearvars goodspans_exc goodspans_sup excsigrois supsigrois
    cd(selpathch1)
end
%% chose folder with Ch2 odor gun files
clearvars 
close all

disp('Select folder with Ch2 Odorant .mat files from Odor Bank')
selpathch2=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR and jRGECO1a');
cd(selpathch2)

fpathch2 = dir2(selpathch2); %this is just getting the name automatically from the folder
for i = 1:length(fpathch2)
cd(fpathch2(i).name);


for k=1
    load('excandsuprois_ch2.mat');
    tk = exist('goodspans_sup');
        if tk ~= 0
        load('excandsuproisunsorted_ch2.mat')
        else
        load('excandsuproiunsorted_ch2.mat');
        end
    load('OdorNumID.mat');
    goodspans_exc=goodspans_exc';
    OdorNumIDrep = repelem(OdorNumID,size(excsigrois,1))';
    excFOVodortraces_ch2 = horzcat(OdorNumIDrep,goodspans_exc,excsigrois);
end

for k=1
    ch = exist('supsigrois');
    if ch ~= 0
        goodspans_sup=goodspans_sup';
        OdorNumIDrep = repelem(OdorNumID,size(supsigrois,1))';
        supFOVodortraces_ch2 = horzcat(OdorNumIDrep,goodspans_sup,supsigrois);
    else
        supFOVodortraces_ch2 = [];
    end
end

allFOVodortraces_ch2 = vertcat(excFOVodortraces_ch2,supFOVodortraces_ch2);%vertically concatenate both

prompt = (['Do you want to save ', num2str(fpathch2(i).name) ' matrix?  Y/N?']);
    savequestion = input(prompt,'s');
    tf=strcmpi(savequestion,'Y');
        for i=1
            if tf==1
                save('acrossFOVodortraces_unsorted_ch2.mat', 'allFOVodortraces_ch2','excFOVodortraces_ch2','supFOVodortraces_ch2','OdorNumID')
            end
        end
     clearvars goodspans_exc goodspans_sup excsigrois supsigrois
    cd(selpathch2)
end