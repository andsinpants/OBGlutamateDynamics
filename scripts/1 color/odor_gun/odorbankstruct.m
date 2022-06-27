clc
clearvars
close all

%% chose Odor Bank folder
disp('Choose Odor Bank folder to make into struct and place OdorNumID into')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR');
cd(selpath)

%% eliminate non-folders from struct
files = dir2;
for i = 1:length(files)
fnames(i) = files(i);
tf(i) = isdir(fnames(i).name);
if tf(i) == 0
    fnames(i) = [];
end
end

%% now place all into struct 
% for i = 1:length(fnames)
for i=1:10
    cd(fnames(i).name)
    load('acrossFOVodortraces_unsorted.mat')
    data.OdorNumID = OdorNumID;
    data.allFOVodortraces = allFOVodortraces;
    data.excFOVodortraces = excFOVodortraces;
    data.supFOVodortraces = supFOVodortraces;
    prompt = (['Do you want to save ', num2str(fnames(i).name) ' struct?  Y/N?']);
    savequestion = input(prompt,'s');
    tf=strcmpi(savequestion,'Y');
        for i=1
            if tf==1
            save('acrossFOVodortraces_unsorted_struct.mat', 'data')
            end
        end
   cd(selpath)     
end

clearvars allFOVodortraces excFOVodortraces supFOVodortraces

%% now place OdorNumID into matrix
% for i = 1:length(fnames)
for i=1:10
    cd(fnames(i).name)
    load('acrossFOVodortraces_unsorted.mat')
    j = 1:size(allFOVodortraces);
    repOdorNumID =repelem(OdorNumID,size(j,2))';
    allFOVodortraces = horzcat(repOdorNumID,allFOVodortraces);
clearvars repOdorNumID
    k = 1:size(excFOVodortraces);
    repOdorNumID =repelem(OdorNumID,size(k,2))';
    excFOVodortraces = horzcat(repOdorNumID,excFOVodortraces);
clearvars repOdorNumID
    l = 1:size(supFOVodortraces,1);
    n = isempty(l);
        if n ==1
        else
        repOdorNumID = repelem(OdorNumID,size(l,2))';
        supFOVodortraces = horzcat(repOdorNumID,supFOVodortraces);
        end
    prompt = (['Do you want to save ', num2str(fnames(i).name) ' matrix?  Y/N?']);
    savequestion = input(prompt,'s');
    tf=strcmpi(savequestion,'Y');
        for i=1
            if tf==1
                save('acrossFOVodortraces_unsorted.mat', 'allFOVodortraces','excFOVodortraces','supFOVodortraces','OdorNumID')
            end
        end
    cd(selpath)
clearvars repOdorNumID
end