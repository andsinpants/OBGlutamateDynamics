clearvars
clc
%% choose folder and load in variables
disp('Select Folder')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\pcd iGluSnFR 3x odors');
cd(selpath);
prompt = 'Odorant (use initials)?';
odornum = input(prompt,'s');
%% now manually load in variables
disp('Load in variables')
vars = uigetfile('*.mat*','MultiSelect','on');

for i=1:length(vars)
    load(vars{i})
end
%% check for missing sup files
for i=1
    if ~exist('normexc1x','var')
        normexc1x=[];
    end
end

for i=1
    if ~exist('normsup1x','var')
        normsup1x=[];
    end
end

for i=1
    if ~exist('normexc3x','var')
        normexc3x=[];
    end
end

for i=1
    if ~exist('normsup3x','var')
        normsup3x=[];
    end
end

for i=1
    if ~exist('normexc10x','var')
        normexc10x=[];
    end
end

for i=1
    if ~exist('normsup10x','var')
        normsup10x=[];
    end
end

%whitespc = zeros(50,length(normexc10x));        

prompt = ('Ready to run? Y/N?');
runinput = input(prompt,'s');
tf=strcmpi(runinput,'Y');

for i=1
if tf==0
    break
else
normall1x = vertcat(normexc1x,normsup1x);%%load in 1x Odo
figure; subplot(1,3,1),imagesc(normall1x);title(['1x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
gca;xlabel({'Frames'});ylabel({'ROIs'});
end
end

for i=1
if tf==0
    break
else
normall3x = vertcat(normexc3x,normsup3x);%%load in 1x Odor 1
hold on; subplot(1,3,2),imagesc(normall3x);title(['3x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);
gca;xlabel({'Frames'});ylabel({'ROIs'});
end
end

for i=1
if tf==0
    break
else
normall10x = vertcat(normexc10x,normsup10x);%%load in 1x Odor 1
hold on; subplot(1,3,3),imagesc(normall10x);title(['10x {\color{red}Odor ID:}', num2str(odornum)]);colormap(gca,nawhimar_auto);%colorbar
% c = colorbar;
% c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});
end
end


%% new odor?
% prompt = ('New odor? Y/N?');
% runinput = input(prompt,'s');
% tf_odor=strcmpi(runinput,'Y');
