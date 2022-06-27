clc
clearvars
close all
%% select folder with odorgun files
disp('Select Odorant Folder from Odor Gun Data')
selpath=uigetdir('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SF-iGluSnFR');
cd(selpath)

file = selpath; %this is just getting the name automatically from the folder
[filepath,name,ext] = fileparts(file);
oname = name;
%% check to see if folder exists, then cocantenate across fields of view
for i=1
    fold = isfolder('FOV1'); 
        if fold==1
            cd FOV1
            load('excandsuprois.mat')
            goodspans_exc = goodspans_exc';
            pluralroi = isfile('excandsuproisunsorted.mat');
                if pluralroi==1 
                load('excandsuproisunsorted.mat') 
                else
                load('excandsuproiunsorted.mat');
                end   
            FOV1_excdfoverftraces = horzcat(goodspans_exc,excsigrois);
            filex = exist('goodspans_sup');
                if filex == 1
                    goodspans_sup = goodspans_sup';
                    FOV1_supdfoverftraces = horzcat(goodspans_sup,supsigrois);
                else 
                    FOV1_supdfoverftraces = [];
                end
        else 
            FOV1_excdfoverftraces = [];
            FOV1_supdfoverftraces = [];
        end
end

cd(selpath)
clearvars fold filex pluralroi goodspans_exc excsigrois goodspans_sup supsigrois

for i=1
    fold = isfolder('FOV2'); 
        if fold==1
            cd FOV2
            load('excandsuprois.mat')
            goodspans_exc = goodspans_exc';
            pluralroi = isfile('excandsuproisunsorted.mat');
                if pluralroi==1 
                load('excandsuproisunsorted.mat') 
                else
                load('excandsuproiunsorted.mat');
                end   
            FOV2_excdfoverftraces = horzcat(goodspans_exc,excsigrois);
            filex = exist('goodspans_sup');
                if filex==1
                    goodspans_sup = goodspans_sup';
                    FOV2_supdfoverftraces = horzcat(goodspans_sup,supsigrois);
                else 
                    FOV2_supdfoverftraces = [];
                end
        else 
            FOV2_excdfoverftraces = [];
            FOV2_supdfoverftraces = [];
        end
end

cd(selpath)
clearvars fold filex pluralroi goodspans_exc excsigrois goodspans_sup supsigrois

for i=1
    fold = isfolder('FOV3'); 
        if fold==1
            cd FOV3
            load('excandsuprois.mat')
            goodspans_exc = goodspans_exc';
            pluralroi = isfile('excandsuproisunsorted.mat');
                if pluralroi==1 
                load('excandsuproisunsorted.mat') 
                else
                load('excandsuproiunsorted.mat');
                end   
            FOV3_excdfoverftraces = horzcat(goodspans_exc,excsigrois);
            filex = exist('goodspans_sup');
                if filex==1
                    goodspans_sup = goodspans_sup';
                    FOV3_supdfoverftraces = horzcat(goodspans_sup,supsigrois);
                else 
                    FOV3_supdfoverftraces = [];
                end
        else
            FOV3_excdfoverftraces = [];
            FOV3_supdfoverftraces = [];
        end
end
cd(selpath)
clearvars goodspans_exc goodspans_sup i selpath supsigrois excsigrois filex fold pluralroi
%% now concantenate across all FOVs for excitatory ROIs

excFOVodortraces = vertcat(FOV1_excdfoverftraces,FOV2_excdfoverftraces,FOV3_excdfoverftraces);
supFOVodortraces = vertcat(FOV1_supdfoverftraces,FOV2_supdfoverftraces,FOV3_supdfoverftraces);

allFOVodortraces = vertcat(excFOVodortraces,supFOVodortraces);

%% find replicate ROIs that are double-labeled as excitatory/suppressive
for i=1
    if ~isempty(supFOVodortraces) 
        la=ismember(allFOVodortraces(:,1),supFOVodortraces(:,1));

        roirepnum = find(la==1)'
        roicontent= allFOVodortraces(la)'
        roinumandcont = vertcat(roirepnum,roicontent)
        fliproirepnum = fliplr(roirepnum);

            for i=1 
                if numel(roirepnum)~=1;
                    roioverlap=[];%% check for overlapping spans in comparison to all
                    for h=1:length(roirepnum)
                     for j=1:size(supFOVodortraces)
                        allroioverlap(h,:) = allFOVodortraces(roirepnum(h),:);
                        suproioverlap(j,:)= allFOVodortraces(fliproirepnum(j),:);
                        roioverlap(j) = isequal(allFOVodortraces(roirepnum(j),:),allFOVodortraces(fliproirepnum(j),:))
                     end        
                    end
                    for i=1
                        if any(roioverlap==1)
                        allFOVodortraces(fliproirepnum(j),:)=[];
                        supFOVodortraces=(flipud(supFOVodortraces));
                        supFOVodortraces((j),:)=[];
                        else
                            break
                        end
                    end
                end
            end
    else
        break
    end
end

%% now plot figures
figure;imagesc(excFOVodortraces(:,2:end));title(['Significant Excitatory, Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(oname)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 850 850],[1 1 1 1], gray);


figure;imagesc(supFOVodortraces(:,2:end));title(['Significant Suppressive, Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(oname)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 850 850],[1 1 1 1], gray);


figure;imagesc(allFOVodortraces(:,2:end));title(['Significant Excitatory & Suppressive Unsorted Glomerular ROIs {\color{red}Odor ID:}', num2str(oname)]);ylabel('Glom-Odor Pairs');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 850 850],[1 1 1 1], gray);

clearvars -except excFOVodortraces supFOVodortraces allFOVodortraces
%% now save arrays
prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
        save('acrossFOVodortraces_unsorted.mat', 'allFOVodortraces','excFOVodortraces','supFOVodortraces')
    end
end


