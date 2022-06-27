load('allodormat_bothodorbanks_noGOpair445.mat')
allexcglurois_bothbanks = allexcglurois;
allsupglurois_bothbanks= allsupglurois;
totglusnfr_bothbanks = vertcat(allexcglurois_bothbanks,allsupglurois_bothbanks);
% expname = 'tbt61';
%% plot all glusnfr rois across experiments with labeled ROIs
figure;imagesc(totglusnfr_bothbanks);title(['Significant Excitatory/Suppressive Unsorted Glomerular ROIs {\color{red}Odor ID:}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);

%% plot normalized all glusnfr rois across experiments with labeled ROIs
todoron=451;%%times when odor is on and off (all matricies are from 3:5 seconds CHANGE IF NEEDED
todoroff=751;

normallexcglurois_bothbanks=[];
for j=1:size(allexcglurois_bothbanks,1) %% in progress
    normallexcglurois_bothbanks(j,:) = normalised_diff(allexcglurois_bothbanks(j,1:todoron:todoroff:end)); %normalize to odor only
end

normallsupglurois_bothbanks = [];
for i=1
    if ~exist('allsupglurois_bothbanks','var') || isempty(allsupglurois_bothbanks)
        break
    end
for j=1:size(allsupglurois_bothbanks,1) %% in progress
    normallsupglurois_bothbanks(j,:) = normalised_diff(allsupglurois_bothbanks(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normallglusnfrrois = vertcat(normallexcglurois_bothbanks,normallsupglurois_bothbanks);
figure;imagesc(normallglusnfrrois);title(['Significant Excitatory/Suppressive Normalized Glomerular ROIs {\color{red}Odor ID:}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
%% plot exc glusnfr rois normalized & sorted to peak  %%include ROI number?
% for k = 1:size(allexcglurois_bothbanks,1) %%%%UNCOMMENT FOR SORTED BUT NOT NORMALIZED RESPONSES!
%    xpeaks(k) = max(allexcglurois_bothbanks(k,todoron:todoroff));
% end
% for k=1:size(allexcglurois_bothbanks,1)
%     peaktimes(k)=find(allexcglurois_bothbanks(k,:)==xpeaks(k));
% end
% peaktimes = peaktimes';
% [r,~] = sortrows(allexcglurois_bothbanks,peaktimes,'ascend');
% 
% figure;imagesc(r);title(['Significant Excitatory/Suppressive Sorted Glomerular ROIs {\color{red}Odor ID:}', num2str(expname)]);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
% c = colorbar;
% c.Label.String = '\DeltaF/F';
% gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% patch([451 451 751 751],[1 1 1 1], gray);

x=1:1502;
times=(x)./150;

xpeaks=[];
for k = 1:size(normallexcglurois_bothbanks,1)
   xpeaks(k) = max(normallexcglurois_bothbanks(k,todoron:todoroff));
end
for k=1:size(normallexcglurois_bothbanks,1)
    peaktimes(k)=find(normallexcglurois_bothbanks(k,:)==xpeaks(k));
end

[sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
normsorte = normallexcglurois_bothbanks(peakind,:);

figure;imagesc(normsorte);title(['Significant Excitatory Norm & Sorted Glomerular ROIs {\color{red}Odor ID:}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
hold on;
plot(peaktimes(peakind),1:length(peaktimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
%% plot sup glusnfr rois normalized & sorted to trough
xtroughs=[];
troughtimes=[];
for i=1
    if isempty(normallsupglurois_bothbanks)
        normsorts = [];
        break
    else
%         xpeaks=[];
        for k = 1:size(normallsupglurois_bothbanks,1)
           xtroughs(k,:) = min(normallsupglurois_bothbanks(k,todoron:todoroff));
        end
%         peaktimes=[];
        for k=1:size(normallsupglurois_bothbanks,1)
            troughtimes(k)=find(normallsupglurois_bothbanks(k,:)==xtroughs(k));
        end
    end


[sortedtroughtimes,troughind] = sort(troughtimes,'descend');
normsorts = normallsupglurois_bothbanks(troughind,:);

figure;imagesc(normsorts);title(['Significant Suppressive Norm & Sorted Glomerular ROIs {\color{red}Odor ID:}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
hold on;
plot(troughtimes(troughind),1:length(troughtimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
end
%% plot all norm and sort to peak
normsortes = vertcat(normsorte,normsorts);
figure;imagesc(normsortes);title(['Significant Excitatory/Suppressive Norm & Sorted Glomerular ROIs {\color{red}Odor ID:}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
patch([751 751 1502 1502],[79 79 79 79], gray);
patch([751 751 1502 1502],[158 158 158 158], gray);
patch([751 751 1502 1502],[237 237 237 237], gray);
patch([751 751 1502 1502],[316 316 316 316], gray);
patch([751 751 1502 1502],[395 395 395 395], gray);
patch([751 751 1502 1502],[474 474 474 474], gray);
patch([751 751 1502 1502],[553 553 553 553], gray);
patch([751 751 1502 1502],[632 632 632 632], gray);
patch([751 751 1502 1502],[711 711 711 711], gray);
patch([751 751 1502 1502],[712 712 712 712], gray);
patch([751 751 1502 1502],[718 718 718 718], gray);
hax = gca;
peaktimes = times(peaktimes);
troughtimes = times(troughtimes);
times=(0:1501)./150;
for i = [1 3 5 7 9]
xtimes(i)=find(times==i);
end
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};