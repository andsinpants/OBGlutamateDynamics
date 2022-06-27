clearvars -except totplot_odorch1hz totplot_odorch1uf totplot_odorch21hz totplot_odorch2uf myplot_times

%% load unfiltered data for display

%% organize responses (may need editing?)
roiinsert = (1:39)';
roich1 = repmat(roiinsert,12,1);
roich2 = repmat(roiinsert,12,1);

totplot_odorch1 = horzcat(roich1,totplot_odorch1uf);
totplot_odorch2 = horzcat(roich2,totplot_odorch2uf);

prompt = 'Odor on?';
odoron = input(prompt)
prompt = 'Odor off?';
odoroff = input(prompt)
timemat_odor = myplot_times;
timemat_odoron = find(timemat_odor(1,:)>=odoron);
timemat_odoroff = find(timemat_odor(1,:)<=odoroff);
todoron = (timemat_odoron(1)+2);
todoroff = (timemat_odoroff(end)+2);
tpreodor = 1:(timemat_odoron(1)+2);

totsize = size(totplot_odorch1,1)*2;
allmat = zeros(totsize,length(totplot_odorch1));

OddI = 1:2:size(allmat,1);
EvenI = 2:2:size(allmat,1);

for j=1:size(totplot_odorch1,1)
    allmat(OddI(j),:) = totplot_odorch1(j,:);
end

for j=1:size(totplot_odorch2,1)
    allmat(EvenI(j),:) = totplot_odorch2(j,:);
end

figure;imagesc(allmat);colorbar;title('All \DeltaF/F waterfall plot of GG-O Pairs with Ch1 and Ch2 1Hz filtered');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%% now normalized to its own max and min
%% normalize all responses of Ch1 and Ch2
normallmat = [];
for i=1
   for j=1:size(allmat,1) %% in progress
    normallmat(j,:) = normalised_diff(allmat(j,1:todoron:todoroff:end)); %normalize to odor only
   end
end


figure;imagesc(normallmat);colorbar;title('All \DeltaF/F norm waterfall of GG-O Pairs with Ch1 and Ch2 1Hz filtered');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%% now find significant responses only in Ch1 and plot those ROIs together%%
totch1odor = totplot_odorch1;
totch2odor = totplot_odorch2;

%% find mean of odor responses (ODOR) ONLY CH1!!!
for g = 1:size(totch1odor)
    meanvals_odorch1(g) = mean(totch1odor(g,todoron:todoroff),2);
    meanvals_preodorch1(g) = mean(totch1odor(g,tpreodor)); 
end
meanvals_odorch1 = meanvals_odorch1';
meanvals_preodorch1 = meanvals_preodorch1';

%% find 95_max of odor responses (ODOR) ONLY CH1!!!
for g = 1:size(totch1odor)
    maxvals_odorch1(g) = max(totch1odor(g,todoron:todoroff));
end
maxvals_odorch1 = maxvals_odorch1';
max95_odorch1 = .95*(maxvals_odorch1);

%% find 5pct of odor responses (ODOR) ONLY CH1!!!
for g = 1:size(totch1odor)
qprctile_odorch1(g,:) = qprctile(totch1odor(g,todoron:todoroff), [5.0 95.0]);
end
pct5_odorch1 = qprctile_odorch1(:,1);

%% find std of odor responses (ODOR) ONLY CH1!!!
for g = 1:size(totch1odor)
std_preodorch1(g,:) = std(totch1odor(g,1:todoron));
end
%close all%remove after writing script
%% zscore data
for g= 1:size(totch1odor)
    zdata_odorch1(g,:) = ((totch1odor(g,:)-meanvals_preodorch1(g))./std_preodorch1(g));
end
 %% smooth responses display ONLY CH1!!!
  allplot_smthch1 = [];
for i=1:length(totch1odor(:,1))
    ys=totch1odor(i,:);
    smoothbch1=fastsmooth(ys,15,3,1);
    allplot_smthch1=[allplot_smthch1;smoothbch1];
end

%% zscore smoothed data ONLY CH1!!!
for g= 1:size(allplot_smthch1,1)
    zdata_odorsmthch1(g,:) = ((allplot_smthch1(g,:)-meanvals_preodorch1(g))./std_preodorch1(g));
end
figure;plot(zdata_odorsmthch1')

%% get std descriptive stats ONLY CH1!!!
for g=1:size(zdata_odorsmthch1,1)
b_odorch1(g,:) = zdata_odorsmthch1(g,todoron:todoroff);
b_odormaxch1(g) = max(b_odorch1(g,:));
b_odorminch1(g) = min(b_odorch1(g,:));
b_odormeanch1(g) = mean(b_odorch1(g,:));
b_odorpct5ch1(g,:) = qprctile(zdata_odorsmthch1(g,todoron:todoroff), [5.0 95.0]);
%b_odorpct5ch1(g) = b_odorpct5ch1(:,1)
end
b_odormaxch1 = b_odormaxch1';
b_odorminch1 = b_odorminch1';
b_odorpct5ch1 = b_odorpct5ch1(:,1);
b_odormeanch1 = b_odormeanch1';
b_odormax95ch1 = b_odormaxch1(:,:)*.95;

%% find significant zscored data by responses above 4SD  ONLY CH1!!!
aboveposThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmthch1(:,1))
    posThreshold = 4; %%change for threshold
    aboveposThreshold(i,:) = (zdata_odorsmthch1(i,todoron:todoroff))>= posThreshold; 
end

for i=1:length(zdata_odorsmthch1(:,1))
    spanLengthpos = regionprops(aboveposThreshold(i,:));
    if isempty(spanLengthpos)
     spanLengthpos_mat(i) = 0;
    else
    spanLengthpos = spanLengthpos.Area;
    spanLengthpos_mat(i) = spanLengthpos;
    end
end

goodspans_excch1 = find(spanLengthpos_mat>=30); %%change for span length
excsigroisch1 = allplot_smthch1(goodspans_excch1,:);

%% find significant zscored data by responses below 4SD ONLY CH1!!!
belownegThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmthch1(:,1))
    negThreshold = -4; %%change for threshold
    belownegThreshold(i,:) = (zdata_odorsmthch1(i,todoron:todoroff))<= negThreshold; 
end

for i=1:length(zdata_odorsmthch1(:,1))
    spanLengthneg = regionprops(belownegThreshold(i,:));
    if isempty(spanLengthneg)
     spanLengthneg_mat(i) = 0;
    else
    spanLengthneg = spanLengthneg.Area;
    spanLengthneg_mat(i) = spanLengthneg;
    end
end

goodspans_supch1 = find(spanLengthneg_mat>=30); %%change for span length
for i=1
    tf = isempty(goodspans_supch1);
    if tf ==0;
    supsigroisch1 = allplot_smthch1(goodspans_supch1,:);
    else
        break
    end
end


%% find mean of odor responses (ODOR) ONLY CH2!!!   <--------start of ch2 significance responses
for g = 1:size(totch2odor)
    meanvals_odorch2(g) = mean(totch2odor(g,todoron:todoroff),2);
    meanvals_preodorch2(g) = mean(totch2odor(g,tpreodor)); 
end
meanvals_odorch2 = meanvals_odorch2';
meanvals_preodorch2 = meanvals_preodorch1';

%% find 95_max of odor responses (ODOR) ONLY CH2!!!
for g = 1:size(totch2odor)
    maxvals_odorch2(g) = max(totch2odor(g,todoron:todoroff));
end
maxvals_odorch2 = maxvals_odorch2';
max95_odorch2 = .95*(maxvals_odorch2);

%% find 5pct of odor responses (ODOR) ONLY CH2!!!
for g = 1:size(totch2odor)
qprctile_odorch2(g,:) = qprctile(totch2odor(g,todoron:todoroff), [5.0 95.0]);
end
pct5_odorch2 = qprctile_odorch2(:,1);

%% find std of odor responses (ODOR) ONLY CH2!!!
for g = 1:size(totch2odor);
std_preodorch2(g,:) = std(totch2odor(g,1:todoron));
end
%close all%remove after writing script
%% zscore data
for g= 1:size(totch2odor)
    zdata_odorch2(g,:) = ((totch2odor(g,:)-meanvals_preodorch2(g))./std_preodorch2(g));
end
 %% smooth responses display ONLY CH2!!!
  allplot_smthch2 = [];
for i=1:length(totch2odor(:,1))
    ys=totch2odor(i,:);
    smoothbch2=fastsmooth(ys,15,3,1);
    allplot_smthch2=[allplot_smthch2;smoothbch2];
end

%% zscore smoothed data ONLY CH2!!!
for g= 1:size(allplot_smthch2,1)
    zdata_odorsmthch2(g,:) = ((allplot_smthch2(g,:)-meanvals_preodorch2(g))./std_preodorch2(g));
end
figure;plot(zdata_odorsmthch2')

%% get std descriptive stats ONLY CH2!!!
for g=1:size(zdata_odorsmthch2,1)
b_odorch2(g,:) = zdata_odorsmthch2(g,todoron:todoroff);
b_odormaxch2(g) = max(b_odorch2(g,:));
b_odorminch2(g) = min(b_odorch2(g,:));
b_odormeanch2(g) = mean(b_odorch2(g,:));
b_odorpct5ch2(g,:) = qprctile(zdata_odorsmthch2(g,todoron:todoroff), [5.0 95.0]);
%b_odorpct5ch1(g) = b_odorpct5ch1(:,1)
end
b_odormaxch2 = b_odormaxch2';
b_odorminch2 = b_odorminch2';
b_odorpct5ch2 = b_odorpct5ch2(:,1);
b_odormeanch2 = b_odormeanch2';
b_odormax95ch2 = b_odormaxch2(:,:)*.95;

%% find significant zscored data by responses above 4SD  ONLY CH2!!!
aboveposThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmthch2(:,1))
    posThreshold = 4; %%change for threshold
    aboveposThreshold(i,:) = (zdata_odorsmthch2(i,todoron:todoroff))>= posThreshold; 
end

for i=1:length(zdata_odorsmthch2(:,1))
    spanLengthpos = regionprops(aboveposThreshold(i,:));
    if isempty(spanLengthpos)
     spanLengthpos_mat(i) = 0;
    else
    spanLengthpos = spanLengthpos.Area;
    spanLengthpos_mat(i) = spanLengthpos;
    end
end

goodspans_excch2 = find(spanLengthpos_mat>=30); %%change for span length
excsigroisch2 = allplot_smthch2(goodspans_excch2,:);

%% find significant zscored data by responses below 4SD ONLY CH2!!!
belownegThreshold = [];
spanLengthpos_mat =[];

for i=1:length(zdata_odorsmthch2(:,1))
    negThreshold = -4; %%change for threshold
    belownegThreshold(i,:) = (zdata_odorsmthch2(i,todoron:todoroff))<= negThreshold; 
end

for i=1:length(zdata_odorsmthch2(:,1))
    spanLengthneg = regionprops(belownegThreshold(i,:));
    if isempty(spanLengthneg)
     spanLengthneg_mat(i) = 0;
    else
    spanLengthneg = spanLengthneg.Area;
    spanLengthneg_mat(i) = spanLengthneg;
    end
end

goodspans_supch2 = find(spanLengthneg_mat>=30); %%change for span length
for i=1
    tf = isempty(goodspans_supch2);
    if tf ==0;
    supsigroisch2 = allplot_smthch2(goodspans_supch2,:);
    else
        break
    end
end


%% for the entire dataset find intersecting ROIs and mark with 1, non intesecting ROIs with 0
clearvars -except todoron todoroff allplotsmthch1 allplotsmthch2 goodspans_excch1 goodspans_excch2 goodspans_supch1 goodspans_supch2 excsigroisch1 excsigroisch2 supsigroisch1 supsigroisch2 totplot_odorch1 totplot_odorch2 myplot_times

allsigrois = horzcat(goodspans_excch1,goodspans_excch2,goodspans_supch1,goodspans_supch2);
[C,ia,ib] = unique(allsigrois);

for i=1:length(C)
    sigroisch1(i,:) = totplot_odorch1(C(i),:);
    sigroisch2(i,:) = totplot_odorch2(C(i),:);
    figure;plot(sigroisch1(i,3:end),'g');hold on; plot(sigroisch2(i,3:end),'m');title(['Row:',num2str(C(i)),'-ROI',num2str(sigroisch1(i,1)),'-Odor',num2str(sigroisch1(i,2))]);
    %saveas(gcf,sprintf('tbt68_%d.png',i));
end

close all
totsize = length(C)*2;
allmat = zeros(totsize,length(totplot_odorch1));
betweenmat = zeros(1,length(totplot_odorch1));

OddI = 1:2:size(allmat,1);
EvenI = 2:2:size(allmat,1);

for j=1:size(sigroisch1,1)
    allmat(OddI(j),:) = sigroisch1(j,:);
end

for j=1:size(sigroisch2,1)
    allmat(EvenI(j),:) = sigroisch2(j,:);
end

for i=2:size(OddI,2)
   akmtest(i,:) = allmat(OddI(i),:);betweenmat; allmat(EvenI(i),:);
end

figure;imagesc(allmat(:,3:end));colorbar;title('All \DeltaF/F heatmap of Significant Glom-Odor Pairs with Ch1 and Ch2 combined');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%saveas(gcf,'tbt68_signonnormalizedheatmaps_OdorBankX.pdf') %change for odor bank specificity

%% normalize allmat of both channels significantly 
normallmat = [];
for i=1
   for j=1:size(allmat,1) %% in progress
    normallmat(j,:) = normalised_diff(allmat(j,3:todoron:todoroff:end)); %normalize to odor only
   end
end

figure;imagesc(normallmat);colorbar;title('All \DeltaF/F normalized heatmap of Significant Glom-Odor Pairs with Ch1 and Ch2 combined');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs
gray = [.9 .9 .9];
patch([451 451 1051 1051],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

%saveas(gcf,'tbt68_signormalizedheatmaps_OdorBankX.pdf')

prompt = ('Do you want to save vars? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf==1
    save( 'tbt68_ch1and2vars.mat', 'allmat', 'normallmat', 'sigroisch1', 'sigroisch2', 'goodspans_excch1', 'goodspans_excch2', 'goodspans_supch1', 'goodspans_supch2', 'excsigroisch1', 'excsigroisch2', 'supsigroisch1', 'supsigroisch2');
else
end







