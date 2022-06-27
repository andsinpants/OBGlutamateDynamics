close all; 
clear;

[prefile,prepath] = uigetfile('*.mat','Choose Predrug File(s) (2Hz filter?)','MultiSelect','off');
if ~prefile; return; end
pre = load(fullfile(prepath,prefile)); clear prepath prefile;
[postfile,postpath] = uigetfile('*.mat','Choose Postdrug File(s) (2Hz filter?)','MultiSelect','off');
if ~postfile; return; end
post = load(fullfile(postpath,postfile)); clear postpath postfile;


bplotbyodor = 1; %1 means figure# same as odor#; 0 means figure# is roi#
% ibase = 1:275; %preodor frames
% iodor = 301:726; %postodor frames 2sec long
ibase = 1:275; %preodor frames
iodor = 301:825; %postodor frames 2sec long

minthresh = -7; maxthresh = 7; %first 6 then 7

odornums = [];
for r = 1:length(pre.myplotdata.avgfile.roi)
    for o = 1:length(pre.myplotdata.avgfile.roi(r).odor)
        odornums = union(odornums,pre.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
    end
end
results = nan(length(pre.myplotdata.avgfile.roi),odornums(end),2);
% results2 = nan(length(ex.myplotdata.avgfile.roi)*24,4); %% use for odor gun with 24 odorants (2 banks)
results2 = nan(length(pre.myplotdata.avgfile.roi)*12,4); %% use for odor gun with 12 odorants (1 bank)

figure(256); hold on; %odor number vs significant min/max values
title('Channel #1');
prebasemeanall = [];
postbasemeanall = [];
index=0;
for r = 1:length(pre.myplotdata.avgfile.roi)
     for o = 1:length(pre.myplotdata.avgfile.roi(r).odor)
        odornum = pre.myplotdata.avgfile.roi(r).odor(o).number;
        
        prebasemean = mean(pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        prebasemeansub = pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-prebasemean;
        prebasemeanall = vertcat(prebasemeansub,prebasemeanall);
        
        postbasemean = mean(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        postbasemeansub = post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-postbasemean;
        postbasemeanall = vertcat(postbasemeansub,postbasemeanall);

     end
     
     
    prebasestd = std(prebasemeanall);
    postbasestd = std(postbasemeanall);
    
    for o = 1:length(pre.myplotdata.avgfile.roi(r).odor)
        odornum = pre.myplotdata.avgfile.roi(r).odor(o).number;
        prebasemean = mean(pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        prebasemeano = mean(pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(iodor));
        prept1 = max(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(iodor(1:150)));
        prept2 = max(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(750:825));
        prepdiff = prept2-prept1;
        %prebasestd = std(pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        prezts = (pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-prebasemean)./prebasestd;
        tmpmax = prctile(prezts(301:601),95);
        postbasemean = mean(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        postbasemeano = mean(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(iodor));
        postpt1 = max(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(iodor(1:150)));
        postpt2 = max(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(750:825));
        postpdiff = postpt2-postpt1;
        %postbasestd = std(post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        postzts = (post.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-postbasemean)./postbasestd;
%         tmpmin = prctile(postzts(iodor),15); %default is 25
        
        
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(0), ', zvalmax' num2str(tmpmax)];   
        index=index+1;

        if tmpmax>maxthresh %only excited
             prepdiff=results(r,odornum,1);
             results_trace(index,1:length(prezts))=prezts;  %<<<<<<<<<<< extract time series (exc)
             results_ts(index,2)= odornum; 
             results_ts(index,1)=r;
             results_ts(index,3)= prepdiff; 
             results_ts(index,4)=0;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(pre.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,prezts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','g');
        elseif tmpmax <= maxthresh
            results2(r,odornum,1)=prepdiff;
            results_trace2(index,1:length(prezts))=postzts;
            results_ts2(index,2)=odornum;
            results_ts2(index,1)=r;
            results_ts2(index,3)=postpdiff; 
            results_ts2(index,4)=0;
        end
    end
end

results_all=horzcat(results_ts, results_trace);
results_all2=horzcat(results_ts2, results_trace2);
%clearvars -except results_all
autoArrangeFigures();

%% now remove arrays with all 0s
results_sig = results_all(any(results_all,2),:);

results_pre = results_sig;
results_postraw = results_all2;

[uovrlp1,ia1,ib1] = union(results_sig(:,1),results_all2(:,1));
[uovrlp2,ia2,ib2] = union(results_sig(:,2),results_all2(:,2));


times = (0:size(results_sig(1,5:end),2)-1)/150; 
%% quick plot of figure
figure;imagesc(results_sig);colormap(nawhimar_auto);
% figure;imagesc(results2new);colormap(nawhimar_auto)
clearvars -except results_pre results_postraw
close all
%% extract exc, sup, biphasic responses and normalize
biph_sige = results_sig(:,3)~=0;
biph_sigs = results_sig(:,4)~=0;
biph_sig = horzcat(biph_sige, biph_sigs);
ind = find(biph_sig(:,1) & biph_sig(:,2)==1);
biph_all = results_sig(ind,:);

for j=1:size(biph_all,1) %% in progress
    if isempty(biph_all)
        normb_all=[];
    else
    normbiph_all(j,:) = normalised_diff(biph_all(j,5:iodor(1)+5:iodor(end)+5:end));
    end
end

if isempty(biph_all)
    normb_all=[];
else
normb_all = horzcat(biph_all(:,1),biph_all(:,2),biph_all(:,3),biph_all(:,4),normbiph_all); %normalize to odor only
end

exc_sig = results_sig(:,3)~=0 & results_sig(:,4)==0;
ind = find(exc_sig==1);
exc_all = results_sig(ind,:);

for j=1:size(exc_all,1) %% in progress
    if isempty(exc_all)
        normexc_all = [];
    else
    normexc_all(j,:) = normalised_diff(exc_all(j,5:iodor(1)+5:iodor(end)+5:end));
    end
end

if isempty(exc_all)
    norme_all=[];
else
norme_all = horzcat(exc_all(:,1),exc_all(:,2),exc_all(:,3),exc_all(:,4),normexc_all); %normalize to odor only
end
%norme_all = vertcat(norme_all(1:8,:),norme_all(10:end,:)); %%% use for tbt61_fov2

for i =1
    if isempty(biph_all)
    break
else
norme_all = vertcat(normb_all,norme_all);
end
end

sup_sig = results_sig(:,4)~=0 & results_sig(:,3)==0;
ind = find(sup_sig==1);
sup_all = results_sig(ind,:);

for j=1:size(sup_all,1) %% in progress
    if isempty(sup_all)
        normsup_all = [];
    else
    normsup_all(j,:) = normalised_diff(sup_all(j,5:iodor(1)+5:iodor(end)+5:end));
    end
end

if isempty(sup_all)
    norms_all=[];
else
norms_all = horzcat(sup_all(:,1),sup_all(:,2),sup_all(:,3),sup_all(:,4),normsup_all); %normalize to odor only
end


%% sort by odor latency to peak
% only excitatory
xpeaks=[];
for k = 1:size(norme_all,1)
   xpeaks(k) = max(norme_all(k,iodor(1)+5:iodor(end)+5));
end
for k=1:size(norme_all,1)
    peaktimes(k)=find(norme_all(k,5:end)==xpeaks(k));
end

[sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
normsorte_l = norme_all(peakind,:);

figure;imagesc(normsorte_l(:,5:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 726 726],[1 1 1 1], gray);
hold on;
plot(peaktimes(peakind),1:length(peaktimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};

%% plot sup glusnfr rois normalized & sorted to trough
% now suppressive
xtroughs=[];
troughtimes=[];
for i=1
    if isempty(norms_all)
        normsorts = [];
        break
    else
%         xpeaks=[];
        for k = 1:size(norms_all,1)
           xtroughs(k,:) = min(norms_all(k,iodor(1)+5:iodor(end)+5));
        end
%         peaktimes=[];
        for k=1:size(norms_all,1)
            troughtimes(k)=find(norms_all(k,5:end)==xtroughs(k));
        end
    end


[sortedtroughtimes,troughind] = sort(troughtimes,'descend');
normsorts_l = norms_all(troughind,:);

figure;imagesc(normsorts_l(:,5:end));title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 726 726],[1 1 1 1], gray);
hold on;
plot(troughtimes(troughind),1:length(troughtimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};
end

%% sort by odor duration
% only excitatory
for k = 1:size(norme_all,1)
   xfwhm(k) = fwhm(times,norme_all(k,5:end));
end

[sorteddurtimes,fwhmind] = sort(xfwhm,'ascend');
normsorte_d = norme_all(fwhmind,:);

figure;imagesc(normsorte_d(:,5:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
clims = [-1 1];
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};

% now suppressive
xfwhm=[];
fwmind=[];
for i=1
    if isempty(norms_all)
        norms_all = [];
    else 
        for k = 1:size(norms_all,1)
            xfwhm(k) = fwhm(times,norms_all(k,5:end));
        end
    end
[sorteddurtimes,fwhmind] = sort(xfwhm,'descend');
normsorts_d = norms_all(fwhmind,:);

figure;imagesc(normsorts_d(:,5:end));title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = '\DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};
%figure;bar(peaktimes(peakind),1:length(peaktimes),'color','k')
end

% clearvars -except normsorte_l normsorte_d normsorts_l normsorts_d results_sig times ext1 ext2 supt1 supt2
% clearvars -except results_sig times 
%% save all sig responses
% prompt = ('Do you want to save? Y/N?');
% savequestion = input(prompt,'s');
% tf=strcmpi(savequestion,'Y');
% 
% for i=1
%     if tf == 1
%        save ('results_timeseries.mat', 'results_sig', 'normsorte_l','normsorts_l','normsorte_d','normsorts_d','times')
%     else
%         break
%     end
% end

%% ask to save
% prompt = ('Do you want to save? Y/N?');
% savequestion = input(prompt,'s');
% tf=strcmpi(savequestion,'Y');
% 
% for i=1
%     if ~exist('newresultszmat_ex','var') || isempty(newresultszmat_sup)
%         save ('excandsuproisunsorted.mat', 'excsigrois')
%      elseif tf==1 && ~isempty(newsresultszmat_sup)
%         save ('excandsuproisunsorted.mat', 'excsigrois', 'supsigrois')
%     end
% end


% %% now excitatory only
% for r = 1:length(ex.myplotdata.file.roi)
%     for o = 1:length(ex.myplotdata.file.roi(r).odor)
%         allts(r*o,:) = ex.myplotdata.file.roi(r).odor(o).avgtrial.series(:);
%     end
% end
% 
% for i = 1:length(ex.myplotdata.file.roi)
%     odor1 = repmat(odornums(1),i,1);
%     odor2 = repmat(odornums(2),i,1);
%     odor3 = repmat(odornums(3),i,1);
%     odor4 = repmat(odornums(4),i,1);
%     odor5 = repmat(odornums(5),i,1);
%     odor6 = repmat(odornums(6),i,1);
%     odor7 = repmat(odornums(7),i,1);
%     odor8 = repmat(odornums(8),i,1);
%     odor9 = repmat(odornums(9),i,1);
%     odor10 = repmat(odornums(10),i,1);
%     odor11 = repmat(odornums(11),i,1);
%     odor12 = repmat(odornums(12),i,1);
% end
% 
% allodors = vertcat(odor1,odor2,odor3,odor4,odor5,odor6,odor7,odor8,odor9,odor10,odor11,odor12);
% rois = repmat(1:length(ex.myplotdata.file.roi),1,12)';
% 
% allts_randoexc = horzcat(rois,allodors,allts);
% 
% allts_randoexcsig = allts_randoexc;
% excthrowaway = find(all(allts_randoexc(:,3:end)==0,2));
% allts_randoexcsig(excthrowaway,:) = [];



%% now suppressive only
% for r = 1:length(sup.myplotdata.file.roi)
%     for o = 1:length(sup.myplotdata.file.roi(r).odor)
%         allts(r*o,:) = sup.myplotdata.file.roi(r).odor(o).avgtrial.series(:);
%     end
% end
% 
% for i = 1:length(sup.myplotdata.file.roi)
%     odor1 = repmat(odornums(1),i,1);
%     odor2 = repmat(odornums(2),i,1);
%     odor3 = repmat(odornums(3),i,1);
%     odor4 = repmat(odornums(4),i,1);
%     odor5 = repmat(odornums(5),i,1);
%     odor6 = repmat(odornums(6),i,1);
%     odor7 = repmat(odornums(7),i,1);
%     odor8 = repmat(odornums(8),i,1);
%     odor9 = repmat(odornums(9),i,1);
%     odor10 = repmat(odornums(10),i,1);
%     odor11 = repmat(odornums(11),i,1);
%     odor12 = repmat(odornums(12),i,1);
% end
% 
% allodors = vertcat(odor1,odor2,odor3,odor4,odor5,odor6,odor7,odor8,odor9,odor10,odor11,odor12);
% rois = repmat(1:length(sup.myplotdata.file.roi),1,12)';
% 
% allts_randosup = horzcat(rois,allodors,allts);
%clear ex sup


% % % hplots = findobj('type','axes');
% % % for p = 1:numel(hplots)
% % %     axes(hplots(p)); ylim([-5 10]);
% % % end
% % exall=0; exsame=0; supall = 0; supsame = 0;
% % for r = 1:length(ch{1}.exMax_95.roi)
% % for o = 1:length(ch{1}.exMax_95.roi(r).odor)
% % if ch{1}.exMax_95.roi(r).odor(o)>maxthresh || ch{2}.exMax_95.roi(r).odor(o)>maxthresh; exall = exall+1; end %excited in 1 ch
% % if ch{1}.supMin_5.roi(r).odor(o)<minthresh || ch{2}.supMin_5.roi(r).odor(o)<minthresh; supall = supall+1; end %suppressed in 1 ch
% % if ch{1}.exMax_95.roi(r).odor(o)>maxthresh && ch{2}.exMax_95.roi(r).odor(o)>maxthresh; exsame = exsame+1; end %excited in both ch
% % if ch{1}.supMin_5.roi(r).odor(o)<minthresh && ch{2}.supMin_5.roi(r).odor(o)<minthresh; supsame = supsame+1; end %suppressed in both ch
% % end
% % end