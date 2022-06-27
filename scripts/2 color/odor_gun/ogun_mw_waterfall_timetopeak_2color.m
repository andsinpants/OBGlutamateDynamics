close all; 
clear;

%channel 1
[exfile,expath] = uigetfile('*.mat','Choose Channel #1 Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ch1_ex = load(fullfile(expath,exfile)); clear expath exfile;
[supfile,suppath] = uigetfile('*.mat','Choose Channel #1 Suppressive File(s) (0.5Hz filter?)','MultiSelect','off');
if ~supfile; return; end
ch1_sup = load(fullfile(suppath,supfile)); clear suppath supfile;

%channel 2
[exfile,expath] = uigetfile('*.mat','Choose Channel #2 Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; clear ibase iodor bplotbyodor; return; end
ch2_ex = load(fullfile(expath,exfile)); clear expath exfile;
[supfile,suppath] = uigetfile('*.mat','Choose Channel #2 Suppressive File(s) (0.5Hz filter?)','MultiSelect','off');
if ~supfile; clear ibase iodor bplotbyodor; return; end
ch2_sup = load(fullfile(suppath,supfile)); clear suppath supfile;


bplotbyodor = 1; %1 means figure# same as odor#; 0 means figure# is roi#
ibase = 1:275; %preodor frames
iodor = 301:825; %postodor frames
minthresh = -7; maxthresh = 7; %first 6 then 7

odornums = [];
for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
    for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornums = union(odornums,ch1_ex.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
    end
end

% odornums2 = [];
% for r = 1:length(ch2_ex.myplotdata.avgfile.roi)
%     for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
%         odornums2 = union(odornums2,ch2_ex.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
%     end
% end
results = nan(length(ch1_ex.myplotdata.avgfile.roi),odornums(end),2);
results2 = nan(length(ch1_ex.myplotdata.avgfile.roi)*24,7);
results(:,:,3:4) = nan(size(results,1),size(results,2),2);   %inceased to add one more column for 5 values %questionable responses
% results_trace = zeros(length(ch1_ex.myplotdata.avgfile.roi)*24,1502);
% results_trace2 = zeros(length(ch2_ex.myplotdata.avgfile.roi)*24,1502);
% results_ts = zeros(length(ch1_ex.myplotdata.avgfile.roi)*24,4);
% results_ts2 = zeros(length(ch2_ex.myplotdata.avgfile.roi)*24,4);

resultsexc1 = nan(length(ch1_ex.myplotdata.avgfile.roi),(odornums(end)),1502);
resultssup1 = nan(length(ch1_sup.myplotdata.avgfile.roi),(odornums(end)),1502);

resultsexc2 = nan(length(ch2_ex.myplotdata.avgfile.roi),(odornums(end)),1502);
resultssup2 = nan(length(ch2_sup.myplotdata.avgfile.roi),(odornums(end)),1502);

figure(256); hold on; %odor number vs significant min/max values
title('Channel #1');
figure(512); hold on; %odor number vs significant min/max values
title('Channel #2');
exbasemeanall = [];
supbasemeanall = [];
exbasemeanall2 = [];
supbasemeanall2 = [];
index=0;
% times = (0:size(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(:,6:end),1)-1)/150;

for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
     for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_ex.myplotdata.avgfile.roi(r).odor(o).number;
        
        exbasemean = mean(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
        
        supbasemean = mean(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
     
     %make sure these variables have unique names for ch2.
        exbasemean2 = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub2 = ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean2;
        exbasemeanall2 = vertcat(exbasemeansub2,exbasemeanall2);
        
        supbasemean2 = mean(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub2 = ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean2;
        supbasemeanall2 = vertcat(supbasemeansub2,supbasemeanall2);
     end
     
     
    exbasestd = std(exbasemeanall);
    supbasestd = std(supbasemeanall);
    
    exbasestd2 = std(exbasemeanall2);  
    supbasestd2 = std(supbasemeanall2);
    
    
    for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %exbasestd = std(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exzts = (ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean)./exbasestd;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %supbasestd = std(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean)./supbasestd;
        tmpmin = prctile(supzts(iodor),15);
        
%         odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
        %this is for the channel 2 calculations
        exbasemean2 = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %exbasestd2 = std(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exzts2 = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2)./exbasestd2;
        tmpmax2 = prctile(exzts2(iodor),95);
        supbasemean2 = mean(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %supbasestd2 = std(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts2 = (ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean2)./supbasestd2;
        tmpmin2 = prctile(supzts2(iodor),15);
        
        rgcorr=corrcoef(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean, ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2);
%         rgcorr_sup = corrcoef(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean, ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2);
        rgcorrcoef=rgcorr(1,2);
%         rgcorrcoeff_sup=rgcorr_sup(1,2);
               
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin),', zvalmax' num2str(tmpmax)];   
        tmplabel2 = ['Ch2, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin2),', zvalmax' num2str(tmpmax2)];
        index=index+1;
  
%         results2(index,7)= rgcorrcoeff;
%% channel 1 sorting into biphasic, excitatory, suppressive responses
        if tmpmax>maxthresh && tmpmin<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
%             results(r,odornum,1)=tmpmax; 
%             results(r,odornum,2)=tmpmin; %results(r,odornum,5)= rgcorrcoeff;
%             results2(index,3)=tmpmax;
%             results2(index,4)=tmpmin;
%             results2(index,2)=odornum;
%             results2(index,1)=r;
            results_trace(index,1:length(exzts))=exzts;   %<<<<<<<<<<< extract time series (biphasic)
            results_ts(index,1)=r; 
            results_ts(index,2)=odornum;
            results_ts(index,3)=tmpmax; 
            results_ts(index,4)=tmpmin; 
            results_ts(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on;
            else
                figure(r); hold on;
            end
            plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'r-','DisplayName',tmplabel);
            plot(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'m-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','r');
            plot(odornum,results(r,odornum,2),'x','Color','r');
        elseif tmpmax>maxthresh %only excited
%              results(r,odornum,1)=tmpmax; 
%              results(r,odornum,2)=0;%results(r,odornum,5)= rgcorrcoeff;
%              results2(index,3)=tmpmax; 
%              results2(index,2)=odornum;
%              results2(index,1)=r;
             results_trace(index,1:length(exzts))=exzts;   %<<<<<<<<<<< extract time series (biphasic)
             results_ts(index,1)=r; 
             results_ts(index,2)=odornum; 
             results_ts(index,3)=tmpmax; 
             results_ts(index,4)=0; 
             results_ts(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','g');
        elseif tmpmin<minthresh %only suppressed
%             results(r,odornum,1)=0; 
%             results(r,odornum,2)=tmpmin; %results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
%             results2(index,4)=tmpmin; 
%             results2(index,2)=odornum;
%             results2(index,1)=r;
            results_trace(index,1:length(supzts))=exzts; %<<<<<<<<<<<<< extract time series (sup)
            results_ts(index,1)=r; 
            results_ts(index,2)=odornum;
            results_ts(index,3)=0; 
            results_ts(index,4)=tmpmin;
            results_ts(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,2),'x','Color','b');
        else
             results_trace(index,1:length(exzts))=exzts;   %<<<<<<<<<<< extract time series (biphasic)
             results_ts(index,1)=r; 
             results_ts(index,2)=odornum; 
             results_ts(index,3)=0; 
             results_ts(index,4)=0;
             results_ts(index,5)=rgcorrcoef;
        end

        % this is for channel 2 if statements.... change labels to '2'
        % versions.2

%% channel 2 sorting into biphasic, excitatory, suppressive responses
       if tmpmax2>maxthresh && tmpmin2<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
%             results(r,odornum,3) = tmpmax2; 
%             results(r,odornum,4) = tmpmin2; 
%             results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
%             results2(index,3)=tmpmax2; 
%             results2(index,4)=tmpmin2; results2(index,2)= odornum; results2(index,1)=r;
            results_trace2(index,1:length(exzts2))= exzts2;   %<<<<<<<<<<< extract time series (biphasic)
            results_ts2(index,1)=r; 
            results_ts2(index,2)= odornum;
            results_ts2(index,3)=tmpmax2; 
            results_ts2(index,4)=tmpmin2;
            results_ts2(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on; %figure(256+odornum); hold on;  <-----use for separate ch2 plots
            else
                figure(r); hold on; %figure(256+r); hold on;  <-----use for separate ch2 plots
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts2,'r--','DisplayName',tmplabel2);
            plot(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts2,'m--','DisplayName',tmplabel2);
            figure(512);
            plot(odornum,results(r,odornum,3),'x','Color','r');
            plot(odornum,results(r,odornum,4),'x','Color','r');
        elseif tmpmax2>maxthresh %only excited
%              results(r,odornum,3)=tmpmax2; results(r,odornum,4)=0; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
%              results2(index,3)=tmpmax2;
%              results2(index,4)=0; 
%              results2(index,2)=odornum; 
%              results2(index,1)=r;
             results_trace2(index,1:length(exzts2))=exzts2;   %<<<<<<<<<<< extract time series (biphasic)
             results_ts2(index,1)=r; 
             results_ts2(index,2)=odornum; 
             results_ts2(index,3)=tmpmax2; 
             results_ts2(index,4)=0;
             results_ts2(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts2,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel2);
            figure(512);
            plot(odornum,results(r,odornum,3),'x','Color','g');
        elseif tmpmin2<minthresh %only suppressed
%             results(r,odornum,3)=0; 
%             results(r,odornum,4)=tmpmin2; 
%             results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
%             results2(index,3)=0; 
%             results2(index,4)=tmpmin2; 
%             results2(index,2)= odornum; 
%             results2(index,1)=r;
            results_trace2(index,1:length(supzts2))=exzts2;   %<<<<<<<<<<< extract time series (biphasic)
            results_ts2(index,1)=r; 
            results_ts2(index,2)=odornum; 
            results_ts2(index,3)=0; 
            results_ts2(index,4)=tmpmin2;
            results_ts2(index,5)=rgcorrcoef;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts2,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel2);
            figure(512);
            plot(odornum,results(r,odornum,4),'x','Color','b');
       else
            results_trace2(index,1:length(exzts2))=exzts2;   %<<<<<<<<<<< extract time series (biphasic)
            results_ts2(index,1)=r; 
            results_ts2(index,2)=odornum; 
            results_ts2(index,3)=0; 
            results_ts2(index,4)=0;
            results_ts2(index,5)=rgcorrcoef;
        end  
    end
%     clear exbasemean exbasestd exzts
%     clear supbasemean supbasestd supzts
%     clear tmpmax tmpmin  tmpcolor tmplabel %odornum
end

%% plot timeseries
autoArrangeFigures();

results_sig = horzcat(results_ts,results_trace);
results_sig2 = horzcat(results_ts2,results_trace2);

times = (0:size(results_sig(1,6:end),2)-1)/150; 

figure;imagesc(results_sig);colormap(nawhimar_auto);title('Channel 1 Results_Sig');
figure;imagesc(results_sig2);colormap(nawhimar_auto);title('Channel 2 Results_Sig');

%% now to find peak response latencies
xmed=[];
% for k = 1:size(results_trace,1)
%    xpeaks(k) = max(results_trace(k,iodor(1):iodor(end))*0.9);
% end

for k = 1:size(results_trace,1)
   xmed(k) = median(results_trace(k,iodor(1):iodor(end)));
end

xmed = xmed';

 peaktimes = zeros(size(results_trace,1),1);
for k = 1:size(results_trace,1)
    if xmed(k,1) > 0 
        peaktimes(k,1) = find(results_trace(k,iodor(1):iodor(end))>= xmed(k,1),1,'first');
    if xmed(k,1) < 0
        peaktimes(k,1) = find(results_trace(k,iodor(1):iodor(end))>= xmed(k,1),1,'first');
    elseif xmed(k,1) == 0 
        peaktimes(k,1) = 0;
    end
    end
end

for k = 1:size(results_trace,1)
ttpeak_ch1(k) = iodor(1)+peaktimes(k);
ttp_raw(k) = times(ttpeak_ch1(k));
ttp_odornorm(k) = ttp_raw(k)-times(iodor(1));
end

xmed2=[];
% for k = 1:size(results_trace2,1)
%    xpeaks2(k) = max(results_trace2(k,iodor(1):iodor(end))*0.9);
% end

for k = 1:size(results_trace2,1)
   xmed2(k) = median(results_trace2(k,iodor(1):iodor(end)));
end

xmed2 = xmed2';

 peaktimes2 = zeros(size(results_trace2,1),1);
for k = 1:size(results_trace2,1)
    if xmed2(k,1) > 0 
        peaktimes2(k,1) = find(results_trace2(k,iodor(1):iodor(end))>= xmed2(k,1),1,'first');
    if xmed2(k,1) < 0
        peaktimes2(k,1) = find(results_trace2(k,iodor(1):iodor(end))>= xmed2(k,1),1,'first');
    elseif xmed2(k,1) == 0 
        peaktimes2(k,1) = 0;
    end
    end
end

for k = 1:size(results_trace2,1)
ttpeak_ch2(k) = iodor(1)+peaktimes2(k);
ttp_raw2(k) = times(ttpeak_ch2(k));
ttp_odornorm2(k) = ttp_raw2(k)-times(iodor(1));
end

ttpeak_ch1 = ttpeak_ch1';
ttpeak_ch2 = ttpeak_ch2';
ttp_raw = ttp_raw';
ttp_raw2 = ttp_raw2';
ttp_odornorm = ttp_odornorm';
ttp_odornorm2 = ttp_odornorm2';

results_sig = horzcat(ttp_odornorm,results_sig);
results_sig2 = horzcat(ttp_odornorm2,results_sig2);

%% now only use those responses that show excitation in both channels

results_fpsig = results_sig(:,1:5);
results_fpsig2 = results_sig2(:,1:5);

for i = 1:length(results_fpsig)
    if results_fpsig(i,4) && results_fpsig2(i,4) > 0
        results_fpexc(i,1:5) = results_fpsig(i,1:5);
        results_fpexc2(i,1:5) = results_fpsig2(i,1:5);
    else
        results_fpexc(i,1:5) = 0;
        results_fpexc2(i,1:5) = 0;
    end
end

onlyexcsig_ch1 = results_fpexc(any(results_fpexc,2),:);
onlyexcsig_ch2 = results_fpexc2(any(results_fpexc2,2),:);

medvals_ch1 = onlyexcsig_ch1(:,1);
medvals_ch2 = onlyexcsig_ch2(:,1);

medvals_comp = horzcat(medvals_ch1,medvals_ch2);

%% now to sort both datasets by peaktimes
% [sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
% resultssig_lsort = results_sig(peakind,:);
% resultssig2_lsort = results_sig2(peakind,:);


%% save all sig responses
prompt = ('Do you want to save? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');

for i=1
    if tf == 1
       save ('results_timeseries.mat', 'results_sig', 'results_sig2', 'onlyexcsig_ch1', 'onlyexcsig_ch2')
    else
        break
    end
end
% % %% start editing here <<<<<<<<<<<<<<<<<<<<<
% % resultszmat_ex1 = zeros(length(ch1_ex.myplotdata.avgfile.roi),size(odornums,2),1502);
% % resultszmat_ex2 = zeros(length(ch2_ex.myplotdata.avgfile.roi),size(odornums,2),1502);
% % 
% % for i=1:size(odornums,2)
% %     resultszmat_ex1(:,i,:) = resultsexc1(:,odornums(i),:);
% % end
% % resultszmat_ex1(isnan(resultszmat_ex1))=0;
% % 
% % for i=1:size(odornums,2)
% %     resultszmat_ex2(:,i,:) = resultsexc2(:,odornums(i),:);
% % end
% % resultszmat_ex2(isnan(resultszmat_ex2))=0;
% % 
% % resultszmat_sup1 = zeros(length(ch1_sup.myplotdata.avgfile.roi),size(odornums,2),1502); %%now suppressive
% % resultszmat_sup2 = zeros(length(ch2_sup.myplotdata.avgfile.roi),size(odornums,2),1502);
% % 
% % for i=1:size(odornums,2)
% %     resultszmat_sup1(:,i,:) = resultssup1(:,odornums(i),:);
% % end
% % resultszmat_sup1(isnan(resultszmat_sup1))=0;
% % 
% % for i=1:size(odornums,2)
% %     resultszmat_sup2(:,i,:) = resultssup2(:,odornums(i),:);
% % end
% % resultszmat_sup2(isnan(resultszmat_sup2))=0;
% % 
% % for i=1:length(resultszmat_ex1)
% % newresultszmat_ex1(:,:,i) = vertcat(odornums,resultszmat_ex1(:,:,i));
% % end
% % 
% % for i=1:length(resultszmat_ex2)
% % newresultszmat_ex2(:,:,i) = vertcat(odornums,resultszmat_ex2(:,:,i));
% % end
% % 
% % for i=1:length(resultszmat_sup1)
% % newresultszmat_sup1(:,:,i) = vertcat(odornums,resultszmat_sup1(:,:,i));
% % end
% % 
% % for i=1:length(resultszmat_sup2)
% % newresultszmat_sup2(:,:,i) = vertcat(odornums,resultszmat_sup2(:,:,i));
% % end
% % 
% % %% convert into 2D matrices and plot waterfalls
% % for j = 2:size(newresultszmat_ex1,1)
% %     for k = 1:size(newresultszmat_ex1,2)
% %        newerresultszmat_ex1=newresultszmat_ex1(j,k,:);
% %        finalresultszmat_ex1(j*k,:) = reshape(newerresultszmat_ex1,1,1502);
% %     end
% % end
% % 
% % ex1_resultszmat = finalresultszmat_ex1(any(finalresultszmat_ex1,2),:);
% % figure;imagesc(ex1_resultszmat);colormap(nawhimar_auto);title('Ch1 Excited')
% % 
% % for j = 2:size(newresultszmat_sup1,1)
% %     for k = 1:size(newresultszmat_sup1,2)
% %        newerresultszmat_sup1 = newresultszmat_sup1(j,k,:);
% %        finalresultszmat_sup1(j*k,:) = reshape(newerresultszmat_sup1,1,1502);
% %     end
% % end
% % 
% % sup1_resultszmat = finalresultszmat_sup1(any(finalresultszmat_sup1,2),:);
% % figure;imagesc(sup1_resultszmat);colormap(nawhimar_auto);title('Ch1 Suppressed')
% % 
% % for j = 2:size(newresultszmat_ex2,1)
% %     for k = 1:size(newresultszmat_ex2,2)
% %        newerresultszmat_ex2=newresultszmat_ex2(j,k,:);
% %        finalresultszmat_ex2(j*k,:) = reshape(newerresultszmat_ex2,1,1502);
% %     end
% % end
% % 
% % ex2_resultszmat = finalresultszmat_ex2(any(finalresultszmat_ex2,2),:);
% % figure;imagesc(ex2_resultszmat);colormap(nawhimar_auto);title('Ch2 Excited')
% % 
% % for j = 2:size(newresultszmat_sup2,1)
% %     for k = 1:size(newresultszmat_sup2,2)
% %        newerresultszmat_sup2 = newresultszmat_sup2(j,k,:);
% %        finalresultszmat_sup2(j*k,:) = reshape(newerresultszmat_sup2,1,1502);
% %     end
% % end
% % 
% % sup2_resultszmat = finalresultszmat_sup2(any(finalresultszmat_sup2,2),:);
% % figure;imagesc(sup2_resultszmat);colormap(nawhimar_auto);title('Ch2 Suppressed')
% % 
% % %% normalize
% % normallexgluroisch1=[];
% % for j=1:size(ex1_resultszmat,1) %% in progress
% %     normallexgluroisch1(j,:) = normalised_diff(ex1_resultszmat(j,1:iodor(1):iodor(end):end)); %normalize to odor only
% % end
% % 
% % normallsupgluroisch1=[];
% % for j=1:size(sup1_resultszmat,1) %% in progress
% %     normallsupgluroisch1(j,:) = normalised_diff(sup1_resultszmat(j,1:iodor(1):iodor(end):end)); %normalize to odor only
% % end
% % 
% % normallexgluroisch2=[];
% % for j=1:size(ex2_resultszmat,1) %% in progress
% %     normallexgluroisch2(j,:) = normalised_diff(ex2_resultszmat(j,1:iodor(1):iodor(end):end)); %normalize to odor only
% % end
% % 
% % normallsupgluroisch2=[];
% % for j=1:size(sup2_resultszmat,1) %% in progress
% %     normallsupgluroisch2(j,:) = normalised_diff(sup2_resultszmat(j,1:iodor(1):iodor(end):end)); %normalize to odor only
% % end
% % %% sort by peak latency
% % % only excitatory ch1 and ch2
% % % ch1
% % xpeaks=[];
% % for k = 1:size(normallexgluroisch1,1)
% %    xpeaks(k) = max(normallexgluroisch1(k,iodor(1):iodor(end)));
% % end
% % for k=1:size(normallexgluroisch1,1)
% %     peaktimes(k)=find(normallexgluroisch1(k,:)==xpeaks(k));
% % end
% % 
% % [sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
% % normsortech1_l = normallexgluroisch1(peakind,:);
% % 
% % figure;imagesc(normsortech1_l);title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,C_green);
% % c = colorbar;
% % c.Label.String = 'norm. z-scored \DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([451 451 751 751],[1 1 1 1], gray);
% % hold on;
% % plot(peaktimes(peakind),1:length(peaktimes),'color','k')
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % 
% % % ch2
% % xpeaks=[];
% % for k = 1:size(normallexgluroisch2,1)
% %    xpeaks(k) = max(normallexgluroisch2(k,iodor(1):iodor(end)));
% % end
% % for k=1:size(normallexgluroisch2,1)
% %     peaktimes(k)=find(normallexgluroisch2(k,:)==xpeaks(k));
% % end
% % 
% % [sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
% % normsortech2_l = normallexgluroisch2(peakind,:);
% % 
% % figure;imagesc(normsortech2_l);title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,C_mag);
% % c = colorbar;
% % c.Label.String = 'norm. z-scored \DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([451 451 751 751],[1 1 1 1], gray);
% % hold on;
% % plot(peaktimes(peakind),1:length(peaktimes),'color','k')
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % 
% % % now suppressive
% % % ch1
% % xtroughs=[];
% % troughtimes=[];
% % for i=1
% %     if isempty(normallsupgluroisch1)
% %         normsorts = [];
% %         break
% %     else
% % %         xpeaks=[];
% %         for k = 1:size(normallsupgluroisch1,1)
% %            xtroughs(k,:) = min(normallsupgluroisch1(k,iodor(1):iodor(end)));
% %         end
% % %         peaktimes=[];
% %         for k=1:size(normallsupgluroisch1,1)
% %             troughtimes(k)=find(normallsupgluroisch1(k,:)==xtroughs(k));
% %         end
% %     end
% % 
% % 
% % [sortedtroughtimes,troughind] = sort(troughtimes,'descend');
% % normsortsch1_l = normallsupgluroisch1(troughind,:);
% % 
% % figure;imagesc(normsortsch1_l);title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([451 451 751 751],[1 1 1 1], gray);
% % hold on;
% % plot(troughtimes(troughind),1:length(troughtimes),'color','k')
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % end
% % 
% % % ch2
% % xtroughs=[];
% % troughtimes=[];
% % for i=1
% %     if isempty(normallsupgluroisch2)
% %         normsorts = [];
% %         break
% %     else
% % %         xpeaks=[];
% %         for k = 1:size(normallsupgluroisch2,1)
% %            xtroughs(k,:) = min(normallsupgluroisch2(k,iodor(1):iodor(end)));
% %         end
% % %         peaktimes=[];
% %         for k=1:size(normallsupgluroisch2,1)
% %             troughtimes(k)=find(normallsupgluroisch2(k,:)==xtroughs(k));
% %         end
% %     end
% % 
% % 
% % [sortedtroughtimes,troughind] = sort(troughtimes,'descend');
% % normsortsch2_l = normallsupgluroisch2(troughind,:);
% % 
% % figure;imagesc(normsortsch2_l);title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([451 451 751 751],[1 1 1 1], gray);
% % hold on;
% % plot(troughtimes(troughind),1:length(troughtimes),'color','k')
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % end
% % %% sort by duration
% % % only excitatory
% % % ch1
% % for k = 1:size(normallexgluroisch1,1)
% %    xfwhm(k) = fwhm(times,normallexgluroisch1(k,:));
% % end
% % 
% % [sorteddurtimes,fwhmind] = sort(xfwhm,'ascend');
% % normsortech1_d = normallexgluroisch1(fwhmind,:);
% % 
% % figure;imagesc(normsortech1_d);title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % clims = [-1 1];
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([450 450 750 750],[1 1 1 1], gray);
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % 
% % % ch2
% % for k = 1:size(normallexgluroisch2,1)
% %    xfwhm(k) = fwhm(times,normallexgluroisch2(k,:));
% % end
% % 
% % [sorteddurtimes,fwhmind] = sort(xfwhm,'ascend');
% % normsortech2_d = normallexgluroisch2(fwhmind,:);
% % 
% % figure;imagesc(normsortech2_d);title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % clims = [-1 1];
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([450 450 750 750],[1 1 1 1], gray);
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % 
% % %now suppressive
% % % ch1
% % xfwhm=[];
% % fwmind=[];
% % for i=1
% %     if isempty(normallsupglurois)
% %         normsortsch1_d = [];
% %     else 
% %         for k = 1:size(normallsupglurois,1)
% %             xfwhm(k) = fwhm(times,normallsupglurois(k,:));
% %         end
% %     end
% % [sorteddurtimes,fwhmind] = sort(xfwhm,'descend');
% % normsortsch1_d = normallsupglurois(fwhmind,:);
% % 
% % figure;imagesc(normsortsch1_d);title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([450 450 750 750],[1 1 1 1], gray);
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % %figure;bar(peaktimes(peakind),1:length(peaktimes),'color','k')
% % end
% % 
% % % ch2
% % xfwhm=[];
% % fwmind=[];
% % for i=1
% %     if isempty(normallsupglurois)
% %         normsortsch2_d = [];
% %     else 
% %         for k = 1:size(normallsupglurois,1)
% %             xfwhm(k) = fwhm(times,normallsupglurois(k,:));
% %         end
% %     end
% % [sorteddurtimes,fwhmind] = sort(xfwhm,'descend');
% % normsortsch2_d = normallsupglurois(fwhmind,:);
% % 
% % figure;imagesc(normsortsch2_d);title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Label.String = '\DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([450 450 750 750],[1 1 1 1], gray);
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % %figure;bar(peaktimes(peakind),1:length(peaktimes),'color','k')
% % end

% vals = [];
% ch1ex_fin = isfinite(results(:,:,1));
% [row,col] = find(ch1ex_fin==1);
% ch1ex_cat = horzcat(row,col);
% for r = 1:size(ch1ex_cat,1)
%     vals(r) = results(ch1ex_cat(r,1),ch1ex_cat(r,2),1);
% end
% vals = vals';
% for i=1
%     if ~isempty(vals)
% ch1ex_cat = horzcat(row,col,vals);
% ch1ex_sortcat = sortrows(ch1ex_cat);
% clearvars row col vals
%     else
%     end
% end
% 
% vals = [];
% ch1sup_fin = isfinite(results(:,:,2));
% [row,col] = find(ch1sup_fin==1);
% ch1sup_cat = horzcat(row,col);
% for r = 1:size(ch1sup_cat,1)
%     vals(r) = results(ch1sup_cat(r,1),ch1sup_cat(r,2),2);
% end
% vals = vals';
% for i=1
%     if ~isempty(vals)
% ch1sup_cat = horzcat(row,col,vals);
% ch1sup_sortcat = sortrows(ch1sup_cat);
% clearvars row col vals
%     else
%     end
% end
% 
% vals = [];
% ch2ex_fin = isfinite(results(:,:,3));
% [row,col] = find(ch2ex_fin==1);
% ch2ex_cat = horzcat(row,col);
% for r = 1:size(ch2ex_cat,1)
%     vals(r) = results(ch2ex_cat(r,1),ch2ex_cat(r,2),3);
% end
% vals = vals';
% for i=1
%     if ~isempty(vals)
% ch2ex_cat = horzcat(row,col,vals);
% ch2ex_sortcat = sortrows(ch2ex_cat);
% clearvars row col vals
%     else
%     end
% end
% 
% vals = [];
% ch2sup_fin = isfinite(results(:,:,4));
% [row,col] = find(ch2sup_fin==1);
% ch2sup_cat = horzcat(row,col);
% for r = 1:size(ch2sup_cat,1)
%     vals(r) = results(ch2sup_cat(r,1),ch2sup_cat(r,2),4);
% end
% vals = vals';
% for i=1
%     if ~isempty(vals)
% ch2sup_cat = horzcat(row,col,vals);
% ch2sup_sortcat = sortrows(ch2sup_cat);
% clearvars row col vals
%     else
%     end
% end
% 
% vals = [];
% [row,col] = find(results(:,:,5));
% chcorr_cat = horzcat(row,col);
% for r = 1:size(chcorr_cat,1)
%     vals(r) = results(chcorr_cat(r,1),chcorr_cat(r,2),5);
% end
% vals = vals';
% for i=1
%     if ~isempty(vals)
% chcorr_cat = horzcat(row,col,vals);
% chcorr_sortcat = sortrows(chcorr_cat);
% clearvars row col vals
%     else
%     end
% end
% 
% akmtestmat = zeros(38*136,7);
% rois = 38;
% for r = 1:38
%     for o = 1:136
%        akmtestmat = horzcat(r,o,results(r,o,1),results(r,o,2),results(r,o,3),results(r,o,4),results(r,o,5));
%     end
% end
% 
% clearvars -except results ch1ex_sortcat ch1sup_sortcat ch2ex_sortcat ch2sup_sortcat chcorr_sortcat
% 
% figure;scatter3(chcorr_sortcat(:,1),chcorr_sortcat(:,2),chcorr_sortcat(:,3),'MarkerFaceColor','k')
% xlabel('Glomerulus');
% ylabel('Odor');
% zlabel('Intensity (zscored dF values)');


% figure(512); hold on; %odor number vs significant min/max values
% title('Channel #2');
% exbasemeanall = [];
% supbasemeanall = [];
% %Note: If you want to plot the channel 2 timeseries in separate windows just do figure(25+whatever)
% % for r = 1:length(ch2_ex.myplotdata.avgfile.roi)
% %     for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
% %         odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
% %         exbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %         exbasemeansub = ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
% %         exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
% %         supbasemean = mean(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %         supbasemeansub = ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
% %         supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
% %     end
% %     exbasestd = std(exbasemeanall);  
% %     supbasestd = std(supbasemeanall);
% %  
% %          for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
% %         odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
% %         exbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %         %exbasestd = std(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %         exzts = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean)./exbasestd;
% %         tmpmax = prctile(exzts(iodor),95);
% %  
% %         supbasemean = mean(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %        %supbasestd = std(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
% %         supzts = (ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean)./supbasestd;
% %         tmpmin = prctile(supzts(iodor),25);
% %         tmplabel = ['Ch2, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin),', zvalmax' num2str(tmpmax)];   
%         
%         
%         if tmpmax>maxthresh && tmpmin<minthresh %plot both excitatory and suppressive
%             %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
%             results(r,odornum,3) = tmpmax; results(r,odornum,4) = tmpmin;
%             if bplotbyodor
%                 figure(odornum); hold on; %figure(256+odornum); hold on;  <-----use for separate ch2 plots
%             else
%                 figure(r); hold on; %figure(256+r); hold on;  <-----use for separate ch2 plots
%             end
%             plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'r--','DisplayName',tmplabel);
%             plot(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'m--','DisplayName',tmplabel);
%             figure(512);
%             plot(odornum,results(r,odornum,3),'x','Color','r');
%             plot(odornum,results(r,odornum,4),'x','Color','r');
%         elseif tmpmax>maxthresh %only excited
%             results(r,odornum,3)=tmpmax;
%             if bplotbyodor
%                 figure(odornum); hold on;
%                 tmpcolor = 'g'; % use 'g' or myColors(r) here
%             else
%                 figure(r); hold on;
%                 tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
%             end
%             plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel);
%             figure(512);
%             plot(odornum,results(r,odornum,3),'x','Color','g');
%         elseif tmpmin<minthresh %only suppressed
%             results(r,odornum,4) = tmpmin;
%             if bplotbyodor
%                 figure(odornum); hold on;
%                 tmpcolor = 'b'; % use 'g' or myColors(r) here
%             else
%                 figure(r); hold on;
%                 tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
%             end
%             plot(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel);
%             figure(512);
%             plot(odornum,results(r,odornum,4),'x','Color','b');
% end
%     clear exbasemean exbasestd exzts
%     clear supbasemean supbasestd supzts
%     clear tmpmax tmpmin odornum tmpcolor tmplabel
% clear o r
% %clear ch2_ex ch2_sup
% clear ibase iodor bplotbyodor
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