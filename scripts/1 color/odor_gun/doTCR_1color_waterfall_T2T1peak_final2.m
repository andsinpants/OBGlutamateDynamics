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
tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmax' num2str(tmpmax)];   
        index=index+1;

        if tmpmax>maxthresh %only excited
             results(r,odornum,1)=prepdiff;
             results_trace(index,1:length(prezts))=prezts;  %<<<<<<<<<<< extract time series (exc)
             results_ts(index,2)= odornum;
             results_ts(index,1)=r;
             results_ts(index,3)=prepdiff; results_ts(index,4)=0;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','k');
        elseif tmpmax <= maxthresh %
           results(r,odornum,2)=postpdiff; 
           results_trace(index,1:length(postzts))=prezts; %<<<<<<<<<<<<< extract time series (sup)
           results_ts(index,2)= odornum;
           results_ts(index,1)=r;
           results_ts(index,3)=postpdiff; results_ts(index,4)=0;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,postzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,2),'x','Color','g');
        end
    end
end
results_all=horzcat(results_ts, results_trace);
%clearvars -except results_all
autoArrangeFigures();

%% now remove arrays with all 0s
results_sig = results_all(any(results_all,2),:);

times = (0:size(results_sig(1,5:end),2)-1)/150; 
%% quick plot of figure
figure;imagesc(results_sig);colormap(nawhimar_auto);