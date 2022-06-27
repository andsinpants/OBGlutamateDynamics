clearvars
clc
close all

[pdfile,pdpath] = uigetfile('*.mat','Choose Predrug file','MultiSelect','off');
if ~pdfile; return; end
pdrug = importdata(fullfile(pdpath,pdfile)); clear pdpath pdfile;
[cgpfile,cgppath] = uigetfile('*.mat','Choose CGP35348 file','MultiSelect','off');
if ~cgpfile; return; end
cgpdrug = importdata(fullfile(cgppath,cgpfile)); clear cgppath cgpfile;
[iglurfile,iglurpath] = uigetfile('*.mat','Choose APV+NBQX file','MultiSelect','off');
if ~iglurfile; return; end
iglurdrug = importdata(fullfile(iglurpath,iglurfile)); clear iglurpath iglurfile;

bplotbyodor = 1;

times = pdrug.file.roi(1).odor.avgtrial.time;

%% to use for prompting 
prompt = 'What is the odor start time (1st sniff)?';
p = input(prompt);
g = find(times >= p,1,'first');
x = g-150;
y = g;

% prompt = 'What is the dF start time (1st sniff)?';
a = g+15;
b = a+150;

prompt = 'What is the odor end time (odor duration)?';
o = input(prompt);
k = find(times >= o,1,'first');

c = x;
d = y;

w = k-525;
z = k;

ibase = x:y; %preodor frames (make an input variable, hard coded for pcd235)
iodor = a:b; %postodor frames 
ibase2 = c:d; %preodor frames 2nd sniff
iodor2 = w:z; %postodor frames 2nd sniff
ibase2alt = 50:108; %alt baseline (preodor)

% % prompt = 'What is the dF start (1st sniff)? ';
% % a = input(prompt);
% % b = a+5;
% uf = 10;
% %the following works well for peak amp and time calculations
% prompt = 'What is the frame of the first sniff onset? ';
% a = input(prompt);
% a = a*uf;   %first sniff onset
% b = a+(25*uf);   %1 sec after first sniff onset = second sniff onset for Hz sniff
% y = a;
% x = a-(5*uf); 
% d = b;   %use for odor ontimefor second sniff
% c= b-(5*uf);  % 5NP frames before second sniff onset
% w = d;
% z = d+(25*uf);   %1 sec after second sniff onset

%% now calculate

maxthresh = 0; 
minthresh = 0; %default 7 (multiplied by 4 as it seems right), use of 6 or 7 could work

figure(256); hold on; %odor number vs significant min/max values
title('Predrug');
exbasemeanall = [];
supbasemeanall = [];
index=0;

for r=1:length(pdrug.file.roi)
    pdrug.data(r,:) = pdrug.file.roi(r).odor.avgtrial.series';
end

for r=1:length(cgpdrug.file.roi)
cgpdrug.data(r,:) = cgpdrug.file.roi(r).odor.avgtrial.series';
end

for r=1:length(iglurdrug.file.roi)
iglurdrug.data(r,:) = iglurdrug.file.roi(r).odor.avgtrial.series';
end

uf = 10;
%times = (0:size(pdrug.data(1,:),2)-1)/150; 
%% now same sorting and plotting of predrug condition data (1st sniff)
for r = 1:size(pdrug.data,1)
     for o = 1:length(1)
        exbasemean = mean(pdrug.data(r,ibase));
        exbasemeansub = pdrug.data(r,ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);  
        
        supbasemean = mean(pdrug.data(r,ibase));
        supbasemeansub = pdrug.data(r,ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);

     end
     
     exbasestd = std(exbasemeanall');
     supbasestd = std(supbasemeanall');
 
     %%%START HERE%%%
    
    for o = 1:length(1)
        odornum = 1;
        exbasemean = mean(pdrug.data(r,ibase(1):ibase(end)));
        exzts = (pdrug.data(r,:)-exbasemean)./exbasestd(1,r);
%       tmpmax = prctile(exzts(130:135),95);
        tmpmean = mean(exzts(iodor(1):iodor(end)));
%       tmpmax = tmpmean;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(pdrug.data(r,ibase(1):ibase(end)));
        %supbasestd = std(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (pdrug.data(r,:)-supbasemean)./supbasestd(1,r);
        tmpmin = prctile(supzts(iodor),15); %default is 25
        
        
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum)];   
        index=index+1;

        if  tmpmax>maxthresh %only excited
             results(r,odornum,1)=tmpmax;
             results_trace(index,1:length(exzts))=pdrug.data(r,:);  %<<<<<<<<<<< extract time series (exc)
             results_ts(index,2)=odornum;
             results_ts(index,1)=r;
             results_ts(index,3)=tmpmax; 
             results_ts(index,4)=0;
            if bplotbyodor 
                figure(odornum); hold on;
                tmpcolor = 'r'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(1,results(r,odornum,1),'x','Color','r');
        elseif tmpmin<minthresh %only suppressed
           results(r,odornum,2)=tmpmin; 
           results_trace(index,1:length(supzts))=pdrug.data(r,:); %<<<<<<<<<<<<< extract time series (sup)
           results_ts(index,2)= odornum;
           results_ts(index,1)=r;
           results_ts(index,3)=0; 
           results_ts(index,4)=tmpmin;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,2),'x','Color','b');
        else
           results(r,odornum,2)=0; 
           results_trace(index,1:length(exzts))=0; 
           results_ts(index,2)= odornum;
           results_ts(index,1)=r;
           results_ts(index,3)=0; 
           results_ts(index,4)=0;
        end
    end
end

%% now same sorting and plotting of CGP condition data (1st sniff)
clearvars exzts exbasemean exbasemeansub exbasemeanll exbasestd tmpmax supzts supbasemean supbasemeansub supbasemeanll supbasestd tmpmin 
index=0;

for r = 1:size(cgpdrug.data,1)
     for o = 1:length(1)
        exbasemean = mean(cgpdrug.data(r,ibase));
        exbasemeansub = cgpdrug.data(r,ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);  
        
        supbasemean = mean(cgpdrug.data(r,ibase));
        supbasemeansub = cgpdrug.data(r,ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);

     end
     
     exbasestd = std(exbasemeanall');
     supbasestd = std(supbasemeanall');
 
     %%%START HERE%%%
    
    for o = 1:length(1)
        odornum = 1;
        exbasemean = mean(cgpdrug.data(r,ibase(1):ibase(end)));
        exzts = (cgpdrug.data(r,:)-exbasemean)./exbasestd(1,r);
%       tmpmax = prctile(exzts(130:135),95);
        tmpmean = mean(exzts(iodor(1):iodor(end)));
%       tmpmax = tmpmean;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(cgpdrug.data(r,ibase(1):ibase(end)));
        %supbasestd = std(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (cgpdrug.data(r,:)-supbasemean)./supbasestd(1,r);
        tmpmin = prctile(supzts(iodor),15); %default is 25
        
        
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum)];   
        index=index+1; 
        
        k(r) = results_ts(r,3)~= 0;
        p(r) = results_ts(r,4)~= 0;

        if k(r) == 1 %only excited
           results_cgp(r,odornum,1)=tmpmax;
           results_trace_cgp(index,1:length(exzts))=cgpdrug.data(r,:);  %<<<<<<<<<<< extract time series (exc)
           results_ts_cgp(index,2)=odornum;
           results_ts_cgp(index,1)=r;
           results_ts_cgp(index,3)=tmpmax; 
           results_ts_cgp(index,4)=0;
            if bplotbyodor 
                figure(odornum); hold on;
                tmpcolor = 'm'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(2,results_cgp(r,odornum,1),'x','Color','r');
        elseif p(r) == 1 %only suppressed
           results_cgp(r,odornum,2)=tmpmin; 
           results_trace_cgp(index,1:length(supzts))=cgpdrug.data(r,:); %<<<<<<<<<<<<< extract time series (sup)
           results_ts_cgp(index,2)=odornum;
           results_ts_cgp(index,1)=r;
           results_ts_cgp(index,3)=0; 
           results_ts_cgp(index,4)=tmpmin;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(2,results_cgp(r,odornum,2),'x','Color','b');
        else
           results_cgp(r,odornum,2)=0; 
           results_trace_cgp(index,1:length(exzts))=0; 
           results_ts_cgp(index,2)=odornum;
           results_ts_cgp(index,1)=r;
           results_ts_cgp(index,3)=0; 
           results_ts_cgp(index,4)=0;
        end
    end
end

%% now same sorting and plotting of APV+NBQX condition data (1st sniff)
clearvars exzts exbasemean exbasemeansub exbasemeanll exbasestd tmpmax supzts supbasemean supbasemeansub supbasemeanll supbasestd tmpmin 
index=0;

for r = 1:size(iglurdrug.data,1)
     for o = 1:length(1)
        exbasemean = mean(iglurdrug.data(r,ibase));
        exbasemeansub = cgpdrug.data(r,ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);  
        
        supbasemean = mean(iglurdrug.data(r,ibase));
        supbasemeansub = iglurdrug.data(r,ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);

     end
     
     exbasestd = std(exbasemeanall');
     supbasestd = std(supbasemeanall');
 
     %%%START HERE%%%
    
    for o = 1:length(1)
        odornum = 1;
        exbasemean = mean(iglurdrug.data(r,ibase(1):ibase(end)));
        exzts = (iglurdrug.data(r,:)-exbasemean)./exbasestd(1,r);
%       tmpmax = prctile(exzts(130:135),95);
        tmpmean = mean(exzts(iodor(1):iodor(end)));
%       tmpmax = tmpmean;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(iglurdrug.data(r,ibase(1):ibase(end)));
        %supbasestd = std(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (iglurdrug.data(r,:)-supbasemean)./supbasestd(1,r);
        tmpmin = prctile(supzts(iodor),15); %default is 25
        
        
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum)];   
        index=index+1;

        k(r) = results_ts(r,3)~= 0;
        p(r) = results_ts(r,4)~= 0;
        
        if  k(r) == 1 %only excited
           results_iglur(r,odornum,1)=tmpmax;
           results_trace_iglur(index,1:length(exzts))=iglurdrug.data(r,:);  %<<<<<<<<<<< extract time series (exc)
           results_ts_iglur(index,2)=odornum;
           results_ts_iglur(index,1)=r;
           results_ts_iglur(index,3)=tmpmax; 
           results_ts_iglur(index,4)=0;
            if bplotbyodor 
                figure(odornum); hold on;
                tmpcolor = 'y'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(3,results_iglur(r,odornum,1),'x','Color','r');
        elseif p(r) == 1 %only suppressed
           results_iglur(r,odornum,2)=tmpmin; 
           results_trace_iglur(index,1:length(supzts))=iglurdrug.data(r,:); %<<<<<<<<<<<<< extract time series (sup)
           results_ts_iglur(index,2)=odornum;
           results_ts_iglur(index,1)=r;
           results_ts_iglur(index,3)=0; 
           results_ts_iglur(index,4)=tmpmin;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'k'; % use 'g' or myColors(r) here
            else 
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(times,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(3,results_iglur(r,odornum,2),'x','Color','b');
        else
           results_iglur(r,odornum,2)=0; 
           results_trace_iglur(index,1:length(exzts))=0; 
           results_ts_iglur(index,2)=odornum;
           results_ts_iglur(index,1)=r;
           results_ts_iglur(index,3)=0; 
           results_ts_iglur(index,4)=0;
        end
    end
end


figure(1);line([times(a),times(a)],[0,10],'Color','red','LineWidth',2);
figure(1);line([times(b),times(b)],[0,10],'Color','black','LineWidth',2);
figure(1);line([times(w),times(w)],[0,10],'Color','blue','LineWidth',2);
figure(1);line([times(z),times(z)],[0,10],'Color','green','LineWidth',2);
pause

%% now calculate sniff 2-sniff 1 across all ROIs in each condition
figure;
for i = 1:size(pdrug.data,1)
    firstsniff_pdrug(i) = mean(results_trace(i,iodor))-mean(results_trace(i,ibase));
    secondsniff_pdrug(i) = mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2));
    secandfirsniff_pdiff(i) = secondsniff_pdrug(i)-firstsniff_pdrug(i);
    concatsniff_pdrug = vertcat(firstsniff_pdrug,secondsniff_pdrug);
    j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0; 
    if j(i) == 1
    plot([1 2],[firstsniff_pdrug(i) secondsniff_pdrug(i)],'--o');hold on;title('Sniff 1 and All Odor Mean (Predrug)');xlabel('Sniff 1 || All Odor Mean');ylabel('{\Delta}F/F');
    end
end

figure;
for i = 1:size(cgpdrug.data,1)
    firstsniff_cgpdrug(i) = mean(results_trace_cgp(i,iodor))-mean(results_trace_cgp(i,ibase));
    secondsniff_cgpdrug(i) = mean(results_trace_cgp(i,iodor2))-mean(results_trace_cgp(i,ibase2));
    concatsniff_cgpdrug = vertcat(firstsniff_cgpdrug,secondsniff_cgpdrug);
    secandfirsniff_cgpdiff(i) = secondsniff_cgpdrug(i)-firstsniff_cgpdrug(i);
    j(i) = concatsniff_cgpdrug(1,i)~= 0 && concatsniff_cgpdrug(2,i)~= 0; 
    if j(i) == 1
    plot([1 2],[firstsniff_cgpdrug(i) secondsniff_cgpdrug(i)],'--o');hold on;title('Sniff 1 and All Odor Mean (CGP)');xlabel('Sniff 1 || All Odor Mean');ylabel('{\Delta}F/F');
    end
end

figure;
for i = 1:size(iglurdrug.data,1)
    firstsniff_iglurdrug(i) = mean(results_trace_iglur(i,iodor))-mean(results_trace_iglur(i,ibase));
    secondsniff_iglurdrug(i) = mean(results_trace_iglur(i,iodor2))-mean(results_trace_iglur(i,ibase2));
    secandfirsniff_iglurdiff(i) = secondsniff_iglurdrug(i)-firstsniff_iglurdrug(i);
    concatsniff_iglurdrug = vertcat(firstsniff_iglurdrug,secondsniff_iglurdrug);
    j(i) = concatsniff_iglurdrug(1,i)~= 0 && concatsniff_iglurdrug(2,i)~= 0; 
    if j(i) == 1
    plot([1 2],[firstsniff_iglurdrug(i) secondsniff_iglurdrug(i)],'--o');hold on;title('Sniff 1 and All Odor Mean (APV+NBQX)');xlabel('Sniff 1 || All Odor Mean');ylabel('{\Delta}F/F');
    end
end

%% now take the same ROI and plot the Sniff 2-Sniff 1 diff across each condition
figure;
for i = 1:size(pdrug.data,1)
    j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0;
    if j(i) == 1
        plot([1 2 3],[secandfirsniff_pdiff(i) secandfirsniff_cgpdiff(i) secandfirsniff_iglurdiff(i)],'--x');hold on;title('{\Delta}F All Odor Mean - {\Delta}F Sniff 1 (Across Conditions)');xlabel('Condition');ylabel('Change in {\Delta}F/F');
    end
end


%% max(iodor)-min(ibase)
figure;
close all
for i = 1:size(pdrug.data,1)
    
    normmaxmin_pdrug_fs(i) = (max(results_trace(i,iodor))-min(results_trace(i,ibase)))./(max(results_trace(i,iodor))-min(results_trace(i,ibase)));
    normmaxmin_cgpdrug_fs(i) = (max(results_trace_cgp(i,iodor))-min(results_trace_cgp(i,ibase)))./(max(results_trace(i,iodor))-min(results_trace(i,ibase)));
    normmaxmin_iglurdrug_fs(i) = (max(results_trace_iglur(i,iodor))-min(results_trace_iglur(i,ibase)))./(max(results_trace(i,iodor))-min(results_trace(i,ibase)));
    
    if results_ts(i,3) ~= 0
        
    
    % temptrace_pdrug=[]; R = []; peakamp1 = []; peakindex = []; %clearvars in loop    
    [peakampl, peakindex]=max(results_trace(i,iodor));
        peaklat_pdrug(i)=peakindex*(64.5/uf);
    temptrace_pdrug=results_trace(i,iodor);
    mean_pre=mean(temptrace_pdrug(1:25));
    %     figure(i);plot(temptrace_pdrug);hold on;
    temptrace_pdrug(peakindex:end)=peakampl;
    figure(i);plot(temptrace_pdrug);hold on;
    [R,LT_pdrug,UT_pdrug,LL_pdrug,UL_pdrug]=risetime(temptrace_pdrug,(15.5*uf),'StateLevels',[mean_pre peakampl]);
     figure(i);line([1,size(temptrace_pdrug,2)],[UL_pdrug, UL_pdrug],'Color','red','LineWidth',2);line([1,size(temptrace_pdrug,2)],[mean_pre, mean_pre],'Color','magenta','LineWidth',2)
%     figure(i);line([1,size(temptrace_pdrug,2)],[LL_pdrug, LL_pdrug],'Color','red','LineWidth',2);line([1,size(temptrace_pdrug,2)],[mean_pre, mean_pre],'Color','magenta','LineWidth',2)
    if isempty(LT_pdrug) || isempty(UT_pdrug)
        LT_pdrug = 0
        UT_pdrug = 0
    end
    onsets_pdrug(i)=LT_pdrug(1)*1000;  %convert seconds to milliseconds
    peaklats90_pdrug(i)=UT_pdrug(1)*1000;
    
    
    %temptrace_cgpdrug=[]; R = []; peakamp1 = []; %peakindex = []; %clearvars in loop
    [peakampl, peakindex]=max(results_trace_cgp(i,iodor));
    peaklat_cgpdrug(i)=peakindex*(64.5/uf);
    temptrace_cgpdrug=results_trace_cgp(i,iodor);
    mean_pre=mean(temptrace_cgpdrug(1:25));
    %     figure(i);plot(temptrace_cgpdrug);hold on;
    temptrace_cgpdrug(peakindex:end)=peakampl;
    figure(i);plot(temptrace_cgpdrug);
    [R,LT_cgpdrug,UT_cgpdrug,LL_cgpdrug,UL_cgpdrug]=risetime(temptrace_cgpdrug,(15.5*uf),'StateLevels',[mean_pre peakampl]);
     figure(i);line([1,size(temptrace_cgpdrug,2)],[UL_cgpdrug, UL_cgpdrug],'Color','blue','LineWidth',2);line([1,size(temptrace_cgpdrug,2)],[mean_pre, mean_pre],'Color','green','LineWidth',2)
%     figure(i);line([1,size(temptrace_cgpdrug,2)],[LL_cgpdrug, LL_cgpdrug],'Color','blue','LineWidth',2);line([1,size(temptrace_cgpdrug,2)],[mean_pre, mean_pre],'Color','green','LineWidth',2)
    if isempty(LT_cgpdrug) || isempty(UT_cgpdrug)
        LT_cgpdrug = 0
        UT_cgpdrug = 0
    end
    onsets_cgp(i)=LT_cgpdrug(1)*1000;  %convert seconds to milliseconds
    peaklats90_cgpdrug(i)=UT_cgpdrug(1)*1000;
    
    %temptrace_iglur=[]; R = []; peakamp1 = []; %peakindex = []; %clearvars in loop
    [peakampl, peakindex]=max(results_trace_iglur(i,iodor));
    peaklat_iglur(i)=peakindex*(64.5/uf);
    temptrace_iglur=results_trace_iglur(i,iodor);
    mean_pre=mean(temptrace_iglur(1:25));
    %     figure(i);plot(temptrace_iglur);hold on;
    temptrace_iglur(peakindex:end)=peakampl;
    figure(i);plot(temptrace_iglur);
    [R,LT_iglur,UT_iglur,LL_iglur,UL_iglur]=risetime(temptrace_iglur,(15.5*uf),'StateLevels',[mean_pre peakampl]);
     figure(i);line([1,size(temptrace_iglur,2)],[UL_iglur, UL_iglur],'Color','cyan','LineWidth',2);line([1,size(temptrace_iglur,2)],[mean_pre, mean_pre],'Color','yellow','LineWidth',2);hold off
%     figure(i);line([1,size(temptrace_iglur,2)],[LL_iglur, LL_iglur],'Color','cyan','LineWidth',2);line([1,size(temptrace_iglur,2)],[mean_pre, mean_pre],'Color','yellow','LineWidth',2);hold off
    if isempty(LT_iglur) || isempty(UT_iglur)
        LT_iglur = 0
        UT_iglur = 0
    end
    onsets_iglur(i)=LT_iglur(1)*1000;  %convert seconds to milliseconds
    peaklats90_iglur(i)=UT_iglur(1)*1000;
      
    else
        R = 0;
        LT_iglur = 0;
        UT_iglur = 0;
    end
    j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0;
%     if j(i) == 1
%         plot([1 2 3],[normmaxmin_pdrug_fs(i) normmaxmin_cgpdrug_fs(i) normmaxmin_iglurdrug_fs(i)],'--x');hold on;title('Sniff 1 dF (max(iodor)-min(ibase)) values normalized (Across Conditions)');xlabel('Condition (Predrug || CGP35348 || APV+NBQX');ylabel('{\Delta}F Normalized to Predrug'); 
%     end
end

%% now to concantenate to compare and check
concat_onsetpeaklat_pdrug = vertcat(onsets_pdrug,peaklats90_pdrug,peaklat_pdrug);
concat_onsetpeaklat_cgpdrug = vertcat(onsets_cgp,peaklats90_cgpdrug,peaklat_cgpdrug);
concat_onsetpeaklat_iglur = vertcat(onsets_iglur,peaklats90_iglur,peaklat_iglur);

%% (mean(iodor2)-mean(ibase2)

figure;
for i = 1:size(pdrug.data,1)
    normmean_pdrug_ss(i) = (mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2)))./(mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2)));
    normmean_cgpdrug_ss(i) = (mean(results_trace_cgp(i,iodor2))-mean(results_trace_cgp(i,ibase2)))./(mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2)));
    normmean_iglurdrug_ss(i) = (mean(results_trace_iglur(i,iodor2))-mean(results_trace_iglur(i,ibase2)))./(mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2)));
j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0;
    if j(i) == 1
        plot([1 2 3],[normmean_pdrug_ss(i) normmean_cgpdrug_ss(i) normmean_iglurdrug_ss(i)],'--x');hold on;title('Sniff 1 dF (mean(iodor)-mean(ibase)) values normalized (Across Conditions)');xlabel('Condition (Predrug || CGP35348 || APV+NBQX');ylabel('{\Delta}F Normalized to Predrug'); 
    end
end

%% bar plot
for i = 1:length(normmaxmin_pdrug_fs)
figure;bar(1,normmaxmin_pdrug_fs(i));hold on;bar(2,normmaxmin_cgpdrug_fs(i));bar(3,normmaxmin_iglurdrug_fs(i));
end

% close all

%% mean(iodor2)-mean(ibase2) NOT NORMALIZED
figure;
for i = 1:size(pdrug.data,1)
    mean_pdrug_ss(i) = (mean(results_trace(i,iodor2))-mean(results_trace(i,ibase2)));
    mean_cgpdrug_ss(i) = (mean(results_trace_cgp(i,iodor2))-mean(results_trace_cgp(i,ibase2)));
    mean_iglurdrug_ss(i) = (mean(results_trace_iglur(i,iodor2))-mean(results_trace_iglur(i,ibase2)));
j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0;
    if j(i) == 1
        plot([1 2 3],[mean_pdrug_ss(i) mean_cgpdrug_ss(i) mean_iglurdrug_ss(i)],'--x');hold on;title('mean dF across odor duration(mean(iodor2)-mean(ibase2)) values non-normalized (Across Conditions)');xlabel('Condition (Predrug || CGP35348 || APV+NBQX');ylabel('{\Delta}F Normalized to Predrug'); 
    end
end

%% max(iodor)-min(ibase) NOT NORMALIZED
figure;
for i = 1:size(pdrug.data,1)
    maxmin_pdrug_fs(i) = (max(results_trace(i,iodor))-min(results_trace(i,ibase)));
    maxmin_cgpdrug_fs(i) = (max(results_trace_cgp(i,iodor))-min(results_trace_cgp(i,ibase)));
    maxmin_iglurdrug_fs(i) = (max(results_trace_iglur(i,iodor))-min(results_trace_iglur(i,ibase)));
j(i) = concatsniff_pdrug(1,i)~= 0 && concatsniff_pdrug(2,i)~= 0;
    if j(i) == 1
        plot([1 2 3],[maxmin_pdrug_fs(i) maxmin_cgpdrug_fs(i) maxmin_iglurdrug_fs(i)],'--x');hold on;title('Sniff 1 dF (max(iodor)-min(ibase)) values non-normalized (Across Conditions)');xlabel('Condition (Predrug || CGP35348 || APV+NBQX');ylabel('{\Delta}F Normalized to Predrug'); 
    end
end

maxmin_pdrug_fs = maxmin_pdrug_fs';
mean_pdrug_ss = mean_pdrug_ss';
mean_cgpdrug_ss = mean_cgpdrug_ss';
maxmin_cgpdrug_fs = maxmin_cgpdrug_fs';
mean_iglurdrug_ss = mean_iglurdrug_ss';
maxmin_iglurdrug_fs = maxmin_iglurdrug_fs';



%% now do peak lats based on 90% max
figure;scatter(peaklats90_pdrug,peaklats90_cgpdrug,'o','g','filled'); 
title('Normalized \DeltaF/F Unity Plot (Predrug vs CGP35348)');
xlabel('Predrug');
ylabel('CGP35348');
ylimit = ylim;
xlimit = xlim;
newmin = min(ylimit(1),xlimit(1)); newmax = max(ylimit(2),xlimit(2));
xlim([newmin newmax]); ylim([newmin newmax]);
lgd = legend('90 perc of Peak Latencies','Unity line');
refline(1,0);

figure;scatter(peaklats90_cgpdrug,peaklats90_iglur,'o','g','filled'); 
title('Normalized \DeltaF/F Unity Plot (CGP35348 vs APV+NBQXAPV/NBQX)');
xlabel('CGP35348');
ylabel('APV+NBQX');
ylimit = ylim;
xlimit = xlim;
newmin = min(ylimit(1),xlimit(1)); newmax = max(ylimit(2),xlimit(2));
xlim([newmin newmax]); ylim([newmin newmax]);
lgd = legend('90 perc of Peak Latencies','Unity line');
refline(1,0);

%% now plot onset latencies

figure;scatter(onsets_pdrug,onsets_cgp,'o','b','filled'); 
title('Normalized \DeltaF/F Unity Plot (Predrug vs CGP35348)');
xlabel('Predrug');
ylabel('CGP35348');
ylimit = ylim;
xlimit = xlim;
newmin = min(ylimit(1),xlimit(1)); newmax = max(ylimit(2),xlimit(2));
xlim([newmin newmax]); ylim([newmin newmax]);
lgd = legend('Onset Latencies','Unity line');
refline(1,0);

figure;scatter(onsets_cgp,onsets_iglur,'o','b','filled'); 
title('Normalized \DeltaF/F Unity Plot (CGP35348 vs APV+NBQX)');
xlabel('CGP35348');
ylabel('APV/NBQX');
ylimit = ylim;
xlimit = xlim;
newmin = min(ylimit(1),xlimit(1)); newmax = max(ylimit(2),xlimit(2));
xlim([newmin newmax]); ylim([newmin newmax]);
lgd = legend('Onset Latencies','Unity line');
refline(1,0);

clearvars -except maxmin_pdrug_fs mean_pdrug_ss mean_cgpdrug_ss mean_iglurdrug_ss maxmin_cgpdrug_fs maxmin_iglurdrug_fs onsets_pdrug onsets_cgp onsets_iglur peaklats90_pdrug peaklats90_cgpdrug peaklats90_iglur results_trace results_trace_cgp results_trace_iglur results_ts results_ts_cgp results_ts_iglur temptrace_pdrug temptrace_cgp temptrace_iglur


%% now transpose latency variables
onsets_pdrug = onsets_pdrug';
onsets_cgp = onsets_cgp';
onsets_iglur = onsets_iglur';
peaklats90_pdrug = peaklats90_pdrug';
peaklats90_cgpdrug = peaklats90_cgpdrug';
peaklats90_iglur = peaklats90_iglur';

for i=1:length(onsets_pdrug)
    normpratdrug_fs(i) = maxmin_pdrug_fs(i)/max(maxmin_pdrug_fs);
    normcgpratdrug_fs(i) = maxmin_cgpdrug_fs(i)/max(maxmin_pdrug_fs);
    normiglurratdrug_fs(i) = maxmin_iglurdrug_fs(i)/max(maxmin_pdrug_fs);
    normpratdrug_ss(i) = mean_pdrug_ss(i)/max(mean_pdrug_ss);
    normcgpratdrug_ss(i) = mean_cgpdrug_ss(i)/max(mean_pdrug_ss);
    normiglurratdrug_ss(i) = mean_iglurdrug_ss(i)/max(mean_pdrug_ss);
end

normpratdrug_fs = normpratdrug_fs';
normcgpratdrug_fs = normcgpratdrug_fs';
normiglurratdrug_fs = normiglurratdrug_fs';
normpratdrug_ss = normpratdrug_ss';
normcgpratdrug_ss = normcgpratdrug_ss';
normiglurratdrug_ss = normiglurratdrug_ss';
    
for i=1:length(onsets_pdrug)
fsrat_pdrugcgp(i) = normcgpratdrug_fs(i)/normpratdrug_fs(i);
fsrat_cgpiglur(i) = normiglurratdrug_fs(i)/normcgpratdrug_fs(i);

ssrat_pdrugcgp(i) = normcgpratdrug_ss(i)/normpratdrug_ss(i);
ssrat_cgpiglur(i) = normiglurratdrug_ss(i)/normcgpratdrug_ss(i);

onsetdiff_pdrugcgp(i) = onsets_cgp(i)-onsets_pdrug(i);
onsetdiff_cgpiglur(i) = onsets_iglur(i)-onsets_cgp(i);

peakdiff_pdrugcgp(i) = peaklats90_cgpdrug(i)-peaklats90_pdrug(i);
peakdiff_cgpiglur(i) = peaklats90_iglur(i)-peaklats90_cgpdrug(i);
end

medom_fs_pdrugcgp = median(fsrat_pdrugcgp(:));
meanom_fs_pdrugcgp = mean(fsrat_pdrugcgp(:));
stdom_fs_pdrugcgp = std(fsrat_pdrugcgp(:));

medom_ss_pdrugcgp = median(ssrat_pdrugcgp(:));
meanom_ss_pdrugcgp = mean(ssrat_pdrugcgp(:));
stdom_ss_pdrugcgp = std(ssrat_pdrugcgp(:));

medom_fs_cgpiiglur = median(fsrat_cgpiglur(:));
meanom_fs_cgpiglur = mean(fsrat_cgpiglur(:));
stdom_fs_cgpiglur = std(fsrat_cgpiglur(:));

medom_ss_cgpiglur = median(ssrat_cgpiglur(:));
meanom_ss_cgpiglur = mean(ssrat_cgpiglur(:));
stdom_ss_cgpiglur = std(ssrat_cgpiglur(:));

medom_onset_pdrugcgp = median(onsetdiff_pdrugcgp(:));
meanom_onset_pdrugcgp = mean(onsetdiff_pdrugcgp(:));
stdom_onset_pdrugcgp = std(onsetdiff_pdrugcgp(:));

medom_onset_cgpiglur = median(onsetdiff_cgpiglur(:));
meanom_onset_cgpiglur = mean(onsetdiff_cgpiglur(:));
stdom_onset_cgpiglur = std(onsetdiff_cgpiglur(:));

medom_peaklat_pdrugcgp = median(peakdiff_pdrugcgp(:));
meanom_peaklat_pdrugcgp = mean(peakdiff_pdrugcgp(:));
stdom_peaklat_pdrugcgp = std(peakdiff_pdrugcgp(:));

medom_peaklat_cgpiglur = median(peakdiff_cgpiglur(:));
meanom_peaklat_cgpiglur = mean(peakdiff_cgpiglur(:));
stdom_peaklat_cgpiglur = std(peakdiff_cgpiglur(:));


autoArrangeFigures();
%pause;
%close all;

figure(23);
for i=1:length(onsets_pdrug)
    plot([1.5 2.5],[fsrat_pdrugcgp(i) fsrat_cgpiglur(i)],'--x');
    title('1^s^t sniff amplitude ratio (Postdrug/Previous Condition)');
    xlabel('CGP35348/Predrug -------- APV+NBQX/CGP35348');
    ylabel('Post/Pre Value');
    xlim([1 3]);
    hold on;
end

figure(24);
for i=1:length(onsets_pdrug)
    plot([1.5 2.5],[ssrat_pdrugcgp(i) ssrat_cgpiglur(i)],'--x');
    title('All Odor Mean ratio (Postdrug/Previous Condition)');
    xlabel('CGP35348/Predrug -------- APV+NBQX/CGP35348');
    ylabel('Post/Pre Value');
    xlim([1 3]);
    hold on;
end

figure(25);
for i=1:length(onsets_pdrug)
    plot([1.5 2.5],[onsetdiff_pdrugcgp(i) onsetdiff_cgpiglur(i)],'--x');
    title('Onset Latency Differences (Postdrug-Previous Condition)');
    xlabel('CGP35348-Predrug -------- APV+NBQX-CGP35348');
    ylabel('Post-Pre Value');
    xlim([1 3]);
    hold on;
end

figure(26);
for i=1:length(onsets_pdrug)
    plot([1.5 2.5],[peakdiff_pdrugcgp(i) peakdiff_cgpiglur(i)],'--x');
    title('Peak Latency Differences (Postdrug-Previous Condition)');
    xlabel('CGP35348-Predrug -------- APV+NBQX-CGP35348');
    ylabel('Post-Pre Value');
    xlim([1 3]);
    hold on;
end

clearvars i

autoArrangeFigures();