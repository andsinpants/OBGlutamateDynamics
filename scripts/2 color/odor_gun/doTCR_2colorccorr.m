close all; 
clear;

[exfile,expath] = uigetfile('*.mat','Choose Channel #1 Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ch1_ex = load(fullfile(expath,exfile)); clear expath exfile;
[supfile,suppath] = uigetfile('*.mat','Choose Channel #1 Suppressive File(s) (2Hz filter?)','MultiSelect','off');
if ~supfile; return; end
ch1_sup = load(fullfile(suppath,supfile)); clear suppath supfile;

%channel 2
[exfile,expath] = uigetfile('*.mat','Choose Channel #2 Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; clear ibase iodor bplotbyodor; return; end
ch2_ex = load(fullfile(expath,exfile)); clear expath exfile;
[supfile,suppath] = uigetfile('*.mat','Choose Channel #2 Suppressive File(s) (2Hz filter?)','MultiSelect','off');
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
results = nan(length(ch1_ex.myplotdata.avgfile.roi),odornums(end),2);
results(:,:,3:4) = nan(size(results,1),size(results,2),2);   %inceased to add one more column for 5 values

figure(256); hold on; %odor number vs significant min/max values
title('Channel #1');
exbasemeanall = [];
supbasemeanall = [];
exbasemeanall2 = [];
supbasemeanall2 = [];
for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
     for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
        
        supbasemean = mean(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
     
     %make usre these variables have unique names for ch2.
        %odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
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
        tmpmin = prctile(supzts(iodor),25);
        
        %this is for the channel 2 calculations
        exbasemean2 = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %exbasestd = std(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exzts2 = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2)./exbasestd2;
        tmpmax2 = prctile(exzts2(iodor),95);
        supbasemean2 = mean(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
       %supbasestd = std(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts2 = (ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean2)./supbasestd2;
        tmpmin2 = prctile(supzts2(iodor),25);
        
        rgcorr=corrcoef(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean, ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2);
%         rgcorr_sup = corrcoef(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean, ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean2);
        rgcorrcoeff=rgcorr(1,2);
%         rgcorrcoeff_sup=rgcorr_sup(1,2);
       
        tmplabel2 = ['Ch2, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin2),', zvalmax' num2str(tmpmax2)];   
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin),', zvalmax' num2str(tmpmax)];   
        
        
        if tmpmax>maxthresh && tmpmin<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
            results(r,odornum,1)=tmpmax; results(r,odornum,2)=tmpmin; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
            
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
            results(r,odornum,1)=tmpmax; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
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
            results(r,odornum,2)=tmpmin; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
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
        end
        
        % this is for channel 2 if statements.... change labels to '2'
        % versions.2
        
       if tmpmax2>maxthresh && tmpmin2<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
            results(r,odornum,3) = tmpmax2; results(r,odornum,4) = tmpmin2; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
            if bplotbyodor
                figure(odornum); hold on; %figure(256+odornum); hold on;  <-----use for separate ch2 plots
            else
                figure(r); hold on; %figure(256+r); hold on;  <-----use for separate ch2 plots
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts2,'r--','DisplayName',tmplabel);
            plot(ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts2,'m--','DisplayName',tmplabel);
            figure(512);
            plot(odornum,results(r,odornum,3),'x','Color','r');
            plot(odornum,results(r,odornum,4),'x','Color','r');
        elseif tmpmax2>maxthresh %only excited
            results(r,odornum,3)=tmpmax2;results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts2,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel);
            figure(512);
            plot(odornum,results(r,odornum,3),'x','Color','g');
        elseif tmpmin2<minthresh %only suppressed
            results(r,odornum,4) = tmpmin2; results(r,odornum,5)= rgcorrcoeff; %results(r,odornum,6)=rgcorrcoef_sup;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts2,'Color',tmpcolor,'LineStyle','--','DisplayName',tmplabel);
            figure(512);
            plot(odornum,results(r,odornum,4),'x','Color','b');
        end 
        
    end
    
    clear exbasemean exbasestd exzts
    clear supbasemean supbasestd supzts
    clear tmpmax tmpmin  tmpcolor tmplabel %odornum
end
clear o 
clearvars -except odornums r results 
%clear ch1_ex ch1_sup


vals = [];
ch1ex_fin = isfinite(results(:,:,1));
[row,col] = find(ch1ex_fin==1);
ch1ex_cat = horzcat(row,col);
for r = 1:size(ch1ex_cat,1)
    vals(r) = results(ch1ex_cat(r,1),ch1ex_cat(r,2),1);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch1ex_cat = horzcat(row,col,vals);
ch1ex_sortcat = sortrows(ch1ex_cat);
clearvars row col vals
    else
    end
end

vals = [];
ch1sup_fin = isfinite(results(:,:,2));
[row,col] = find(ch1sup_fin==1);
ch1sup_cat = horzcat(row,col);
for r = 1:size(ch1sup_cat,1)
    vals(r) = results(ch1sup_cat(r,1),ch1sup_cat(r,2),2);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch1sup_cat = horzcat(row,col,vals);
ch1sup_sortcat = sortrows(ch1sup_cat);
clearvars row col vals
    else
    end
end

vals = [];
ch2ex_fin = isfinite(results(:,:,3));
[row,col] = find(ch2ex_fin==1);
ch2ex_cat = horzcat(row,col);
for r = 1:size(ch2ex_cat,1)
    vals(r) = results(ch2ex_cat(r,1),ch2ex_cat(r,2),3);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch2ex_cat = horzcat(row,col,vals);
ch2ex_sortcat = sortrows(ch2ex_cat);
clearvars row col vals
    else
    end
end

vals = [];
ch2sup_fin = isfinite(results(:,:,4));
[row,col] = find(ch2sup_fin==1);
ch2sup_cat = horzcat(row,col);
for r = 1:size(ch2sup_cat,1)
    vals(r) = results(ch2sup_cat(r,1),ch2sup_cat(r,2),4);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch2sup_cat = horzcat(row,col,vals);
ch2sup_sortcat = sortrows(ch2sup_cat);
clearvars row col vals
    else
    end
end

vals = [];
[row,col] = find(results(:,:,5));
chcorr_cat = horzcat(row,col);
for r = 1:size(chcorr_cat,1)
    vals(r) = results(chcorr_cat(r,1),chcorr_cat(r,2),5);
end
vals = vals';
for i=1
    if ~isempty(vals)
chcorr_cat = horzcat(row,col,vals);
chcorr_sortcat = sortrows(chcorr_cat);
clearvars row col vals
    else
    end
end

akmtestmat = zeros(38*136,7);
rois = 38;
for r = 1:38
    for o = 1:136
       akmtestmat = horzcat(r,o,results(r,o,1),results(r,o,2),results(r,o,3),results(r,o,4),results(r,o,5));
    end
end

clearvars -except results ch1ex_sortcat ch1sup_sortcat ch2ex_sortcat ch2sup_sortcat chcorr_sortcat

figure;scatter3(chcorr_sortcat(:,1),chcorr_sortcat(:,2),chcorr_sortcat(:,3),'MarkerFaceColor','k')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');


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