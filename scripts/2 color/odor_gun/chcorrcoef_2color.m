close all
clear all

[exfile,expath] = uigetfile('*.mat','Choose Channel #1 File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ch1_ex = load(fullfile(expath,exfile)); clear expath exfile;

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
sigresults = nan(length(ch1_ex.myplotdata.avgfile.roi),odornums(end),2);

figure(256); hold on; %odor number vs significant min/max values
title('Channel #1');
exbasemeanall = [];
supbasemeanall = [];
for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
     for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
         
     end
end


% clear o r
ch1ex_fin = isfinite(results(:,:,1));
[row,col] = find(ch1ex_fin==1);

vertcat
% for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
%     for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
%         ch1ex(:,o) = exzts
%     end
%     ch1ex_tot(:,r) = ch1ex
% end

[exfile,expath] = uigetfile('*.mat','Choose Channel #2 File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ch2_ex = load(fullfile(expath,exfile)); clear expath exfile;

bplotbyodor = 1; %1 means figure# same as odor#; 0 means figure# is roi#
ibase = 1:275; %preodor frames
iodor = 301:825; %postodor frames
minthresh = -7; maxthresh = 7; %first 6 then 7

odornums = [];
for r = 1:length(ch2_ex.myplotdata.avgfile.roi)
    for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
        odornums = union(odornums,ch2_ex.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
    end
end
results = nan(length(ch2_ex.myplotdata.avgfile.roi),odornums(end),2);

figure(256); hold on; %odor number vs significant min/max values
title('Channel #2');
exbasemeanall = [];
supbasemeanall = [];
for r = 1:length(ch2_ex.myplotdata.avgfile.roi)
     for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
        
        supbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
    end
    exbasestd = std(exbasemeanall);  
    supbasestd = std(supbasemeanall);
    
    
    for o = 1:length(ch2_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch2_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        %exbasestd = std(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exzts = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean)./exbasestd;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
       %supbasestd = std(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean)./supbasestd;
        tmpmin = prctile(supzts(iodor),25);
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin),', zvalmax' num2str(tmpmax)];      
        if tmpmax>maxthresh && tmpmin<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
            results(r,odornum,1)=tmpmax; results(r,odornum,2)=tmpmin;
            if bplotbyodor
                figure(odornum); hold on;
            else
                figure(r); hold on;
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'r-','DisplayName',tmplabel);
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'m-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','r');
            plot(odornum,results(r,odornum,2),'x','Color','r');
            bimodalch2=ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time;
        elseif tmpmax>maxthresh %only excited
            results(r,odornum,1)=tmpmax;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,1),'x','Color','k');
            excite_ch2 = (ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time); %%
        elseif tmpmin<minthresh %only suppressed
            results(r,odornum,2)=tmpmin;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ch2_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,results(r,odornum,2),'x','Color','c');
            suppress_ch2 = ch2_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time;
        end
    end
    clear exbasemean exbasestd exzts
    clear supbasemean supbasestd supzts
    clear tmpmax tmpmin odornum tmpcolor tmplabel
end
clear o r
