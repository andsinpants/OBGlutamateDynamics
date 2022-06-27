close all; 
% clear;

[exfile,expath] = uigetfile('*.mat','Choose Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ex = load(fullfile(expath,exfile)); %clear expath exfile;

[supfile,suppath] = uigetfile('*.mat','Choose Suppressive File(s) (0.5Hz filter?)','MultiSelect','off');
if ~supfile; return; end
sup = load(fullfile(suppath,supfile)); %clear suppath supfile;

bplotbyodor = 1; %1 means figure# same as odor#; 0 means figure# is roi#
ibase = 1:275; %preodor frames
% iodor = 301:876; %postodor frames %default
iodor = 301:726; %now for 2sec responses
% minthresh = -0; maxthresh = 0; %0.1 to catch nsvalues 
minthresh = -7; maxthresh = 7; %first 6 then 7

%% use for avg set of file
odornums = [];
for r = 1:length(ex.myplotdata.avgfile.roi)
    for o = 1:length(ex.myplotdata.avgfile.roi(r).odor)
        odornums = union(odornums,ex.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
    end
end


%% use for a single file
% % % % odornums = [];  
% % % % for r = 1:length(ex.myplotdata.file.roi)
% % % %     for o = 1:length(ex.myplotdata.file.roi(r).odor)
% % % %         odornums = union(odornums,ex.myplotdata.file.roi(r).odor(o).number,'sorted');
% % % %     end
% % % % end
results = nan(length(ex.myplotdata.avgfile.roi),(odornums(end)),1802);
resultsexc = nan(length(ex.myplotdata.avgfile.roi),(odornums(end)),1802);
resultssup = nan(length(ex.myplotdata.avgfile.roi),(odornums(end)),1802);

figure(256); hold on; %odor number vs significant min/max values
title('Excited and Suppressed');
exbasemeanall = [];
supbasemeanall = [];
for r = 1:length(ex.myplotdata.avgfile.roi)
     for o = 1:length(ex.myplotdata.avgfile.roi(r).odor)
        
        odornum = ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
        
        supbasemean = mean(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
    end
    exbasestd = std(exbasemeanall);  
    supbasestd = std(supbasemeanall);
    
    
    for o = 1:length(ex.myplotdata.avgfile.roi(r).odor)
        odornum = ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasestd = std(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exzts = (ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean)./exbasestd;
        tmpmax = prctile(exzts(iodor),95);
        supbasemean = mean(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasestd = std(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supzts = (sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean)./supbasestd;
        tmpmin = prctile(supzts(iodor),25);
        tmplabel = ['Ch1, roi#' num2str(r) ', odor#' num2str(odornum) ', zvalmin' num2str(tmpmin),', zvalmax' num2str(tmpmax)];      
        
        if tmpmax>maxthresh && tmpmin<minthresh %plot both excitatory and suppressive
            %ch{f}.supMin_5.roi(r).odor(o)<minthresh && ch{f}.exMax_95.roi(r).odor(o)>maxthresh 
            resultsexc(r,odornum,1:length(exzts))=exzts; resultssup(r,odornum,1:length(supzts))=supzts;
            if bplotbyodor
                figure(odornum); hold on;
            else
                figure(r); hold on;
            end
            plot(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'r-','DisplayName',tmplabel);
            plot(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'m-','DisplayName',tmplabel);

            figure(256);
            plot(odornum,resultsexc(r,odornum,1),'x','Color','r');
            plot(odornum,resultssup(r,odornum,2),'x','Color','r');
        elseif tmpmax>maxthresh %only excited
            resultsexc(r,odornum,1:length(exzts))=exzts;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'g'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,exzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,resultsexc(r,odornum,1),'x','Color','g');
        elseif tmpmin<minthresh %only suppressed
            resultssup(r,odornum,1:length(supzts))=supzts;
            if bplotbyodor
                figure(odornum); hold on;
                tmpcolor = 'b'; % use 'g' or myColors(r) here
            else
                figure(r); hold on;
                tmpcolor = myColors(odornum); % use 'g' or myColors(odornum) here
            end
            plot(sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.time,supzts,'Color',tmpcolor,'LineStyle','-','DisplayName',tmplabel);
            figure(256);
            plot(odornum,resultssup(r,odornum,2),'x','Color','b');
        end
    end
end
   
   % clear exbasemean exbasestd exzts
   % clear supbasemean supbasestd supzts
   % clear tmpmax tmpmin odornum tmpcolor tmplabel
%  clear o r
% results(any(~any(results, 1), 3)) = [];
% tmpresults = tmpresults;
%akmtest = find(all(tmpresults(:,:,:)==0,1))
%tmpresults(akmtest) = [];
autoArrangeFigures();

resultszmat_ex = zeros(length(ex.myplotdata.avgfile.roi),size(odornums,2),1802);

for i=1:size(odornums,2)
    resultszmat_ex(:,i,:) = resultsexc(:,odornums(i),:);
end
resultszmat_ex(isnan(resultszmat_ex))=0;

resultszmat_sup = zeros(length(sup.myplotdata.avgfile.roi),size(odornums,2),1802); %%now suppressive

for i=1:size(odornums,2)
    resultszmat_sup(:,i,:) = resultssup(:,odornums(i),:);
end
resultszmat_sup(isnan(resultszmat_sup))=0;

% roi = (1:length(sup.myplotdata.file.roi))';
% 
% for i=1:length(resultszmat_ex)
% newresultszmat_ex(:,:,i) = horzcat(roi,resultszmat_ex(:,:,i));
% end

for i=1:length(resultszmat_ex)
newresultszmat_ex(:,:,i) = vertcat(odornums,resultszmat_ex(:,:,i));
end

for i=1:length(resultszmat_sup)
newresultszmat_sup(:,:,i) = vertcat(odornums,resultszmat_sup(:,:,i));
end

%% now convert to 2D arrays and waterfall plots
for j = 2:size(newresultszmat_ex,1)
    for k = 1:size(newresultszmat_ex,2)
       newerresultszmat_ex=newresultszmat_ex(j,k,:);
       finalresultszmat_ex(j*k,:) = reshape(newerresultszmat_ex,1,1802);
    end
end

ex_resultszmat = finalresultszmat_ex(any(finalresultszmat_ex,2),:);
figure;imagesc(ex_resultszmat);colormap(nawhimar_auto);title('Ch1 Excited')

for j = 2:size(newresultszmat_sup,1)
    for k = 1:size(newresultszmat_sup,2)
       newerresultszmat_sup = newresultszmat_sup(j,k,:);
       finalresultszmat_sup(j*k,:) = reshape(newerresultszmat_sup,1,1802);
    end
end

sup_resultszmat = finalresultszmat_sup(any(finalresultszmat_sup,2),:);
figure;imagesc(sup_resultszmat);colormap(nawhimar_auto);title('Ch1 Suppressed')

times = (0:size(ex_resultszmat,2)-1)/150; 
%% now sort by duration

%% sort by odor latency

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