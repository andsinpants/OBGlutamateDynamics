%%% sort out deep MTC and superfical %%%
% close all; clear all;
%load('tc/TCsomatanoCCK.mat');
% % load('cck/TCsomataCCK.mat'); %original matrix
% load('tbt100_00003_ethylbutyrate_avgSTAdata_struct_1to3sec_bgsub')
myplotdata=A;
[rows num_of_rois ] = size(myplotdata);%%% calculate the num of ROI per file

odor=[];exp=[];name=[];SVch1=[];SVch2=[];num_of_roisCh1=[];
num_of_roisCh2=[];roisYCh1= [];roisYCh2=[];roilableCh2=[];
roilableCh1=[];roisXCh1=[];roisXCh2=[];expCh1=[];expCh1a=[];
roia=[];roiCh1=[];odorch1=[];mice=[];name=[];Ch1Odor2=[];
resp_value2=[];resp_value=[];Ch1rig=[];Ch1rigind=[];

liqdilch1=[];
SVch1_sv=[];
for i=1:num_of_rois
    roi=length(myplotdata{1,i}.roi);
    for ii=1:roi
        %roia=myplotdata(1,i).roi(1,ii).number; %
        %%roiCh1=[roiCh1;roia];
        %expCh1=myplotdata(1,i).exp; expCh1a=[expCh1a;expCh1];
        %mouse=myplotdata(1,i).mouse; mice=[mice;mouse];
        
        roisYindCh1=myplotdata{1,i}.roi(ii).odor.avgtrial.series';
        roisYCh1=catpad(1,roisYCh1,roisYindCh1);
        
        roisXindCh1=myplotdata{1,i}.roi(ii).odor.avgtrial.time;
        roisXCh1=catpad(1,roisXCh1,roisXindCh1);
        
        qq=find(roisXindCh1>1.15&roisXindCh1<1.5);
        qqq=roisYindCh1(qq); resp_value2=mean(qqq);
        resp_value=[resp_value;resp_value2];
        
        %roi1lab=myplotdata{1,i}.roi(ii).number;
        roi1lab=0;
        roilableCh1=[roilableCh1;roi1lab];
        % %             %                 SatvapCh1=myplotdata(1,i).liqcal; %
        %                 SVch1=[SVch1;SatvapCh1]; %             % %
        %                 SatvapCh1=myplotdata(1,i).dil; %             %
        % liqdilch1=[SVch1;SatvapCh1];
        
        %SatvapCh1=myplotdata(1,i).SV;%recomment
        %SVch1_sv=[SVch1_sv;SatvapCh1];%recomment
        
        name1=myplotdata{1,1}.name;%recomment
        name=[name;cellstr(name1)];%recomment
        %odor22=myplotdata(1,i).odornum;%recomment
        %Ch1Odor2=[Ch1Odor2;odor22];%recomment
        %Ch1rigind2=myplotdata(1,i).rig;
        %Ch1rigind=[Ch1rigind;Ch1rigind2];
        
        
        
    end %%%% uncomment for TC cells/other cells
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% fit Ch1 data to a double gaussian%%%
min98peak=[];min98peakXt=[];minCh198=[];minTCh198=[];

peakTCh115Exp=[];peakTCh115Odor=[];peakTCh115SV=[];peakTCh198Exp=[];
peakTCh198Odor=[];peakTCh198=[];SV=[];minTCh198Exp=[];minTCh198Odor=[];
minTCh198SV=[];peakTCh115Exp=[];peakTCh115Odor=[];peakTCh115SV=[];

latencyY=[];latencyX=[];latencyYinh=[];latencyXinh=[];roi_indinh=[];

Ch1Expinh=[];Ch1SVinh=[];roiCh1inh=[];

roi_ind=[];roi_indin=[];earlyXlat=[];earlyYlat=[];

roiname=[];roinameinh=[];Ch1Exp=[];Ch1Odor=[];Ch1SV=[];mice_min=[];
mice_ex=[];mice_ex198=[];

ylat3graph=[];xlat3graph=[];ylat4graph=[];

Consect5lat=[];initLat=[];inhcells=[];inhcellsrois=[];yinh=[];peak_LatY=[];peak_LatX=[];min_LatY=[]; min_LatX=[];

latencyMininh=[];liq=[];svonly=[];liqinh=[];
svonlyinh=[];Ch1Odorinh=[];Ch1Odor=[];
decayminCh150=[];decayminTCh150=[]; minCh150=[];minTCh150=[];  minCh1=[];minTCh1=[];
RiseslopesCh1=[];halfRiseslopesCh1=[];peakCh1=[];peakTCh1=[];peakCh175=[];peakTCh175=[];peakCh125=[];peakTCh125=[];halfpeakCh1=[];halfpeakTch1=[];PolyfitYCh1=[];halfSlopeCh1=[];linregCh1=[];peakCh198=[];peakTCh198=[];peakCh115=[];peakTCh115=[];minrespCh1=[];minrespTCh1=[];minCh198=[];minTCh198=[];decayCh150=[];decayTCh150=[];peakCh150=[];peakTCh150=[];halfwidth=[];peakCh50=[];peakTCh50=[];decayCh50=[];decayTCh50=[];
yinflect=[];
[num_of_roisCh1 col] = size(roisYCh1);%%% calculate the num of ROI per Ch2
Ch1mice=[];

fastsuppression=[];
delayedexcit=[];

halfwidth_pca=[]; halfwidthmin_pca=[];
halfwidthmin_pca1=[];halfwidth_pca1=[];


respwindow=[];
inc=[];
inc95=[];
quantin=[];
quantinmin=[];
dec=[];
min05=[];
odoresp=[];
odorebound=[];
rebound=[];


for i=1:num_of_roisCh1
    x=roisXCh1(i,:); roisX=roisXCh1(i,:); time=roisXCh1(i,:);
    y=roisYCh1(i,:);
    %
    %     [miny, minIn]=min(y);
    %     minTime=x(minIn);
    %     minrespCh1=[minrespCh1;miny];
    %     minrespTCh1=[minrespTCh1;minTime];
    %
    %     %     %USE to save variabls to calculate group STATS
    %     %     exp2=expCh1(i);
    %     %     odor2=odorch1(i);
    %     %     sv2=SVch1(i);
    %     %     Ch1Exp=[Ch1Exp;exp2];
    %     %     Ch1Odor=[Ch1Odor;odor2];
    %     %     Ch1SV=[Ch1SV;sv2];
    %
    %     if minTime>0
    %         if miny<-0.01
    %             %Calculate 98%peak and 98%peakt time
    %             min98peakInd=find(y==1*miny);
    %             min98peakXt=x(min98peakInd);
    %             min98peak=y(min98peakInd);
    %             [min98peakInd2]=find(min98peakXt>1,1);
    %             min98peak=min98peak(min98peakInd2);
    %             min98peakXt=min98peakXt(min98peakInd2);
    %             %linreg2=find(x==min98peakXt);
    %             minCh198=[minCh198;min98peak];
    %             minTCh198=[minTCh198;min98peakXt];
    %
    %             %         %STATS
    %             % %         exp2=expCh1(i);
    %             % %         odor2=odorch1(i);
    %             % %         sv2=SVch1(i);
    %             % %         minTCh198Exp=[minTCh198Exp;exp2];
    %             % %         minTCh198Odor=[minTCh198Odor;odor2];
    %             % %         minTCh198SV=[minTCh198SV;sv2];
    %
    %             %mice_min2=mice(i);
    %
    %         end
    %         % mice_min=[mice_min,mice_min2];
    %
    %     end
    
    
    %     %% Calculate Derivate of MaxSlope
    %     f1 = diff(y);
    %     [Maxslope_onset Max_slope_dfind]=max(f1);
    %     xinflect=x(Max_slope_dfind);
    %     yinflect1=max(f1);
    %     yinflect=[yinflect;yinflect1];
    
    
    %%%%%SET stimulus window and calculate basic stats
    stimoff=find(x>2,1);%x is the x value timeseries data
    stimon=find(x>1.005,1);
    stimoffearly=find(x>1.75,1);
    prestim=(1:stimon);
    
    respwindow=y(stimon:stimoff);
    offrespwindow=y(stimoffearly:stimoff);
    if mean(y(stimon:stimoffearly))>0
        inc=max(respwindow);
        inc95=quantile(respwindow,0.95);
        quantind=find(respwindow>inc95);
        max95=mean(respwindow(quantind));
        resp=max95;
        rebound=NaN;
        
    elseif mean(y(stimon:stimoffearly))<0
        dec=min(respwindow);
        inc5=quantile(respwindow, 0.05);
        quantindmin=find(respwindow<inc5);
        min05=mean(respwindow(quantindmin));
        resp=min05;
        if mean(y(stimoffearly:stimoff))>0
            inc=max(offrespwindow);
            maxind=find(offrespwindow==inc);
            minin=find(offrespwindow==dec);
            if maxind>minin
                
                quantindmaxin=find(respwindow>inc95);
                max95maxin=mean(respwindow(quantind));
                if max95maxin>0
                    rebound=max95maxin;
                end
            end
            if maxind<minin
                rebound=NaN;
            end
        end
        
        
        
    elseif mean(y(stimon:stimoffearly))==0
        dec=min(respwindow);
        inc50=quantile(respwindow, 0.05);
        quantindmin0=find(respwindow<inc50);
        min050=mean(respwindow(quantindmin0));
        resp=min050;
        if mean(y(stimoffearly:stimoff))>0
            inc=max(offrespwindow);
            maxind=find(offrespwindow==inc);
            minin=find(offrespwindow==dec);
            if maxind>minin
                
                quantindmaxin=find(respwindow>inc95);
                max95maxin=mean(respwindow(quantind));
                if max95maxin>0
                    rebound=max95maxin;
                end
            end
            
        end
        if maxind<minin
            rebound=NaN;
        end
        
    end
    
    
    odoresp=[odoresp;resp];
    odorebound=[odorebound;rebound];
    
    
    %%%% Use threshold to set response onset
    
    prestimstd=std(y(1:stimon));%y is the y value timeseries data
    prestimean=mean(y(1:stimon));
    
    
    %%% z-score each value by prestim period mean ans STD %%%
    ylat2=y-prestimean;
    ylat=ylat2./(prestimstd*1);
%     
%     %     %%% apply moving average, horiszontal (2), every 10 unit %%%%% no
%     %     longer being use, replaced with consecutive thres
%     %     ylat=movingmean(ylat3,10,2);
%     
%     
%     ylat4graph=[ylat4graph;ylat];
%     xlat4graph=[xlat3graph;x];
%     
%     poststimY=ylat(stimon:stimoff);
%     
%     prestimperiod=ylat(1:stimon);
%     figure()
    
    prestimperiod=y(1:stimon);
    prestimperiodx=x(1:stimon);
%     plot(prestimperiodx,prestimperiod,'LineWidth',2)
    % Fit line to data using polyfit
    c = polyfit(prestimperiodx,prestimperiod,1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,prestimperiodx);
    % Add trend line to plot
%     hold on
%     plot(prestimperiodx,y_est,'r--','LineWidth',2)
%     hold off
    
    linesub_prestimperiod=prestimperiod-y_est;
    prestimstdlat=std(linesub_prestimperiod);
    lastptfit=y_est(end);
    
%     %%recalc STD on new trace zscored data
%     prestimstdlat=std(ylat(1:stimon));%should be 1
    
    
    earlystim=y(stimon:stimoffearly);
    earlystimpost=y(stimoffearly:stimoff);
    
    poststimrealy=y(stimon:end);%% save real y post stim for graphing
    poststimX=x(stimon:end); %collect x's post stim
    
     poststimY=y(stimon:end);
    thresstd=4;lengthwindow=40; 
    %%%5 consecutive values above the threshold
    aboveThreshold = (poststimY > lastptfit+prestimstdlat*thresstd);
    thresholdChange = [aboveThreshold(1) diff(aboveThreshold)];
    thresholdChange(thresholdChange==-1) = 0;
    spanLocs = cumsum(thresholdChange);
    spanLocs(~aboveThreshold) = 0;
    aboveThresholdIndex = find(aboveThreshold==1);
    notConsecutiveIndex = [true diff(aboveThresholdIndex) ~= 1];
    sumIndex = cumsum(notConsecutiveIndex);
    spanLength = histc(sumIndex, 1:sumIndex(end));
    goodSpans = find(spanLength>=lengthwindow);%orig>=5
    allInSpans = find(ismember(spanLocs, goodSpans));
    latencysd_in=min(allInSpans);
    
    %save index of 5 consective value above a threshold
    Consect5lat=[Consect5lat;latencysd_in];
    
    %%%save index of initial latency value
    latency2sd_in=find(poststimY>lastptfit+prestimstdlat*thresstd,1);%%% mult prestimstdlat by # of SD threshold
    initLat=[initLat;latency2sd_in];
    
    
    latencyY2=poststimrealy(latencysd_in); %%set latency value to real data pts for graphs
    latencyX2=poststimX(latencysd_in);
    maxpoststim=max(poststimrealy);
    %    if latencyX2>1.005
    %         if maxpoststim>0.15
    if sum(ismember(y,latencyY2))>0
        latencyY=[latencyY;latencyY2];
        latencyX=[latencyX;latencyX2];%save latency values
        
        roi_ind=[roi_ind;i];% save max values associated with onset lat. caclucated
        peak_LatY=[peak_LatY;max(poststimrealy)];
        
        peak_LatX=[peak_LatX;x(find(y==max(poststimrealy)))];
        
        
        
        
%         %%%latency specific experiment variables saved
%         Ch1Exp=expCh1a(roi_ind);
%         Ch1SV=SVch1_sv(roi_ind);
        roiCh1=roilableCh1(roi_ind);
        roiname=name(roi_ind);
        %             liq=liqdilch1(roi_ind);
        %             svonly=SVch1_sv(roi_ind);
%         Ch1Odor=Ch1Odor2(roi_ind);
%         Ch1mice=mice(roi_ind);
        %Ch1rig=Ch1rigind(roi_ind);
        
%         x=roisXCh1(i,:); roisX=roisXCh1(i,:); time=roisXCh1(i,:);
%         y=roisYCh1(i,:);
        
%         h=figure();
%         plot(x,y)
%         hold on
%         scatter(latencyX2,latencyY2,'ro')
%         saveas(h,sprintf('MC_Ex6_40Pts4SD/FIG%d.png',i));
        
        %Calculate peak and peak time, Halfwidth
        ysmooth2 = smoothdata(y,'gaussian',50)
        [peak, peakIn]=max(ysmooth2);
        peakTime=x(peakIn);
        peakCh1=[peakCh1;peak];
        peakTCh1=[peakTCh1;peakTime];
        Per150peakInd=find(ysmooth2>.5*peak);
        Per150peakXt=x(Per150peakInd);
        Per150peak=ysmooth2(Per150peakInd);
        [Per150peakInd2]=find(Per150peakXt>1,1);
        Per150peak=Per150peak(Per150peakInd2);
        Per150peakXt=Per150peakXt(Per150peakInd2);
        %linreg1=find(x==Per150peakXt);
        peakCh150=[peakCh150;Per150peak];
        peakTCh150=[peakTCh150;Per150peakXt];
        
        
        %Calculate 50%decay and 50%decayt time
        Per50decayInd=find(ysmooth2>(.5*peak));
        Per50decayIndq=Per50decayInd(end);
        Per50decayXt=x(Per50decayIndq);
        Per50decay=ysmooth2(Per50decayIndq);
        decayCh150=[decayCh150;Per50decay];
        decayTCh150=[decayTCh150;Per50decayXt];
        
        halfwidth_pca1=Per50decayXt-Per150peakXt;
        halfwidth_pca=[halfwidth_pca;halfwidth_pca1];
        
        h=figure();
        plot(x,y)
        hold on
        plot(x,ysmooth2)
%         scatter(peak_LatX,peak_LatY,'go')
%         hold off
        scatter(latencyX2,latencyY2,'ro')
        
%         
%         h=figure();
%         plot(x,y)
%         hold on
%         plot(x,ysmooth2,'r')
%         hold on
%         scatter(Per50decayXt,Per50decay,'r*')
%         hold on
%         scatter(Per150peakXt,Per150peak,'r*')
%         hold on 
%         scatter(peakTime,peak,'g*')
%         saveas(h,sprintf('MC_Ex6_100gauss/FIG%d.png',i));
        
        %         end
        %
        %     end
    end
    
    
    
    %%%% NOW repeat for inhibitory events
    %%%5 consecutive values above the threshold
    
            %%% z-score each value by prestim period mean ans STD %%%
    ylat2=y-prestimean;
    ylat=ylat2./(prestimstd*1);
    
    %     %%% apply moving average, horiszontal (2), every 10 unit %%%%% no
    %     longer being use, replaced with consecutive thres
    %     ylat=movingmean(ylat3,10,2);
    
    
    ylat4graph=[ylat4graph;ylat];
    xlat4graph=[xlat3graph;x];
    
    poststimY2=ylat(stimon:stimoff);
    
%     prestimperiod=ylat(1:stimon);
    
    
    
        %%recalc STD on new trace zscored data
    prestimstdlat=std(ylat(1:stimon));%should be 1
    
    minbaseline=min(ylat(1:stimon));
    maxbaseline=max(ylat(1:stimon));
    mintotal=min(ylat);
    maxtotal=max(ylat);
    belowThreshold = (poststimY2 < -prestimstdlat*3);%% originally set to 1STD
    thresholdChangeinh = [belowThreshold(1) diff(belowThreshold)];
    thresholdChangeinh(thresholdChangeinh==-1) = 0;
    spanLocs = cumsum(thresholdChangeinh);
    spanLocs(~belowThreshold) = 0;
    belowThresholdIndex = find(belowThreshold==1);
    notConsecutiveIndex = [true diff(belowThresholdIndex) ~= 1];
    sumIndex = cumsum(notConsecutiveIndex);
    spanLength = histc(sumIndex, 1:sumIndex(end));
    goodSpans = find(spanLength>=lengthwindow);
    allInSpansbelow = find(ismember(spanLocs, goodSpans));
    latencysd_inINHB=min(allInSpansbelow);
    
    
    
    fastsuppression=[fastsuppression;mean(earlystim)];
    delayedexcit=[delayedexcit;mean(earlystimpost)];
    
   
        
        %if abs(maxbaseline-abs(minbaseline))<1.5
        latencyY2inh=poststimrealy(latencysd_inINHB); %%set latency value to real data pts for graphs
        latencyX2inh=poststimX(latencysd_inINHB);
        minpoststim=min(poststimrealy);
        if latencyX2inh>1.005
            if latencyX2inh<3
                if minpoststim<0
                    
                    latencyYinh=[latencyYinh;latencyY2inh];
                    latencyXinh=[latencyXinh;latencyX2inh];%save latency values
                    latencyMininh=[latencyMininh;minpoststim];%save min 'peak'
                    roi_indinh=[roi_indinh;i];
                    
                % save max values associated with onset lat. caclucated
                
                
                min_LatY=[min_LatY;min(poststimrealy)];
                min_LatX=[min_LatX;x(find(y==min(poststimrealy)))];
                
%                 %%%latency specific experiment variables saved
%                 Ch1Expinh=expCh1a(roi_indinh);
%                 Ch1SVinh=SVch1_sv(roi_indinh);
                roiCh1inh=roilableCh1(roi_indinh);
                roinameinh=name(roi_indinh);
%                 Ch1Odorinh=Ch1Odor2(roi_indinh);
%                 
                
                %Calculate suppressions and time
                [minsupp, suppIn]=min(y);
                minTime=x(suppIn);
                minCh1=[minCh1;minsupp];
                minTCh1=[minTCh1;minTime];
                Per150minInd=find(y<.333*minsupp);
                Per150minXt=x(Per150minInd);
                Per150min=y(Per150minInd);
                [Per150minInd2]=find(Per150minXt>1,1);
                Per150min=Per150min(Per150minInd2);
                Per150minXt=Per150minXt(Per150minInd2);
                %linreg1=find(x==Per150minXt);
                minCh150=[minCh150;Per150min];
                minTCh150=[minTCh150;Per150minXt];
                
                %Calculate 50%decay and 50%decayt time
                Min50decayInd=find(y<(.333*minsupp));
                Min50decayIndq=Min50decayInd(end);
                
                Min50decayXt=x(Min50decayIndq);
                Min50decay=y(Min50decayIndq);
                decayminCh150=[decayminCh150;Min50decay];
                decayminTCh150=[decayminTCh150;Min50decayXt];
                
                
                halfwidthmin_pca1=Min50decayXt-Per150minXt;
                halfwidthmin_pca=[halfwidthmin_pca;halfwidthmin_pca1];
              
        end
        
     end
end

end

%end
autoArrangeFigures();

halfwidthmin=[]%decayminTCh150-minTCh150;

%Calculate half width
halfwidth=decayTCh150-peakTCh150;

% % Calculate Rise Time
risetime=[]%peakTCh198-peakTCh115;



risetime_normPOP=[];

for i=min(mice_ex):max(mice_ex)
    indM=find(mice_ex==i);
    M1=zscore(risetime(indM));
    risetime_normPOP=[risetime_normPOP;M1];
end


peakTCh115_normPOP=[];
for i=min(mice_ex):max(mice_ex)
    indM=find(mice_ex==i);
    M1=zscore(peakCh115(indM));
    peakTCh115_normPOP=[peakTCh115_normPOP;M1];
end


%%%% Normalize by the mean of the population by mouse

% risetime_normPOP=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=risetime(indM)- mean(risetime(indM));
%     risetime_normPOP=[risetime_normPOP;risetime_normpop];
% end


halfwidth_normPOP=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=halfwidth(indM)/mean(halfwidth(indM));
%     halfwidth_normPOP=[halfwidth_normPOP;risetime_normpop];
% end
%
%
peakCh115_normPOP=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=peakCh115(indM)- mean(peakCh115(indM));
%     peakCh115_normPOP=[peakCh115_normPOP;risetime_normpop];
% end

% peakTCh115_normPOP=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=peakTCh115(indM)-mean(peakTCh115(indM));
%     peakTCh115_normPOP=[peakTCh115_normPOP;risetime_normpop];
% end


peakTCh198_normPOP=[];
peakTCh198_normPOP2=[];
% mice_rise=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=peakTCh198(indM)-mean(peakTCh198(indM));
%     %risetime2_normpop=peakTCh198(indM).-mean(peakTCh198(indM));
%     peakTCh198_normPOP=[peakTCh198_normPOP;risetime_normpop];
%     %peakTCh198_normPOP2=[peakTCh198_normPOP2;risetime2_normpop];
%     mice_rise=[mice_rise;mean(peakTCh198(indM))];
% end


peakCh198_normPOP=[];
% for i=1:max(mouse)
%     indM=find(mice_ex==i);
%     risetime_normpop=peakCh198(indM)-mean(peakCh198(indM));
%     peakCh198_normPOP=[peakCh198_normPOP;risetime_normpop];
% end

%
% %%%Z score
% risetime_normPOP=zscore(risetime);
% halfwidth_normPOP=zscore(halfwidth);
%
% peakCh115_normPOP=zscore(peakCh115);
% peakTCh115_normPOP=zscore(peakTCh115);
%
% peakTCh198_normPOP=zscore(peakTCh198);
% peakCh198_normPOP=zscore(mean(peakCh198);

% %% divid mean
% risetime_normPOP=risetime./(mean(risetime));
% halfwidth_normPOP=halfwidth./(mean(halfwidth));
%
% peakCh115_normPOP=peakCh115./(mean(peakCh115));
% peakTCh115_normPOP=peakTCh115./(mean(peakTCh115));
%
% peakTCh198_normPOP=peakTCh198./(mean(peakTCh198));
% peakCh198_normPOP=peakCh198./(mean(peakCh198));
%
% %% Subtract mean%
% risetime_normPOP=risetime-(mean(risetime));
% halfwidth_normPOP=halfwidth-(mean(halfwidth));
%
% peakCh115_normPOP=peakCh115-(mean(peakCh115));
% peakTCh115_normPOP=peakTCh115-(mean(peakTCh115));
%
% peakTCh198_normPOP=peakTCh198-(mean(peakTCh198));
% peakCh198_normPOP=peakCh198-(mean(peakCh198));


fastsup=[];
for i=1:length(delayedexcit)
    if delayedexcit(i)>0  && fastsuppression(i)<0
        qtst=1;
    else
        qtst=0;
    end
    fastsup=[fastsup;qtst];
end





%TBTsoma_super,TBTsoma
%TCsoma, MCsoma, OMP,
%CCKsoma2,GADsoma,PCDsoma,THYsoma,OMP35,DATsoma,
%MTCsoma_thyTC,MTCsoma_pcdTC,MTCsoma_tbtTC,
%MTCsoma_thyMC,MTCsoma_tbtMC,MTCsoma_pcdMC
%redpg,greentc
filesv=(num2str(['TCnoCCK']));
%% Write all matrices of intests with data file label%%%
eval([ 'risetime_normPOP' num2str(filesv) ' = risetime_normPOP;' ]);
eval([ 'halfwidth_normPOP' num2str(filesv) ' = halfwidth_normPOP;' ]);
eval([ 'peakCh115_normPOP' num2str(filesv) ' = peakCh115_normPOP;' ]);
eval([ 'peakTCh115_normPOP' num2str(filesv) ' = peakTCh115_normPOP;' ]);
eval([ 'peakTCh198_normPOP' num2str(filesv) ' = peakTCh198_normPOP;' ]);
eval(['peakCh198_normPOP' num2str(filesv) ' = peakCh198_normPOP;' ]);
eval(['name' num2str(filesv) ' = name;' ]);



eval([ 'roinameinh' num2str(filesv) ' = roinameinh;' ]);
eval([ 'roiname' num2str(filesv) ' = roiname;' ]);
eval([ 'Ch1Expinh' num2str(filesv) ' = Ch1Expinh;' ]);
eval([ 'Ch1SVinh' num2str(filesv) ' = Ch1SVinh;' ]);
eval([ 'roiCh1inh' num2str(filesv) ' = roiCh1inh;' ]);
eval([ 'Ch1Odorinh' num2str(filesv) ' = Ch1Odorinh;' ]);
eval([ 'Ch1Odor2' num2str(filesv) ' = Ch1Odor2;' ]);

eval([ 'fastsuppression' num2str(filesv) ' = fastsuppression;' ]);
eval([ 'delayedexcit' num2str(filesv) ' = delayedexcit;' ]);
eval([ 'fastsup' num2str(filesv) ' = fastsup;' ]);


eval([ 'halfwidthmin' num2str(filesv) ' = halfwidthmin;' ]);


eval([ 'halfwidth_pca' num2str(filesv) ' = halfwidth_pca;' ]);
eval([ 'halfwidthmin_pca' num2str(filesv) ' = halfwidthmin_pca;' ]);



eval([ 'Ch1rig' num2str(filesv) ' = Ch1rig;' ]);
eval([ 'Ch1mice' num2str(filesv) ' = Ch1mice;' ]);

eval([ ' resp_value' num2str(filesv) ' =  resp_value;' ]);
eval([ 'expCh1a' num2str(filesv) ' = expCh1a;' ]);% all exp for resp value
eval([ 'roilableCh1' num2str(filesv) ' = roilableCh1;' ]);% all exp for resp value
eval([ 'SVch1_sv' num2str(filesv) ' = SVch1_sv;' ]);% all exp for resp value


eval([ 'latencyMininh' num2str(filesv) ' = latencyMininh;' ]);
eval([ 'peak_LatY' num2str(filesv) ' = peak_LatY;' ]);
eval(['peak_LatX' num2str(filesv) ' = peak_LatX;' ]);

eval([ 'latencyY' num2str(filesv) ' = latencyY;' ]);
eval(['latencyX' num2str(filesv) ' = latencyX;' ]);
eval([ 'latencyYinh' num2str(filesv) ' = latencyYinh;' ]);
eval(['latencyXinh' num2str(filesv) ' = latencyXinh;' ]);

eval([ 'liqinh' num2str(filesv) ' = liqinh;' ]);
eval(['svonlyinh' num2str(filesv) ' = svonlyinh;' ]);
eval([ 'liq' num2str(filesv) ' = liq;' ]);
eval(['svonly' num2str(filesv) ' = svonly;' ]);


eval(['SVch1_sv' num2str(filesv) ' = SVch1_sv;' ]);


eval([ 'minrespCh1' num2str(filesv) ' = minrespCh1;' ]);
eval([ 'minrespTCh1' num2str(filesv) ' = minrespTCh1;' ]);

eval([ 'minrespCh1' num2str(filesv) ' = minrespCh1;' ]);
eval([ 'minrespTCh1' num2str(filesv) ' = minrespTCh1;' ]);
eval([ 'minCh198' num2str(filesv) ' = minCh198;' ]);
eval([ 'minTCh198' num2str(filesv) ' = minTCh198;' ]);
eval([ 'peakCh125' num2str(filesv) ' = peakCh125;' ]);
eval(['peakTCh125' num2str(filesv) ' = peakTCh125;' ]);
eval(['peakCh175' num2str(filesv) ' = peakCh175;' ]);
eval(['peakTCh175' num2str(filesv) ' = peakTCh175;' ]);
eval(['peakCh1' num2str(filesv) ' = peakCh1;' ]);
eval(['peakTCh1' num2str(filesv) ' = peakTCh1;' ]);
eval(['peakCh115' num2str(filesv) ' = peakCh115;' ]);
eval(['peakTCh115' num2str(filesv) ' = peakTCh115;' ]);
eval(['yinflect' num2str(filesv) ' = yinflect;' ]);
eval(['linregCh1' num2str(filesv)  ' = linregCh1;' ]);%% first value is the slope, second value is the y-intercept
eval(['decayCh150' num2str(filesv)  ' = decayCh150;' ]);
eval(['decayTCh150' num2str(filesv)  ' = decayTCh150;' ]);
eval(['peakCh150' num2str(filesv)  ' = peakCh150;' ]);
eval(['peakTCh150' num2str(filesv)  ' = peakTCh150;' ]);
eval(['halfwidth' num2str(filesv)  ' = halfwidth;' ]);
eval(['peakCh198' num2str(filesv)  ' = peakCh198;' ]);
eval(['peakTCh198' num2str(filesv)  ' = peakTCh198;' ]);
eval(['risetime' num2str(filesv)  ' = risetime;' ]);

eval(['peakTCh115Exp' num2str(filesv)  ' = peakTCh115Exp;' ]);
eval(['peakTCh115Odor' num2str(filesv)  ' = peakTCh115Odor;' ]);
eval(['peakTCh115SV' num2str(filesv)  ' = peakTCh115SV;' ]);

eval(['peakTCh198Exp' num2str(filesv)  ' = peakTCh198Exp;' ]);
eval(['peakTCh198Odor' num2str(filesv)  ' = peakTCh198Odor;' ]);
% eval(['peakTCh198SV' num2str(filesv)  ' = peakTCh198SV;' ]);

eval(['Ch1Exp' num2str(filesv)  ' = Ch1Exp;' ]);
eval(['Ch1Odor' num2str(filesv)  ' = Ch1Odor;' ]);
eval(['Ch1SV' num2str(filesv)  ' = Ch1SV;' ]);
eval(['roiCh1' num2str(filesv)  ' = roiCh1;' ]);
eval([ 'Ch1Odor' num2str(filesv) ' = Ch1Odor;' ]);

eval([ 'risetime_normPOP' num2str(filesv) ' = risetime_normPOP;' ]);
eval([ 'halfwidth_normPOP' num2str(filesv) ' = halfwidth_normPOP;' ]);
eval([ 'peakCh115_normPOP' num2str(filesv) ' = peakCh115_normPOP;' ]);
eval([ 'peakTCh115_normPOP' num2str(filesv) ' = peakTCh115_normPOP;' ]);
eval([ 'peakTCh198_normPOP' num2str(filesv) ' = peakTCh198_normPOP;' ]);
eval(['peakCh198_normPOP' num2str(filesv) ' = peakCh198_normPOP;' ]);

eval([ 'odoresp' num2str(filesv) ' = odoresp;' ]);
eval(['odorebound' num2str(filesv) ' = odorebound;' ]);

eval([ 'min_LatY' num2str(filesv) ' = min_LatY;' ]);
eval(['min_LatX' num2str(filesv) ' = min_LatX;' ]);



%%%%%above retired %%%%%%%%%

%OMP_popnorm_VARSExInONLYEX
%TCsoma2_popnorm_VARSExInONLYEX
%MCsoma2_popnorm_VARSExInONLYEX
%CCK_popnorm_VARSExInONLYEX
%GADsoma_popnorm_VARSExInONLYEX
%DATsoma_popnorm_VARSExInONLYEX



%TBTsomaSuper_popnorm_VARSExInONLYEX') %superficial only
%MTCsoma_pcdTC_popnorm_VARSExInONLYEX')
%MTCsoma_thyTC_popnorm_VARSExInONLYEX')


%MTCsoma_tbtMC_popnorm_VARSExInONLYEX')%mitral only
%MTCsoma_pcdMC_popnorm_VARSExInONLYEX')
%MTCsoma_thyMC_popnorm_VARSExInONLYEX')

%'redpg_popnorm_VARSExInONLYEX'
%'greentc_popnorm_VARSExInONLYEX
% GADsaremove_popnorm_VARSExInONLYEX
%save('OMP_popnorm_VARSExInONLYEX',['resp_value' num2str(filesv)],['Ch1Odorinh' num2str(filesv)],['Ch1Odor' num2str(filesv)],['liqinh' num2str(filesv)],['liq' num2str(filesv)],['svonlyinh' num2str(filesv)],['svonly' num2str(filesv)],['roinameinh' num2str(filesv)],['roiname' num2str(filesv)],['Ch1Expinh' num2str(filesv)],['Ch1SVinh' num2str(filesv)],['roiCh1inh' num2str(filesv)],['latencyMininh' num2str(filesv)],['peak_LatX' num2str(filesv)],['roiCh1' num2str(filesv)],['peak_LatY' num2str(filesv)],['latencyYinh' num2str(filesv)],['latencyXinh' num2str(filesv)],['latencyY' num2str(filesv)],['latencyX' num2str(filesv)],['peakTCh115Exp' num2str(filesv)],['peakTCh115Odor' num2str(filesv)],['peakTCh115SV' num2str(filesv)],['peakTCh198Exp' num2str(filesv)],['peakTCh198Odor' num2str(filesv)],['peakTCh198SV' num2str(filesv)],['Ch1Exp' num2str(filesv)],['Ch1Odor' num2str(filesv)],['Ch1SV' num2str(filesv)],['risetime' num2str(filesv)],['minrespCh1' num2str(filesv)],[ 'minrespTCh1' num2str(filesv)],['minCh198' num2str(filesv)],[ 'minTCh198' num2str(filesv)], ['peakCh125' num2str(filesv)],['peakTCh125' num2str(filesv)],['peakTCh175' num2str(filesv)],['peakCh1' num2str(filesv)],['peakTCh1' num2str(filesv)], ['peakCh115' num2str(filesv)],['peakTCh115' num2str(filesv)],['yinflect' num2str(filesv)],['linregCh1' num2str(filesv)],['decayCh150' num2str(filesv)],['decayTCh150' num2str(filesv)],['peakCh150' num2str(filesv)],['peakTCh150' num2str(filesv)],['halfwidth' num2str(filesv)],['peakCh198' num2str(filesv)],['peakTCh198' num2str(filesv)],['risetime_normPOP' num2str(filesv)],['halfwidth_normPOP' num2str(filesv)],['peakCh115_normPOP' num2str(filesv)],['peakTCh115_normPOP' num2str(filesv)],['peakTCh198_normPOP' num2str(filesv)],['peakCh198_normPOP' num2str(filesv)]);

%glompgthy44tr7_13_popnorm_VARSExInONLYEX', thy44tr14_19glomGv
%thy44tr14_19glomGv_popnorm_VARSExInONLYEX,thy44tr14_19glomRv_popnorm_VARSExInONLYEX



%save('TCnocck_tempSig',['halfwidth_pca' num2str(filesv)],['halfwidthmin_pca' num2str(filesv)],['min_LatY' num2str(filesv)],['min_LatX' num2str(filesv)],['name' num2str(filesv)],['SVch1_sv' num2str(filesv)],['Ch1Odor2' num2str(filesv)],['odoresp' num2str(filesv)],['odorebound' num2str(filesv)],['fastsup' num2str(filesv)],['fastsuppression' num2str(filesv)],['delayedexcit' num2str(filesv)],['halfwidthmin' num2str(filesv)],['Ch1rig' num2str(filesv)],['Ch1mice' num2str(filesv)],['Ch1mice' num2str(filesv)],['roilableCh1' num2str(filesv)],['expCh1a' num2str(filesv)],['resp_value' num2str(filesv)],['Ch1Odorinh' num2str(filesv)],['Ch1Odor' num2str(filesv)],['liqinh' num2str(filesv)],['liq' num2str(filesv)],['svonlyinh' num2str(filesv)],['svonly' num2str(filesv)],['roinameinh' num2str(filesv)],['roiname' num2str(filesv)],['Ch1Expinh' num2str(filesv)],['Ch1SVinh' num2str(filesv)],['roiCh1inh' num2str(filesv)],['latencyMininh' num2str(filesv)],['peak_LatX' num2str(filesv)],['roiCh1' num2str(filesv)],['peak_LatY' num2str(filesv)],['latencyYinh' num2str(filesv)],['latencyXinh' num2str(filesv)],['latencyY' num2str(filesv)],['latencyX' num2str(filesv)],['peakTCh115Exp' num2str(filesv)],['peakTCh115Odor' num2str(filesv)],['peakTCh115SV' num2str(filesv)],['peakTCh198Exp' num2str(filesv)],['peakTCh198Odor' num2str(filesv)],['Ch1Exp' num2str(filesv)],['Ch1Odor' num2str(filesv)],['Ch1SV' num2str(filesv)],['risetime' num2str(filesv)],['minrespCh1' num2str(filesv)],[ 'minrespTCh1' num2str(filesv)],['minCh198' num2str(filesv)],[ 'minTCh198' num2str(filesv)], ['peakCh125' num2str(filesv)],['peakTCh125' num2str(filesv)],['peakTCh175' num2str(filesv)],['peakCh1' num2str(filesv)],['peakTCh1' num2str(filesv)], ['peakCh115' num2str(filesv)],['peakTCh115' num2str(filesv)],['yinflect' num2str(filesv)],['linregCh1' num2str(filesv)],['decayCh150' num2str(filesv)],['decayTCh150' num2str(filesv)],['peakCh150' num2str(filesv)],['peakTCh150' num2str(filesv)],['halfwidth' num2str(filesv)],['peakCh198' num2str(filesv)],['peakTCh198' num2str(filesv)],['risetime_normPOP' num2str(filesv)],['halfwidth_normPOP' num2str(filesv)],['peakCh115_normPOP' num2str(filesv)],['peakTCh115_normPOP' num2str(filesv)],['peakTCh198_normPOP' num2str(filesv)],['peakCh198_normPOP' num2str(filesv)]);

figure()
scatter(latencyX-1,peak_LatY)
title('Onset Latency vs Peak Magnitude')
xlabel('response onset (inhalation onset= 0ms)')
ylabel('peak magnitude(deltF/F)')
saveas(gcf, ['MTC Onset Latency vs Peak Magnitude'], 'fig')



% figure()
% scatter(latencyXinh,latencyYinh)
% title('Channel 1 Peak latency Inhibitory Responses')
% %saveas(gcf, ['FIGURES/' 'Channel 1 Peak x Peak Time' num2str(file)], 'fig')
% 
% figure()
% scatter(latencyX,latencyY)
% title('Channel 1 Peak latency Excitatory Responses')
% %saveas(gcf, ['FIGURES/' 'Channel 1 Peak x Peak Time' num2str(file)], 'fig')
% 
% 
% figure()
% scatter(peakTCh1,peakCh1)
% title('Channel 1 Peak x Peak Time')
% %saveas(gcf, ['FIGURES/' 'Channel 1 Peak x Peak Time' num2str(file)], 'fig')

% figure()
% scatter(peakTCh198,peakCh198)
% title('Channel 1 95% Peak x Peak Time')
% saveas(gcf, ['FIGURES/' 'Channel 1 95% Peak x Peak Time' num2str(file)], 'fig')

% figure()
% scatter(linregCh1(:,1),risetime)
% title('Channel 1  slope x Peak Time')
% % saveas(gcf, ['FIGURES/' 'Channel 1  slope x Peak Time' num2str(file)], 'fig')
%
% figure()
% scatter(peakTCh1,peakCh1)
% title('Channel 1 Peak x Peak Time')
% saveas(gcf, ['FIGURES/' 'Channel 1 Peak x Peak Time' num2str(file)], 'fig')
%
% figure()
% scatter(minTCh198,minCh198)
% title('Channel 1 MIN respon x Min Time')
% saveas(gcf, ['FIGURES/' 'Channel 1 Peak x Peak Time' num2str(file)], 'fig')



%     figure()
%     scatter(peakTCh198,SVch1)
%     title('CH1: Saturated Vapor vs. Peak Amp Onset')
%     saveas(gcf, ['FIGURES/' 'CH1: Saturated Vapor vs. Peak Amp Onset' num2str(file)], 'fig')

%     figure()
%     scatter(peakCh198,SVch1)
%     title('CH1: Peak Amp vs. Saturated Vapor ')
%     saveas(gcf, ['FIGURES/' 'CH1: Peak Amp vs. Saturated Vapor' num2str(file)], 'fig')

%     figure()
%     scatter(linregCh1(:,1),SVch1)
%     title('CH1: slope vs. Saturated Vapor')
%     saveas(gcf, ['FIGURES/' 'CH1: slope vs. Saturated Vapor' num2str(file)], 'fig')

% %%%%% fit Ch2 data to a double gaussian%%%
% RiseslopesCh2=[];
% halfRiseslopesCh2=[];
% peakCh2=[];
% peakTCh2=[];
% peakCh275=[];
% peakTCh275=[];
% peakCh225=[];
% peakTCh225=[];
% halfpeakCh2=[];
% halfpeakTch2=[];
% PolyfitYCh2=[];
% halfSlopeCh2=[];
% linregCh2=[];
% peakCh298=[];
% peakTCh298=[];
% [num_of_roisCh2 col] = size(roisYCh2);%%% calculate the num of ROI per Ch2
% for i=1:num_of_roisCh2
%     x=time; roisX=time;
%     y=roisYCh2(i,:);
%     %% Calculate Derivate of MaxSlope
%     f1 = diff(y);
%     [Maxslope_onset Max_slope_dfind]=max(f1);
%     xinflect=x(Max_slope_dfind);
%     yinflect=max(f1);
%
%     %Calculate peak and peak time
%     [peak, peakIn]=max(y);
%     peakTime=x(peakIn);
%     peakCh2=[peakCh2;peak];
%     peakTCh2=[peakTCh2;peakTime];
%
%     %Calculate 98%peak and 98%peakt time
%     Per98peakInd=find(y>.98*peak,1);
%     Per98peakXt=x(Per98peakInd);
%     Per98peak=y(Per98peakInd);
%     peakCh298=[peakCh298;Per98peak];
%     peakTCh298=[peakTCh298;Per98peakXt];
%
%     %Calculate 75%peak and 75%peakt time
%     Per75peakInd=find(y>.75*peak,1);
%     Per75peakXt=x(Per75peakInd);
%     Per75peak=y(Per75peakInd);
%     peakCh275=[peakCh275;Per75peak];
%     peakTCh275=[peakTCh275;Per75peakXt];
%
%     %Calculate 25%peak and 25%peakt time
%     Per25peakInd=find(y>.25*peak,1);
%     Per25peakXt=x(Per25peakInd);
%     Per25peak=y(Per25peakInd);
%     peakCh225=[peakCh225;Per25peak];
%     peakTCh225=[peakTCh225;Per25peakXt];
%
%     %% linear regression of 25% to 95% to peak amp
%
%     xreg=x(Per25peakInd:Per75peakInd);
%     yreg=y(Per25peakInd:Per75peakInd);
%     %     b1=x(Per25peakInd:Per75peakInd)/y(Per25peakInd:Per75peakInd);
%     %     yCalc1=b1*xreg;
%     linregpoly=polyfit(xreg,yreg,1); %% first value is the slope, second value is the y-intercept
%
%
%     x2=linspace(min(xreg),max(xreg));
%     y2 = polyval(linregpoly,x2);
%     linregCh2=[linregCh2;linregpoly];%% first value is the slope, second value is the y-intercept
%

%     figure()
%     scatter(xreg,yreg)
%     hold on
%     plot(x2,y2,'r')
%
%     %% POLYFIT Func %% not great
%     p=polyfit(time,roisYCh2(i,:),40);
%     x1 = linspace(0,max(time));
%     y1 = polyval(p,x1);
%
%
%     maxP=max(y1);
%     maxPhalf=.5*(max(y1));
%
%     %Calculate 50%peak and 50%peakt time of fitted data
%     [halfP, halfpeakfitInd]=find(y1>maxPhalf,1);
%     halfpeakTime=x1(halfpeakfitInd);
%     halfpeak=y1(halfpeakfitInd);
%
%     halfpeakCh2=[halfpeakCh2;halfpeak];
%     halfpeakTch2=[halfpeakTch2;halfpeakTime];
%
%     q = polyder(p);
%     dy=polyval(q,halfpeakTch1(i));
%     xpt=halfpeakTch2(i);
%     ypt=halfpeakCh2(i);
%
%
%
%     %ypt-b=dy*xpt
%     %ypt=-1*(dy*xpt-ypt)
%     intp=-1*(dy*xpt-ypt);%-7.95%ypt/(dy*xpt);% solve eq for line
%
%
% %     m=dy;b=intp;
% %     figure()
% %     fplot(@(x)m*x+b, [0 4]);
% %
% %     hold on
% %
% %     plot(time,roisYCh1(1,:),'o')
% %     hold on
% %     plot(x1,y1,'r')
% %     hold on
% %     scatter(xpt, ypt, 'ro')
%
%     PolyfitYCh2=[PolyfitYCh2;y1];
%     halfSlopeCh2=[halfSlopeCh2;dy];
%
% end
%
%
% figure()
% scatter(peakTCh2,peakCh2)
% title('Channel 2 PG Input Peak x Peak Time')
%
% figure()
% scatter(peakTCh298,peakCh298)
% title('Channel 2 PG Input 98% Peak x Peak Time')
%
%
% figure()
% scatter(linregCh2(:,1),peakCh298)
% title('Channel 2 PG Input slope x Peak Time')
%
% figure()
% scatter(peakTCh298,SVch2)
% title('Channel 2 PG: Saturated Vapor vs. Peak Amp Onset')
%
% figure()
% scatter(peakCh298,SVch2)
% title('Channel 2 Tufted: Peak Amp vs. Saturated Vapor ')
%
% figure()
% scatter(linregCh2(:,1),SVch2)
% title('Channel 2 PG: slope vs. Saturated Vapor')
% %



% figure()
% scatter(peakCh295,SVch2)
% title('PG: Saturated Vapor vs. Peak Amp')
%
% peak_SV_Ch2=[peakCh295,SVch2];
% peak_SV_Ch1=[peakCh175,SVch1];
%
%
% high10Ch2=peak_SV_Ch2(:,2)==10;
% Ch2SV10=peak_SV_Ch2(high10Ch2,:);
% Ave10Ch2=mean(Ch2SV10);
% std10Ch2=std(Ch2SV10);
%
% high10Ch1=peak_SV_Ch1(:,2)==10;
% Ch1SV10=peak_SV_Ch1(high10Ch1,:);
% Ave10Ch1=mean(Ch1SV10);
% std10Ch1=std(Ch1SV10);
%
% high5Ch2=peak_SV_Ch2(:,2)==5;
% Ch2SV5=peak_SV_Ch2(high5Ch2,:);
% Ave5Ch2=mean(Ch2SV5);
% std5Ch2=std(Ch2SV5);
%
% high5Ch1=peak_SV_Ch1(:,2)==5;
% Ch1SV5=peak_SV_Ch1(high5Ch1,:);
% Ave5Ch1=mean(Ch1SV5);
% std5Ch1=std(Ch1SV5);
%
% high1Ch2=peak_SV_Ch2(:,2)==1;
% Ch2SV1=peak_SV_Ch2(high1Ch2,:);
% Ave1Ch2=mean(Ch2SV1);
% std1Ch2=std(Ch2SV1);
%
% high1Ch1=peak_SV_Ch1(:,2)==1;
% Ch1SV1=peak_SV_Ch1(high1Ch1,:);
% Ave1Ch1=mean(Ch1SV1);
% std1Ch1=std(Ch1SV1);
%
% high03Ch2=peak_SV_Ch2(:,2)==.3;
% Ch2SV03=peak_SV_Ch2(high03Ch2,:);
% Ave03Ch2=mean(Ch2SV03);
% std03Ch2=std(Ch2SV03);
%
% high03Ch1=peak_SV_Ch1(:,2)==.3;
% Ch1SV03=peak_SV_Ch1(high03Ch1,:);
% Ave03Ch1=mean(Ch1SV03);
% std03Ch1=std(Ch1SV03);
%
% meanRespCh2=[Ave10Ch2;Ave5Ch2;Ave1Ch2;Ave03Ch2];
% SdevRespCh2=[std10Ch2;std5Ch2;std1Ch2;std03Ch2];
%
%
% meanRespCh1=[Ave10Ch1;Ave5Ch1;Ave1Ch1;Ave03Ch1];
% SdevRespCh1=[std10Ch1;std5Ch1;std1Ch1;std03Ch1];
%
% figure()
% hold on
% bar(1:4,meanRespCh2(:,1))
% errorbar(1:4,meanRespCh2(:,1),SdevRespCh2(:,1),'.')
% title('mean PG population response magnitude across odor intensity')
%
% figure()
% hold on
% bar(1:4,meanRespCh1(:,1))
% errorbar(1:4,meanRespCh1(:,1),SdevRespCh1(:,1),'.')
% title('mean tufted population response magnitude across odor intensity')
%
%
% %% STORE FITTED MAX, MIN, and their indicies %%%
% %%%% GENERATE FIGURES of Channel 2 %%%%
% roisY=roisYCh2;
% figure()
% imagesc(roisY);
% set(gca,'YTickLabel',{num2str(SVch1)})
% colorbar
% title('RAW traces: unsorted, not normalized')
%
% %oder by strongest response
% a=roisYCh2;
% [ROIpeak, ROIpeakInd]=max(roisY,[],2); %vector of max value for each ROI
% [dummy,index]=sort(ROIpeak);
% PeakAmp=a(index,:);
%
% figure()
% imagesc(PeakAmp)
% set(gca,'Ydir','Normal')
% colorbar
% title('order by strongest responder')
%
%
% %oder by strongest response
% a=roisYCh1;
% [ROIpeak, ROIpeakInd]=max(roisY,[],2); %vector of max value for each ROI
% [dummy,index]=sort(ROIpeakInd);
% PeakAmp=a(index,:);
%
% figure()
% imagesc(PeakAmp)
% set(gca,'Ydir','Normal')
% colorbar
% title('order by peak onset')
%
%
%
% b=PeakAmp;
% m = max(abs(b), [], 2);
% %m=mHIGH
% b2 = bsxfun(@rdivide, b, m);%normalize each row(ROI) by its own maximum value
% PeakAmp=b2;
% figure()
% imagesc(b2)
% colormap(bluewhitered)
% colorbar
% title('normalizes to and sorted by max')
%
%
% %sort by % max onset
% thrtenp=[];
% tenperc_ROIpeak_top=.5*ROIpeak;
% for i=1:length(tenperc_ROIpeak_top)
%     [thrX thrY]=find(b2(i,:)>=tenperc_ROIpeak_top(i,:),1);
%     thrtenp=[thrtenp;thrY];
% end
%
% figure()
% title('CD plot of Tufted(g) and Pg(r) glomerular response')
% cdfplot(peakTCh175)
% hold on; h = cdfplot(peakTCh175);
% set(h,'color','g')
% hold on
% cdfplot(peakTCh295)
% hold on; h = cdfplot(peakTCh175);
% set(h,'color','r')
%
%
%
%
% %
% %

% %
% %
% % %%% Isolate Excitatory Responses %%%
% % thrtenp=[];
% % tenperc_ROIpeak_top=.25*ROIpeak;
% % for i=1:length(tenperc_ROIpeak_top)
% %     [thrX thrY]=find(roisEX(i,:)>=tenperc_ROIpeak_top(i,:),1);
% %     thrtenp=[thrtenp;thrY];
% % end
% % [dummy,index] = sort(thrtenp);%%%%sort by threshold
% % roisEX=roisEX(index,:);
% %
% % m = max(abs(roisEX), [], 2);
% % roisEX = bsxfun(@rdivide, roisEX, m);%normalize each row(ROI) by its own maximum value
% %
% %
% % %%% Isolate Inhibitory Responses %%%
% % thrtenp=[];
% % tenperc_ROIpeak_top=.25*ROIpeak;
% % for i=1:length(tenperc_ROIpeak_top)
% %     if sum(roisIn)<0
% %     elseif sum(roisIn)>0
% %         [thrX thrY]=find(roisIn(i,:)>=tenperc_ROIpeak_top(i,:),1);
% %         thrtenp=[thrtenp;thrY];
% %     end
% % end
% % [dummy,index] = sort(thrtenp);%%%%sort by threshold
% % roisIn=roisIn(index,:);
% %
% % m = max(abs(roisIn), [], 2);
% % %m=mHIGH
% % roisIn = bsxfun(@rdivide, roisIn, m);%normalize each row(ROI) by its own maximum value
% %
% % combine=[roisEX;roisIn];
% % figure(9)
% % imagesc(combine)
% %
% % colorbar
% % title('separate excitatory from inhibitory')
% %
% %
% % [dummy,index] = sort(SV);%%%%sort by threshold
% % PA=roisYorig(index,:);
% %
% % %PA=roisYorig
% % m = max(abs(PA), [], 2);
% % %m=mHIGH
% % PA = bsxfun(@rdivide, PA, m);%norm
% % figure(10)
% %
% % imagesc(PA)
% % colormap(bluewhitered)
% % colorbar
% % title('sort by SV')
% % NumTicks = length(SV);
% % L = get(gca,'YLim');
% % set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% %
% % set(gca,'YTickLabel',{num2str(dummy)})
% %
% %
% % % %
% % % % %%% SPLIT EXCITATORY MATRIX IN N Groups%%%
% % % % Aex = roisEX;
% % % % split=round(num_of_rois/6);
% % % % first=1;
% % % % second=split;
% % % % third=2*split;
% % % % fourth=3*split;
% % % % fifth=4*split;
% % % % sixth=5*split;
% % % %
% % % % Q1 = Aex((1:split),:);
% % % % Q2 = Aex((split:2*split),:);
% % % % Q3 = Aex((2*split:3*split),:);
% % % % Q4 = Aex((3*split:4*split),:);
% % % % Q5 = Aex((4*split:5*split),:);
% % % % Q6 = Aex((5*split:end),:);
% % % %
% % % % meanQ1=mean(Q1);
% % % % meanQ2=mean(Q2);
% % % % meanQ3=mean(Q3);
% % % % meanQ4=mean(Q4);
% % % % meanQ5=mean(Q5);
% % % % meanQ6=mean(Q6);
% % % %
% % % % figure(10)
% % % % plot(meanQ1)
% % % % hold on
% % % % plot(meanQ2)
% % % % hold on
% % % % plot(meanQ3)
% % % % hold on
% % % % plot(meanQ4)
% % % % hold on
% % % % plot(meanQ5)
% % % % hold on
% % % % plot(meanQ6)
% % % %
% % % %
% % % % Quarters=vertcat(meanQ1,meanQ2, meanQ3,meanQ4,meanQ5, meanQ6);
% % % %
% % % %
% % % %
% % % % %%% SPLIT EXCITATORY MATRIX IN N Groups%%%
% % % % Ain = roisIn;
% % % % split=round(num_of_rois/6);
% % % % first=1;
% % % % second=split;
% % % % third=2*split;
% % % % fourth=3*split;
% % % % fifth=4*split;
% % % % sixth=5*split;
% % % %
% % % % Q1 = Ain((1:split),:);
% % % % Q2 = Ain((split:2*split),:);
% % % % Q3 = Ain((2*split:3*split),:);
% % % % Q4 = Ain((3*split:4*split),:);
% % % % Q5 = Ain((4*split:5*split),:);
% % % % Q6 = Ain((5*split:end),:);
% % % %
% % % % meanQ1=mean(Q1);
% % % % meanQ2=mean(Q2);
% % % % meanQ3=mean(Q3);
% % % % meanQ4=mean(Q4);
% % % % meanQ5=mean(Q5);
% % % % meanQ6=mean(Q6);
% % % %
% % % % figure(11)
% % % % plot(meanQ1)
% % % % hold on
% % % % plot(meanQ2)
% % % % hold on
% % % % plot(meanQ3)
% % % % hold on
% % % % plot(meanQ4)
% % % % hold on
% % % % plot(meanQ5)
% % % % hold on
% % % % plot(meanQ6)
% % % %
% % % %
% % % % Quarters=vertcat(meanQ1,meanQ2, meanQ3,meanQ4,meanQ5, meanQ6);


