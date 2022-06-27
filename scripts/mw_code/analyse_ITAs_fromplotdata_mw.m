%matt making merged STa analysis script
%just work on current myplotdaa set that is in memory
t10=0;
t90=0;
allresults=[];
decayTCh150=[];
decayCh150=[];
count=0;
todoron=150;  %hardcode in time of stim on in frames (150 Hz pseudoframe rate)
duration=250;   %this used to find max and look for onset latency. Not for FWHM.
todoroff=todoron+duration;  %define measurement/analysis time
sniffdelay=0;
preodorduration=125;
posThreshold = 4; %%change for threshold
frameinterval=1/150;
numfiles=length(myplotdata.file);
numrois=length(myplotdata.file(1).roi);
for j=1:numrois
    for i = 1:numfiles
        count=count+1;
        allresults.trial(count).filename=myplotdata.file(i).name;
        allresults.trial(count).roi=myplotdata.file(i).roi(j).number;
        allresults.trial(count).odorname=myplotdata.file(i).type;  %change this later to be odorname
        trace=myplotdata.file(i).roi(j).odor.avgtrial.series;  %pull out trace in question
        %this part subtracts trend line and then calcualtes stdeviation for
        %z-scoring
        prestimperiod=trace(todoron-preodorduration:todoron);
        prestimperiodx=find(prestimperiod);
        c = polyfit(prestimperiodx,prestimperiod,1);
        y_est = polyval(c,prestimperiodx);
        linesub_prestimperiod=prestimperiod-y_est;
        prenoise=std(linesub_prestimperiod);
        lastptfit=y_est(end);
        trace=trace-lastptfit;
        %poststimY=trace(todoron:todoroff);        
       % prenoise=std(trace(todoron-preodorduration:todoron)); %noise of baseline
        %ztrace=trace-mean(trace(todoron-preodorduration:todoron));
        ztrace=trace./prenoise;  %z-score trace based on stdev of pre-time after detrending
        maxz=max(ztrace);
        %allresults.trial(count).peakamp=max(trace(todoron:todoroff)); %peak amp = max of trace
        allresults.trial(count).peakampz=maxz; 
        allresults.trial(count).peakamp=prctile(trace(todoron:todoroff),90); %find 90th percentile
        % note this is from SHaina's code to find onset latency
        poststimY=ztrace(todoron:todoroff);    
        aboveThreshold = (poststimY > posThreshold)';
        thresholdChange = [aboveThreshold(1) diff(aboveThreshold)];
        thresholdChange(thresholdChange==-1) = 0;
        spanLocs = cumsum(thresholdChange);
        spanLocs(~aboveThreshold) = 0;
        aboveThresholdIndex = find(aboveThreshold==1);
        notConsecutiveIndex = [true diff(aboveThresholdIndex) ~= 1];
        sumIndex = cumsum(notConsecutiveIndex);
        spanLength = histc(sumIndex, 1:sumIndex(end));
        goodSpans = find(spanLength>=30);%number of consecutive points required is hard-coded here...
        allInSpans = find(ismember(spanLocs, goodSpans));
        latencysd_in=min(allInSpans)*frameinterval;
        if latencysd_in
            allresults.trial(count).onsetlat=(latencysd_in-sniffdelay)*1000;  % note will be empty if never gets above threshold for long enough
        else allresults.trial(count).onsetlat=0;
        end
        if latencysd_in   %only do all this if signal gets above threshold for long enough.
        %now Shaina's code to do smoothing and find peak times and FWHM
            
             ysmooth2 = smoothdata(ztrace(1:todoroff),'gaussian',40);   %what is this width? Assume that width = 6sigma, so, length of 40 = 8.3. which = sigma of 45 msec at .067 frame interval
             [peak, peakIn]=max(ysmooth2);
             [R,t10,t90,LL,UL]=risetime(ysmooth2(todoron:todoroff),150,'StateLevels',[0 peak]);  %adding in AKM's way of measureing t10 and t90
             peakTime=(peakIn-todoron)*frameinterval; %time of smoothed peak.
             Per150peakInd=find(ysmooth2>.5*peak);  
             Per150peakXt=Per150peakInd;%x(Per150peakInd);
             Per150peak=ysmooth2(Per150peakInd);
             [Per150peakInd2]=find(Per150peakXt>1,1);
             Per150peak=Per150peak(Per150peakInd2);
             Per150peakXt=Per150peakXt(Per150peakInd2);
             %Calculate 50%decay and 50%decayt time
             Per50decayInd=find(ysmooth2>(.5*peak));
             Per50decayIndq=Per50decayInd(end);
             Per50decayXt=Per50decayIndq; %x(Per50decayIndq);
             Per50decay=ysmooth2(Per50decayIndq);
             decayCh150=[decayCh150;Per50decay];  %not using these for now
             decayTCh150=[decayTCh150;Per50decayXt];
             halfwidth_pca1=(Per50decayXt-Per150peakXt)*frameinterval;
        else
            peakTime=0;
            halfwidth_pca1=0;
            t10=0;
            t90=0;
            ysmooth2 = smoothdata(ztrace(1:todoroff),'gaussian',40); 
        end
        allresults.trial(count).t10=(t10-sniffdelay)*1000;
        allresults.trial(count).peaktime=(peakTime-sniffdelay)*1000;
        allresults.trial(count).t90=(t90-sniffdelay)*1000;        
        allresults.trial(count).fwhm=halfwidth_pca1;
        allresults.trial(count).ztrace=ztrace;
        allresults.trial(count).smtrace=ysmooth2;
        
    end
end
%allresultsfilt=allresults;
% for i = 1:length(allresultsfilt.trial)
%     onsetlatarray(i)=allresultsfilt.trial(i).onsetlat;
%     t10array(i)=allresultsfilt.trial(i).t10;
% end
% nr=find(onsetlatarray==0);
% allresultsfilt.trial(nr)=[];