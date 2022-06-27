sampleperiodms=1000/150;
pretime=50
sta1=mean(snifftrig_roi1);
sta1=sta1-min(sta1);
sta2=mean(snifftrig_roi2);
sta2=sta2-min(sta2);
sta3=mean(snifftrig_roi3);
sta3=sta3-min(sta3);
% sta4=mean(snifftrig_roi4);
% sta4=sta4-min(sta4);
% sta5=mean(snifftrig_roi5);
% sta5=sta5-min(sta5);
% sta6=mean(snifftrig_roi6);
% sta6=sta6-min(sta6);
% sta7=mean(snifftrig_roi7);
% sta7=sta7-min(sta7);
stas=vertcat(sta1, sta2, sta3);

for i = 1:3
    [peak, peaktime]=max(stas(i,pretime:end));
    [trough, mintime]=min(stas(i,:));
    sta_stats(i,1)=peak;
    sta_stats(i,2)=(mintime-pretime)*sampleperiodms;
    sta_stats(i,3)=(peaktime)*sampleperiodms;
end
    