tbt132g_iglusnfrtraces=[];
tbt132g_rgecotraces=[];
count=0;
for i = 1:length(results3(:,1));
   trialnum=results3(i,2);
   if find(results3(i,3:6))
       roi = results3(i,1);
       trial=results3(i,2);  
       count=count+1;
       iglusnfrtrace=iglusnfrplotdata.file(trial).roi(roi).odor.avgtrial.series;
       rgecotrace=rgecoplotdata.file(trial).roi(roi).odor.avgtrial.series;
       tbt132g_iglusnfrtraces(count,:)=iglusnfrtrace;
       tbt132g_rgecotraces(count,:)=rgecotrace;       
   end   
end