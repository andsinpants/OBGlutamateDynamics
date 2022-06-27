onsetlatarray=[];
roiarray=[];
allresultsfilt=allresults;
for i = 1:length(allresults.trial)
    onsetlatarray(i)=allresults.trial(i).onsetlat;
    roiarray(i)=allresults.trial(i).roi;
end
nr=find(onsetlatarray==0);
%roiarray(nr)=[];
% omit=find(roiarray==3);
% nr=[nr, omit];
% omit=find(roiarray==9);
% nr=[nr, omit];
allresultsfilt.trial(nr)=[];