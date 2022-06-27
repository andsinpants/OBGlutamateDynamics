respmag=[];
preodor=1:550; %odor on time = 600, on for 8 sec
odoron=600:1800;
numtrials = 195;
gscaledtraces=iglusnfrtraces; %duplicate matrix for scaling
rscaledtraces=rgecotraces;
for i = 1:numtrials;
    gtrace=iglusnfrtraces(i,:);
    maxval=max(gtrace(odoron));
    minval=min(gtrace(odoron));
    if abs(minval)>maxval; maxval=abs(minval); end;
    gscaledtraces(i,:)=gtrace./maxval;
    respmag(i)=sum(gscaledtraces(i,odoron));
    rtrace=rgecotraces(i,:);
    maxval=max(rtrace(odoron));
    minval=min(rtrace(odoron));
    if abs(minval)>maxval; maxval=abs(minval); end;
    rscaledtraces(i,:)=rtrace./maxval;    
end
[sorted, sortorder]=sort(respmag);
sortedgtraces=gscaledtraces(sortorder,:);
sortedrtraces=rscaledtraces(sortorder,:);
figure;
imagesc(sortedgtraces);
figure;
imagesc(sortedrtraces);



