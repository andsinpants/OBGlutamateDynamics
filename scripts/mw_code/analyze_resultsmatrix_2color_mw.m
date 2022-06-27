
% analyze thresholded results matrix
totalnum=length(results3(:,1));   % number of glom-odor pairs in analysis set
%calculate prevalanece and ocncordance of excitatory and suppressive responses
numresp(1)=length(find(results3(:,3))); %numresp=ex1, ex2, sup1, sup2
numresp(2)=length(find(results3(:,4)));
numresp(3)=length(find(results3(:,5)));
numresp(4)=length(find(results3(:,6)));
concord_ex=results3(:,3) & results3(:,4);
numbothex=length(find(concord_ex));
concord1_2=numresp(1)/numbothex; %this is fraction of glom-odor prs showing green signal that also show red signal
concord2_1=numresp(2)/numbothex; %this is fraction of glom-odor prs showing red signal that also show green signal 
concord_sup=results3(:,5) & results3(:,6);
numbothsupp=length(find(concord_sup));
concord_supp1_2=numresp(3)/numbothsupp;
concord_supp2_1=numresp(4)/numbothsupp;
exc1_sup2=results3(:,3) & results3(:,6); 
numexc1_sup2=length(find(exc1_sup2)); %this is number of g-o pairs that show suppression but not excitation.
t2t1tmaxvals=results3(find(concord_ex), 8:9); %trt1tmax vals for green and red channels.






