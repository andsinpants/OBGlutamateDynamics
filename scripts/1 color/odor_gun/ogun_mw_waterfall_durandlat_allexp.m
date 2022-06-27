%%load normsort e_l,d and normsort_s_l,d into workspace
clearvars
close all
load('allexp_resultssig_vars_1color.mat');
%% sort by odor latency to peak
% only excitatory
iodor = 301:726;
times = (0:size(normsorte_l,2)-1)/150; 
xpeaks=[];
for k = 1:size(normsorte_l,1)
   xpeaks(k) = max(normsorte_l(k,iodor(1):iodor(end)));
end

xpeaks = xpeaks';

peaktimes = [];
for k=1:size(normsorte_l,1)
    peaktimes(k)=find(normsorte_l(k,5:end)==xpeaks(k));
end

[sortedpeaktimes,peakind] = sort(peaktimes,'ascend');
normsorte_l_sort = normsorte_l(peakind,:);

figure;imagesc(normsorte_l_sort(:,5:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'Norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
hold on;
plot(peaktimes(peakind),1:length(peaktimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};
hold off;
%% plot sup glusnfr rois normalized & sorted to trough
% now suppressive
xtroughs=[];
troughtimes=[];
for i=1
    if isempty(normsorts_l)
        normsorts = [];
        break
    else
%         xpeaks=[];
        for k = 1:size(normsorts_l,1)
           xtroughs(k,:) = min(normsorts_l(k,iodor(1):iodor(end)));
        end
%         peaktimes=[];
        for k=1:size(normsorts_l,1)
            troughtimes(k)=find(normsorts_l(k,5:end)==xtroughs(k));
        end
    end
end
[sortedtroughtimes,troughind] = sort(troughtimes,'descend');
normsorts_l_sort = normsorts_l(troughind,:);

figure;imagesc(normsorts_l_sort(:,5:end));title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{blue}peak latency}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'Norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([451 451 751 751],[1 1 1 1], gray);
hold on;
plot(troughtimes(troughind),1:length(troughtimes),'color','k')
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};
hold off;
%% sort by odor duration
% only excitatory
for k = 1:size(normsorte_d,1)
   xfwhm(k) = fwhm(times,normsorte_d(k,5:end));
end

[sorteddurtimes,fwhmind] = sort(xfwhm,'ascend');
normsorte_d_sort = normsorte_d(fwhmind,5:end);

figure;imagesc(normsorte_d_sort(:,5:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs sorted by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. z-scored \DeltaF/F';
clims = [-1 1];
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};

% now suppressive
xfwhm=[];
fwmind=[];
for i=1
    if isempty(normsorts_d)
        normsorts_d = [];
    else 
        for k = 1:size(normsorts_d,1)
            xfwhm(k) = fwhm(times,normsorts_d(k,5:end));
        end
    end
[sorteddurtimes,fwhmind] = sort(xfwhm,'descend');
normsorts_d_sort = normsorts_d(fwhmind,5:end);

figure;imagesc(normsorts_d_sort(:,5:end));title(['Significant Suppressive Norm & Sorted Glomerular ROIs by {\color{red}duration}']);ylabel('Glom-Odor Pairs');xlabel('Frames');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([450 450 750 750],[1 1 1 1], gray);
hax=gca;
hax.XTick = [151,451,751,1051,1351,1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};
%figure;bar(peaktimes(peakind),1:length(peaktimes),'color','k')
end

%% sort by peak response magnitude during odor times
results_sig_responsive = [];
netresponse=[];
for i = 1:size(results_sig,1)
    sigtest=[results_sig(i,3), results_sig(i,4)];
    if sum(sigtest) ~= 0
        results_sig_responsive=vertcat(results_sig_responsive,results_sig(i,:));
        netexcit=sigtest(1);
        netsupp=sigtest(2); 
        if netexcit==0
            tmpresp=netsupp;
        else tmpresp=netexcit;
        end
     netresponse=vertcat(netresponse,tmpresp);
    end     
end

resp_sig_cat = horzcat(netresponse,results_sig_responsive);

resp_sortsig_cat = sortrows(resp_sig_cat,1);

for j=1:size(resp_sortsig_cat,1) %% in progress
    normsort_all(j,:) = normalised_diff(resp_sortsig_cat(j,6:iodor(1)+6:iodor(end)+6:end));
end

figure;imagesc(normsort_all(:,6:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{green}peak response magnitude}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'Norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

%% now find sum of responses between iodor times (iodor(1):iodor(end))
for j=1:size(resp_sortsig_cat,1) %% in progress
%     normosum_all(j,:) = normalised_diff(results_sig_responsive(j,5:end));
    if abs(results_sig_responsive(j,4)) > results_sig_responsive(j,3)
      min_range_value = results_sig_responsive(j,4);
      max_range_value = -min_range_value;
  else
      max_range_value = results_sig_responsive(j,3);
      min_range_value = -max_range_value;
  end
% norm_value = 2 .* data ./ (max_range_value - min_range_value);
normosum_all(j,:)=  2 .* results_sig_responsive(j,5:end) ./ (max_range_value - min_range_value);  
end

sigtest = [];
for i = 1:size(normosum_all,1)
    sigtest(i,1:length(iodor)+450)=normosum_all(i,iodor(1):iodor(end)+450);%% change to amount of time to find sums by and sort 600, 450, etc
    sumch(i) = sum(sigtest(i,:));
end

sumch = sumch';

sumodor_sig_cat = horzcat(results_sig_responsive(:,1:4),sumch,normosum_all);

sumodor_sortsig_cat = sortrows(sumodor_sig_cat,5);

figure;imagesc(sumodor_sortsig_cat(:,6:end));title(['Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{magenta}sum of iodor times}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'Norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351 1651];
hax.XTickLabel = {'1';'3';'5';'7';'9';'11'};

% clearvars -except normsorte_l normsorte_d normsorts_l normsorts_d normsorte_l_sort normsorte_d_sort normsorts_l_sort normsorts_d_sort