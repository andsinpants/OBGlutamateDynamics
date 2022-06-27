% close all
clear all
load('allexp_resultssig_vars_2color.mat') %change for either including 
% load('allexp_resultssig_vars_2color_corrcoef')
%% keep all rows with significant values in either matrix
iodor = 301:825;

%% find sum of z-scores and sort by responses
results_sig_responsive1 = [];
results_sig_responsive2 = [];
netresponse=[];
for i = 1:size(results_sig_all,1)
    sigtest=[results_sig_all(i,3), results_sig_all(i,4), results_sig2_all(i,3), results_sig2_all(i, 4)];
    if sum(sigtest) ~= 0
        results_sig_responsive1=vertcat(results_sig_responsive1,results_sig_all(i,:));
        results_sig_responsive2=vertcat(results_sig_responsive2, results_sig2_all(i,:));
        netexcit=sigtest(1) + sigtest(3);
        netsupp=sigtest(2) + sigtest(4);
        if netexcit==0
            tmpresp=netsupp;
        else tmpresp=netexcit;
        end
     netresponse=vertcat(netresponse, tmpresp);
    end     
end

resp_sig_cat1 = horzcat(netresponse,results_sig_responsive1);
resp_sig_cat2 = horzcat(netresponse,results_sig_responsive2);

resp_sortsig_cat1 = sortrows(resp_sig_cat1,1);
resp_sortsig_cat2 = sortrows(resp_sig_cat2,1);

for j=1:size(resp_sortsig_cat1,1) %% in progress
    normsortch1_all(j,:) = normalised_diff(resp_sortsig_cat1(j,7:iodor(1)+7:iodor(end)+7:end));
end

figure;imagesc(normsortch1_all(:,7:end));title(['Ch1 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{green}response magnitude}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

for j=1:size(resp_sortsig_cat2,1) %% in progress
    normsortch2_all(j,:) = normalised_diff(resp_sortsig_cat2(j,7:iodor(1)+7:iodor(end)+7:end));
end


figure;imagesc(normsortch2_all(:,7:end));title(['Ch2 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{green}response magnitude}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Limits = [-1 1];
c.Label.String = 'norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

% for i = 1:size(results_sig2_all,1)
%     tf = results_sig_all(:,3)>0 & results_sig_all(:,4)==0;
%     gf = results_sig_all(:,4)<0 & results_sig_all(:,3)==0;
%     tgf = results_sig_all(:,3)>0 & results_sig_all(:,4)<0;
%        if tf(i)==1
%            exc_all = results_sig_all(tf,:);
%        elseif gf(i)==1
%            sup_all = results_sig_all(gf,:);
%        else tgf(i)==1;
%            biph_all = results_sig_all(tgf,:);
%end
% end

% for i = 1:size(results_sig2_all,1)
%     tf2 = results_sig2_all(:,3)>0 & results_sig2_all(:,4)==0;
%     gf2 = results_sig2_all(:,4)<0 & results_sig2_all(:,3)==0;
%     tgf2 = results_sig2_all(:,3)>0 & results_sig2_all(:,4)<0;
%        if tf2(i)==1
%            exc_all2 = results_sig2_all(tf2,:);
%        elseif gf2(i)==1
%            sup_all2 = results_sig2_all(gf2,:);
%        else tgf2(i)==1;
%            biph_all2 = results_sig2_all(tgf2,:);
%        end
% end

times = (0:size(results_sig_all,2)-1)/150; 
%% now find sum of responses between iodor times (iodor(1):iodor(end))
% normosumch1_all=[];
% normosumch2_all=[];
% % for j=1:size(resp_sortsig_cat1,1) %% in progress
% %     %normosumch1_all(j,:) = normalised_diff(results_sig_responsive1(j,5:iodor(1)+5:iodor(end)+5:end));
% %   if abs(resp_sortsig_cat1(j,5)) > resp_sortsig_cat1(j,4)
% %       min_range_value = resp_sortsig_cat1(j,5);
% %       max_range_value = -min_range_value;
% %       normosumch1_all(j,:)=  2 .* resp_sortsig_cat1(j,6:end) ./ (max_range_value - min_range_value);
% %   elseif (resp_sortsig_cat1(j,5)) < resp_sortsig_cat1(j,4)
% %       max_range_value = resp_sortsig_cat1(j,4);
% %       min_range_value = -max_range_value;
% %       normosumch1_all(j,:)=  2 .* resp_sortsig_cat1(j,6:end) ./ (max_range_value - min_range_value);
% %   else
% %       normosumch1_all(j,:) = normalised_diff(resp_sortsig_cat1(j,6:end));
% %   end
% % end
% % 
% % 
% % clearvars min_range_value max_range_value j
% % 
% % for j=1:size(resp_sortsig_cat2,1) %% in progress
% %     %normosumch1_all(j,:) = normalised_diff(results_sig_responsive1(j,5:iodor(1)+5:iodor(end)+5:end));
% %   if abs(resp_sortsig_cat2(j,5)) > resp_sortsig_cat2(j,4)
% %       min_range_value = resp_sortsig_cat2(j,5);
% %       max_range_value = -min_range_value;
% %       normosumch2_all(j,:)=  2 .* resp_sortsig_cat2(j,6:end) ./ (max_range_value - min_range_value);
% %   elseif (resp_sortsig_cat2(j,5)) < resp_sortsig_cat2(j,4)
% %       max_range_value = resp_sortsig_cat2(j,4);
% %       min_range_value = -max_range_value;
% %       normosumch2_all(j,:)=  2 .* resp_sortsig_cat2(j,6:end) ./ (max_range_value - min_range_value);
% %   else
% %       normosumch2_all(j,:) = normalised_diff(resp_sortsig_cat2(j,6:end));
% %   end
% % end

% normosumch1_all = results_sig_responsive1(:,5:end);
% normosumch2_all = results_sig_responsive2(:,5:end);

normosumch1_all = normsortch1_all;
normosumch2_all = normsortch2_all;

for i = 1:size(normosumch1_all,1)
    sigtestch1(i,1:length(iodor)+450)=normosumch1_all(i,iodor(1):iodor(end)+450);
    sigtestch2(i,1:length(iodor)+450)=normosumch2_all(i,iodor(1):iodor(end)+450);
    sumch1(i) = sum(sigtestch1(i,:));
    sumch2(i) = sum(sigtestch2(i,:));
    end

sumch1 = sumch1';
sumch2 = sumch2';
sumodor_sig_cat1 = horzcat(resp_sortsig_cat1(:,4),resp_sortsig_cat1(:,5),sumch1,resp_sortsig_cat1(:,6),normosumch1_all);
sumodor_sig_cat2 = horzcat(resp_sortsig_cat2(:,4),resp_sortsig_cat2(:,5),sumch1,resp_sortsig_cat2(:,6),normosumch2_all);

sumodor_sortsig_cat1 = sortrows(sumodor_sig_cat1,3);
sumodor_sortsig_cat2 = sortrows(sumodor_sig_cat2,3);

% % figure;imagesc(sumodor_sortsig_cat1(:,4:end));title(['Ch1 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{magenta}sum of iodor times}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Limits = [-1 1];
% % c.Label.String = 'norm. z-scored \DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([301 301 825 825],[1 1 1 1], gray);
% % hold on;
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};
% % 
% % figure;imagesc(sumodor_sortsig_cat2(:,4:end));title(['Ch2 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{magenta}sum of iodor times}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar_auto);
% % c = colorbar;
% % c.Limits = [-1 1];
% % c.Label.String = 'norm. z-scored \DeltaF/F';
% % gray = [.9 .9 .9]; %odor onset and offset signal overlay 
% % patch([301 301 825 825],[1 1 1 1], gray);
% % hold on;
% % hax=gca;
% % hax.XTick = [151,451,751,1051,1351];
% % hax.XTickLabel = {'1';'3';'5';'7';'9'};

% % % for j=1:size(resp_sortsig_cat1,1)
% % %     if sumodor_sortsig_cat1(j,1)==0 && sumodor_sortsig_cat1(j,2)==0
% % %         sumodor_sortsig_cat1(j,5:end)= 0;
% % %     end
% % % end
% % % 
% % % 
% % % for j=1:size(resp_sortsig_cat2,1)
% % %     if sumodor_sortsig_cat2(j,1) ==0 && sumodor_sortsig_cat2(j,2)==0
% % %         sumodor_sortsig_cat2(j,5:end)= 0;
% % %     end
% % % end
% [gcmap]=buildcmap('kwg')
% [mcmap]=buildcmap('kwm')

ch1corrcoef = sumodor_sortsig_cat1(:,4);
ch2corrcoef = sumodor_sortsig_cat2(:,4);
index = 1:541;
ch1corrcoef_flip = flip(ch1corrcoef);
findex = flip(index)';
figure;

sz = 6;
for i=1:1:81
scatter(ch1corrcoef(i),findex(i),sz,'bo','filled'); hold on;
end
for i=82:541
scatter(ch1corrcoef(i),findex(i),sz,'ro','filled'); hold on;
end

figure;imagesc(sumodor_sortsig_cat1(:,4:end));title(['Ch1 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{magenta}sum of iodor times}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar);
c = colorbar;
% % c.Limits = [-1 1];
cmin = -1;
cmax = 1;
caxis([cmin cmax]);
%c.Limits = [-7 7];
c.Label.String = 'norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};

figure;imagesc(sumodor_sortsig_cat2(:,4:end));title(['Ch2 Significant Excitatory Norm & Sorted Glomerular ROIs by {\color{magenta}sum of iodor times}']);ylabel('Glom-Odor Pairs');xlabel('Time (s)');colormap(gca,nawhimar);
c = colorbar;
% % c.Limits = [-1 1];
cmin = -1;
cmax = 1;
caxis([cmin cmax]);
c.Label.String = 'norm. z-scored \DeltaF/F';
gray = [.9 .9 .9]; %odor onset and offset signal overlay 
patch([301 301 825 825],[1 1 1 1], gray);
hold on;
hax=gca;
hax.XTick = [151,451,751,1051,1351];
hax.XTickLabel = {'1';'3';'5';'7';'9'};
% clearvars -except exc_all exc_all2 sup_all sup_all2 biph_all biph_all2
%% pull out Ch2 responses that show signficance in Ch1
% % for i = 1:size(results_sig_all,1)
% %     tf = results_sig_all(:,3)>0;
% %      if tf(i)==1
% %          ch1_exsig = results_sig_all(tf,:);
% %          ch2_expartner = results_sig2_all(tf,:);
% %      end
% % end
% % 
% % for i = 1:size(results_sig2_all,1)
% %     tf2 = results_sig2_all(:,3)>0;
% %      if tf2(i)==1
% %          ch1_expartner = results_sig_all(tf2,:);
% %          ch2_exsig = results_sig2_all(tf2,:);
% %      end
% % end
% % 
% % for i = 1:size(results_sig_all,1)
% %     gf = results_sig_all(:,4)<0;
% %      if gf(i)==1
% %          ch1_supsig = results_sig_all(gf,:);
% %          ch2_suppartner = results_sig2_all(gf,:);
% %      end
% % end
% % 
% % for i = 1:size(results_sig2_all,1)
% %     gf2 = results_sig2_all(:,4)<0;
% %      if gf2(i)==1
% %          ch1_suppartner = results_sig_all(gf2,:);
% %          ch2_supsig = results_sig2_all(gf2,:);
% %      end
% % end



