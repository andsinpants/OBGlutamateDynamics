clearvars
close all
load('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SFiGluSnFR\tbt_SFiGluSnFR_ogun analysis\allexp_resultssig_vars_1color_forsparsenessanalysis')

sort_tbt51_fov1_sig = sortrows(results_sig_tbt51_fov1,2);
sort_tbt51_fov2_sig = sortrows(results_sig_tbt51_fov2,2);
sort_tbt57_fov1_sig = sortrows(results_sig_tbt57_fov1,2);
sort_tbt57_fov2_sig = sortrows(results_sig_tbt57_fov2,2);
sort_tbt61_lobfov1_sig = sortrows(results_sig_tbt61_lobfov1,2);
sort_tbt61_lobfov2_sig = sortrows(results_sig_tbt61_lobfov2,2);
sort_tbt61_rob_sig = sortrows(results_sig_tbt61_rob,2);

%% tbt51_fov1
tbt51_fov1roinum = max(sort_tbt51_fov1_sig(:,1));

tbt51_fov1mat = zeros(tbt51_fov1roinum,136);

for i = 1:size(sort_tbt51_fov1_sig,1)
    roi=sort_tbt51_fov1_sig(i,1);
    odor=sort_tbt51_fov1_sig(i,2);
    tbt51_fov1mat(roi,odor)=sort_tbt51_fov1_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt51_fov1_real=tbt51_fov1mat(:, rep_odor)';
tmpdata=tbt51_fov1_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt51fov1(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt51fov1); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt51fov1(:,deletecols) = [];

meantbt51_fov1_LS = mean(sparsenessval_tbt51fov1);
stdtbt51_fov1_LS = std(sparsenessval_tbt51fov1);
 
% clearvars -except tbt51_fov1_real sparsenessval_tbt51fov1 meantbt51_fov1_LS  stdtbt51_fov1_LS
 
%% tbt51_fov2
tbt51_fov2roinum = max(sort_tbt51_fov2_sig(:,1));

tbt51_fov2mat = zeros(tbt51_fov2roinum,136);

for i = 1:size(sort_tbt51_fov2_sig,1)
    roi=sort_tbt51_fov2_sig(i,1);
    odor=sort_tbt51_fov2_sig(i,2);
    tbt51_fov2mat(roi,odor)=sort_tbt51_fov2_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt51_fov2_real=tbt51_fov2mat(:, rep_odor);
tmpdata=tbt51_fov2_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt51fov2(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt51fov2); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt51fov2(:,deletecols) = [];

meantbt51_fov2_LS = mean(sparsenessval_tbt51fov2);
stdtbt51_fov2_LS = std(sparsenessval_tbt51fov2);
 
% clearvars -except tbt51_fov2_real sparsenessval_tbt51fov2 meantbt51_fov2_LS  stdtbt51_fov2_LS
 
%% tbt57_fov1
tbt57_fov1roinum = max(sort_tbt57_fov1_sig(:,1));

tbt57_fov1mat = zeros(tbt57_fov1roinum,136);

for i = 1:size(sort_tbt57_fov1_sig,1)
    roi=sort_tbt57_fov1_sig(i,1);
    odor=sort_tbt57_fov1_sig(i,2);
    tbt57_fov1mat(roi,odor)=sort_tbt57_fov1_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt57_fov1_real=tbt57_fov1mat(:, rep_odor);
tmpdata=tbt57_fov1_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt57fov1(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt57fov1); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt57fov1(:,deletecols) = [];

meantbt57_fov1_LS = mean(sparsenessval_tbt57fov1);
stdtbt57_fov1_LS = std(sparsenessval_tbt57fov1);
 
% clearvars -except tbt57_fov1_real sparsenessval_tbt57fov1 meantbt57_fov1_LS  stdtbt57_fov1_LS

%% tbt57_fov2
tbt57_fov2roinum = max(sort_tbt57_fov2_sig(:,1));

tbt57_fov2mat = zeros(tbt57_fov2roinum,136);

for i = 1:size(sort_tbt57_fov2_sig,1)
    roi=sort_tbt57_fov2_sig(i,1);
    odor=sort_tbt57_fov2_sig(i,2);
    tbt57_fov2mat(roi,odor)=sort_tbt57_fov2_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt57_fov2_real=tbt57_fov2mat(:, rep_odor);
tmpdata=tbt57_fov2_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt57fov2(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt57fov2); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt57fov2(:,deletecols) = [];

meantbt57_fov2_LS = mean(sparsenessval_tbt57fov2);
stdtbt57_fov2_LS = std(sparsenessval_tbt57fov2);
 
% clearvars -except tbt57_fov2_real sparsenessval_tbt57fov2 meantbt57_fov2_LS  stdtbt57_fov2_LS

%% tbt61_fov1
tbt61_lobfov1roinum = max(sort_tbt61_lobfov1_sig(:,1));

tbt61_lobfov1mat = zeros(tbt61_lobfov1roinum,136);

for i = 1:size(sort_tbt61_lobfov1_sig,1)
    roi=sort_tbt61_lobfov1_sig(i,1);
    odor=sort_tbt61_lobfov1_sig(i,2);
    tbt61_lobfov1mat(roi,odor)=sort_tbt61_lobfov1_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt61_lobfov1_real=tbt61_lobfov1mat(:, rep_odor);
tmpdata=tbt61_lobfov1_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt61lobfov1(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt61lobfov1); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt61lobfov1(:,deletecols) = [];

meantbt61_lobfov1_LS = mean(sparsenessval_tbt61lobfov1);
stdtbt61_lobfov1_LS = std(sparsenessval_tbt61lobfov1);
 
% clearvars -except tbt61_lobfov1_real sparsenessval_tbt61fov1 meantbt61_lobfov1_LS  stdtbt61_lobfov1_LS
%% tbt61_fov2
tbt61_lobfov2roinum = max(sort_tbt61_lobfov2_sig(:,1));

tbt61_lobfov2mat = zeros(tbt61_lobfov2roinum,136);

for i = 1:size(sort_tbt61_lobfov2_sig,1)
    roi=sort_tbt61_lobfov2_sig(i,1);
    odor=sort_tbt61_lobfov2_sig(i,2);
    tbt61_lobfov2mat(roi,odor)=sort_tbt61_lobfov2_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt61_lobfov2_real=tbt61_lobfov2mat(:, rep_odor);
tmpdata=tbt61_lobfov2_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt61lobfov2(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt61lobfov2); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt61lobfov2(:,deletecols) = [];

meantbt61_lobfov2_LS = mean(sparsenessval_tbt61lobfov2);
stdtbt61_lobfov2_LS = std(sparsenessval_tbt61lobfov2);
 
% clearvars -except tbt61_lobfov2_real sparsenessval_tbt61fov2 meantbt61_lobfov2_LS  stdtbt61_lobfov2_LS
%% tbt61_robfov1
tbt61_robroinum = max(sort_tbt61_rob_sig(:,1));

tbt61_robmat = zeros(tbt61_robroinum,136);

for i = 1:size(sort_tbt61_rob_sig,1)
    roi=sort_tbt61_rob_sig(i,1);
    odor=sort_tbt61_rob_sig(i,2);
    tbt61_robmat(roi,odor)=sort_tbt61_rob_sig(i,3);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 107 110 2 3 5 10 11 20 25 27 43 56 102 136];
tbt61_rob_real=tbt61_robmat(:, rep_odor);
tmpdata=tbt61_rob_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt61rob(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt61rob); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt61rob(:,deletecols) = [];

meantbt61_rob_LS = mean(sparsenessval_tbt61rob);
stdtbt61_rob_LS = std(sparsenessval_tbt61rob);
 
% clearvars -except tbt61_rob_real sparsenessval_tbt61rob meantbt61_rob_LS  stdtbt61_rob_LS
clearvars -except tbt51_fov1_real sparsenessval_tbt51fov1 meantbt51_fov1_LS stdtbt51_fov1_LS tbt51_fov2_real sparsenessval_tbt51fov2 meantbt51_fov2_LS  stdtbt51_fov2_LS tbt57_fov1_real sparsenessval_tbt57fov1 meantbt57_fov1_LS  stdtbt57_fov1_LS tbt57_fov2_real sparsenessval_tbt57fov2 meantbt57_fov2_LS  stdtbt57_fov2_LS tbt61_lobfov1_real sparsenessval_tbt61fov1 meantbt61_lobfov1_LS stdtbt61_lobfov1_LS tbt61_lobfov2_real sparsenessval_tbt61fov2 meantbt61_lobfov2_LS  stdtbt61_lobfov2_LS tbt61_rob_real sparsenessval_tbt61rob meantbt61_rob_LS  stdtbt61_rob_LS

% % allfov = horzcat(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);
% % allmean = mean(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);
% % allstd = std(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);