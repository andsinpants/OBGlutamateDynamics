clearvars
close all
load('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SFiGluSnFR_jRGECO1a\tbt_jRGECO_SFiGluSnFR ogun only\tbt_jRGECO_SFiGluSnFR_ogun analysis\allexp_resultssig_vars_2color_forsparsenessanalysis')

sort_tbt68_sig2 = sortrows(results_sig2_tbt68,2);
sort_tbt78_robfov1_sig2 = sortrows(results_sig2_tbt78_robfov1,2);
sort_tbt78_robfov2_sig2 = sortrows(results_sig2_tbt78_robfov2,2);
sort_tbt78_lobfov1_sig2 = sortrows(results_sig2_tbt78_lobfov1,2);
sort_tbt82_sig2 = sortrows(results_sig2_tbt82,2);


%% tbt68_fov
tbt68_roinum = max(sort_tbt68_sig2(:,1));

tbt68_mat = zeros(tbt68_roinum,136);

for i = 1:size(sort_tbt68_sig2,1)
    roi=sort_tbt68_sig2(i,1);
    odor=sort_tbt68_sig2(i,2);
    tbt68_mat(roi,odor)=sort_tbt68_sig2(i,4);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136];
tbt68_real=tbt68_mat(:, rep_odor)'; %% transpose back to calculate population sparseness;
tmpdata=tbt68_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt68(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt68); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt68(:,deletecols) = [];

meantbt68_LS = mean(sparsenessval_tbt68);
stdtbt68_LS = std(sparsenessval_tbt68);
 
% clearvars -except tbt68_real sparsenessval_tbt68 meantbt68_LS  stdtbt68_LS
 
%% tbt78_robfov1
tbt78_robfov1roinum = max(sort_tbt78_robfov1_sig2(:,1));

tbt78_robfov1mat = zeros(tbt78_robfov1roinum,136);

for i = 1:size(sort_tbt78_robfov1_sig2,1)
    roi=sort_tbt78_robfov1_sig2(i,1);
    odor=sort_tbt78_robfov1_sig2(i,2);
    tbt78_robfov1mat(roi,odor)=sort_tbt78_robfov1_sig2(i,4);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136];
tbt78_robfov1_real=tbt78_robfov1mat(:, rep_odor)'; %% transpose back to calculate population sparseness;
tmpdata=tbt78_robfov1_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt78robfov1(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt78robfov1); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt78robfov1(:,deletecols) = [];

meantbt78_robfov1_LS = mean(sparsenessval_tbt78robfov1);
stdtbt78_robfov1_LS = std(sparsenessval_tbt78robfov1);
 
% clearvars -except tbt78_robfov1_real sparsenessval_tbt78robfov1 meantbt78_robfov1_LS  stdtbt78_robfov1_LS
 
%% tbt78_robfov2
tbt78_robfov2roinum = max(sort_tbt78_robfov2_sig2(:,1));

tbt78_robfov2mat = zeros(tbt78_robfov2roinum,136);

for i = 1:size(sort_tbt78_robfov2_sig2,1)
    roi=sort_tbt78_robfov2_sig2(i,1);
    odor=sort_tbt78_robfov2_sig2(i,2);
    tbt78_robfov2mat(roi,odor)=sort_tbt78_robfov2_sig2(i,4);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136];
tbt78_robfov2_real=tbt78_robfov2mat(:, rep_odor)'; %% transpose back to calculate population sparseness;
tmpdata=tbt78_robfov2_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt78robfov2(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt78robfov2); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt78robfov2(:,deletecols) = [];

meantbt78_robfov2_LS = mean(sparsenessval_tbt78robfov2);
stdtbt78_robfov2_LS = std(sparsenessval_tbt78robfov2);
 
% clearvars -except tbt78_robfov2_real sparsenessval_tbt78robfov2 meantbt78_robfov2_LS  stdtbt78_robfov2_LS

%% tbt78_lobfov1
tbt78_lobfov1roinum = max(sort_tbt78_lobfov1_sig2(:,1));

tbt78_lobfov1mat = zeros(tbt78_lobfov1roinum,136);

for i = 1:size(sort_tbt78_lobfov1_sig2,1)
    roi=sort_tbt78_lobfov1_sig2(i,1);
    odor=sort_tbt78_lobfov1_sig2(i,2);
    tbt78_lobfov1mat(roi,odor)=sort_tbt78_lobfov1_sig2(i,4);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136];
tbt78_lobfov1_real=tbt78_lobfov1mat(:, rep_odor)'; %% transpose back to calculate population sparseness;
tmpdata=tbt78_lobfov1_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt78lobfov1(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt78lobfov1); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt78lobfov1(:,deletecols) = [];

meantbt78_lobfov1_LS = mean(sparsenessval_tbt78lobfov1);
stdtbt78_lobfov1_LS = std(sparsenessval_tbt78lobfov1);
 
% clearvars -except tbt78_lobfov1_real sparsenessval_tbt78lobfov1 meantbt78_lobfov1_LS  stdtbt78_lobfov1_LS

%% tbt82_robfov1
tbt82_roinum = max(sort_tbt82_sig2(:,1));

tbt82_mat = zeros(tbt82_roinum,136);

for i = 1:size(sort_tbt82_sig2,1)
    roi=sort_tbt82_sig2(i,1);
    odor=sort_tbt82_sig2(i,2);
    tbt82_mat(roi,odor)=sort_tbt82_sig2(i,4);
end

rep_odor = [1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136];
tbt82_real=tbt82_mat(:, rep_odor)'; %% transpose back to calculate population sparseness;
tmpdata=tbt82_real;

for o = 1:size(tmpdata,2)
    sparsenessval_tbt82(o) = (1- (sum(tmpdata(:,o)./size(tmpdata,1))^2/sum(tmpdata(:,o).^2./size(tmpdata,1))) )/(1-1/size(tmpdata,1));
end

checksparse = isnan(sparsenessval_tbt82); %remove NaN numbers from LS calculation
deletecols = find(checksparse);
sparsenessval_tbt82(:,deletecols) = [];

meantbt82_LS = mean(sparsenessval_tbt82);
stdtbt82_LS = std(sparsenessval_tbt82);
 
% clearvars -except tbt82_real sparsenessval_tbt82 meantbt82_LS  stdtbt82_LS

 
clearvars -except tbt68_real sparsenessval_tbt68 meantbt68_LS  stdtbt68_LS tbt78_robfov1_real sparsenessval_tbt78robfov1 meantbt78_robfov1_LS  stdtbt78_robfov1_LS tbt78_robfov2_real sparsenessval_tbt78robfov2 meantbt78_robfov2_LS  stdtbt78_robfov2_LS tbt78_lobfov1_real sparsenessval_tbt78lobfov1 meantbt78_lobfov1_LS  stdtbt78_lobfov1_LS tbt82_real sparsenessval_tbt82 meantbt82_LS  stdtbt82_LS


% % allfov = horzcat(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);
% % allmean = mean(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);
% % allstd = std(meantbt51_fov1_LS,meantbt51_fov2_LS,meantbt57_fov1_LS,meantbt57_fov2_LS,meantbt61_lobfov1_LS,meantbt61_lobfov2_LS,meantbt61_rob_LS);