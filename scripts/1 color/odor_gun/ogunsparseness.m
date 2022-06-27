% clearvars
close all
load('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\odor_gun\SFiGluSnFR\tbt_SFiGluSnFR_ogun analysis\allexp_resultssig_vars_1color_forsparsenessanalysis')

tbt51 = vertcat(results_sig_tbt51_fov1,results_sig_tbt51_fov2);
tbt57 = vertcat(results_sig_tbt57_fov1,results_sig_tbt57_fov2);
tbt61 = vertcat(results_sig_tbt61_lobfov1,results_sig_tbt61_lobfov2,results_sig_tbt61_rob);

rep_odor = {1 6 9 13 23 24 28 33 40 45 51 81 2 3 5 10 11 20 25 27 43 56 102 136};
rep_expt = [{tbt51};{tbt57};{tbt61}]; 
N_roi = 748;
z_exc = allresults_sig(:,3);
z_sup = allresults_sig(:,4);

sort_tbt51_fov1_sig = sortrows(results_sig_tbt51_fov1,2);
sort_tbt51_fov2_sig = sortrows(results_sig_tbt51_fov2,2);
sort_tbt57_fov1_sig = sortrows(results_sig_tbt57_fov1,2);
sort_tbt61_lobfov1_sig = sortrows(results_sig_tbt61_lobfov1,2);
sort_tbt61_lobfov2_sig = sortrows(results_sig_tbt61_lobfov2,2);
sort_tbt61_rob_sig = sortrows(results_sig_tbt61_rob,2);

tbt51_fov1roinum = max(sort_tbt51_fov1_sig(:,1));
tbt51_fov1odornum = max(sort_tbt51_fov1_sig(:,2));

tbt51_fov1mat = zeros(tbt51_fov1roinum,tbt51_fov1odornum);

for i = 1:size(tbt51_fov1mat,2)
    for j = 1:size(tbt51_fov1mat,1)
    tbt51_fov1mat(i,j) = tbt51_fov1mat(i,2);
    end
end

 
 for r=1:length(rep_odor)  %%% for particular odorants interested in
%     for r=1:size(allresults_sig,1)
        for E=1:size(rep_expt{r},1) %%% for particular experiments interested in
        e={rep_expt{r}(E)}; %%% index of experiment number
        e = mat2cell(e);
        %%% begin actual sparseness calculation
        int1=0;
        int2=0;
        for g=1:N_roi(e) %%% cycle through ROIs
            if (z_exc{e}(g,rep_odor(r)))>0 %%% negative responses not included in this calculation
                int1=int1+z_exc{e}(g,rep_odor(r))/N_roi(e);
                int2=int2+z_exc{e}(g,rep_odor(r)^2)/N_roi(e);
%             elseif z_sup{e}(g,rep_odor(r))<0 %%added by AKM
%                 int3=int3+z_sup{e}(g,rep_odor(r))/N_roi(e);
%                 int4=int4+(z_sup{e}(g,rep_odor(r))^2)/N_roi(e);
            end
        end
        LS{r}(E)=(1-int1^2/int2)/(1-(1/N_roi(e)));
        end
    mean_LS(r)=mean(LS{r});
    std_LS(r)=std(LS{r});
end