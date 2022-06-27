for r=1:length(rep_odor)  %%% for particular odorants interested in
    for E=1:length(rep_expt{r}) %%% for particular experiments interested in
        e=rep_expt{r}(E); %%% index of experiment number

        %%% begin actual sparseness calculation
        int1=0;
        int2=0;
        for g=1:N_roi(e) %%% cycle through ROIs
            if z{e}(g,rep_odor(r))>0 %%% negative responses not included in this calculation
                int1=int1+z{e}(g,rep_odor(r))/N_roi(e);
                int2=int2+(z{e}(g,rep_odor(r))^2)/N_roi(e);
            end
        end
        LS{r}(E)=(1-int1^2/int2)/(1-(1/N_roi(e)));
    end
    mean_LS(r)=mean(LS{r});
    std_LS(r)=std(LS{r});
end