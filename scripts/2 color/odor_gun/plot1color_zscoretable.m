clearvars
close all 
% clc

%% load results files
[twocfile,twocpath] = uigetfile('*.mat','Choose Results file for twocolor imaging','MultiSelect','off');
if ~twocfile; return; end
load(twocfile); %clear twocfile twocpath;

iglumat_ex = results(:,:,1);
igluexmat = isnan(iglumat_ex);
iglui = find(igluexmat~=1);
for i=1:length(iglui)
    igluexmat_num(i) = iglumat_ex(iglui(i));
end

clearvars iglui

iglumat_sup = results(:,:,2);
iglusupmat = isnan(iglumat_sup);
iglui = find(iglusupmat~=1);
for i=1:length(iglui)
    iglusupmat_num(i) = iglumat_sup(iglui(i));
end

%% now compare arrays

igluexmat_ind = isfinite(iglumat_ex); %pull out any that isn't nan
iglusupmat_ind = isfinite(iglumat_sup);


igluexmat_ind = find(igluexmat_ind); %find location of values for all finite values
iglusupmat_ind = find(iglusupmat_ind);


[Csect,isecta,isectb] = intersect(igluexmat_ind,iglusupmat_ind); %do the locations overlap?
[Cdiff,idiff] = setdiff(igluexmat_ind,iglusupmat_ind);


iglusnfr_onlyex =iglumat_ex(Cdiffa); %iglusnfr only
iglusnfr_onlysup =iglumat_sup(Cdiffa);

%clearvars -except iglusnfr_onlyex iglusnfr_onlysup iglusnfr_jrgecoex iglusnfr_jrgecosup jrgeco_onlyex jrgeco_onlysup % use for clearing but don't keep code here
 
%now to give tables of solely excited, suppressed, and coactivated gloms
 
a = numel(iglusnfr_onlyex);
b = numel(iglusnfr_onlysup);

rois =  (size(results,1)-1);
odors = 24;
sumallgo = (rois*odors);

%% start table
iglusnfr_percexcited = [a/sumallgo];
iglusnfr_percsuppressed = [b/sumallgo];

% exc = [iglusnfr_percexcited,rgeco_percexcited,iglusnfrgeco_percexcited];
% sup = [iglusnfr_percsuppressed,rgeco_percsuppressed,iglusnfrgeco_percsuppressed];
% 
% T = table(exc,sup);
T = table(iglusnfr_percexcited,iglusnfr_percsuppressed);

