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

jrgecomat_ex = results(:,:,3);
jrgecoexmat = isnan(jrgecomat_ex);
jrgecoi = find(jrgecoexmat~=1);
for i=1:length(jrgecoi)
    jrgecoexmat_num(i) = jrgecomat_ex(jrgecoi(i));
end

clearvars jrgecoi

jrgecomat_sup = results(:,:,4);
jrgecosupmat = isnan(jrgecomat_sup);
jrgecoi = find(jrgecosupmat~=1);
for i=1:length(jrgecoi)
    jrgecosupmat_num(i) = jrgecomat_sup(jrgecoi(i));
end

hold on; legend();

%% now compare arrays

igluexmat_ind = isfinite(iglumat_ex); %pull out any that isn't nan
iglusupmat_ind = isfinite(iglumat_sup);
jregcoexmat_ind = isfinite(jrgecomat_ex);
jrgecosupmat_ind = isfinite(jrgecomat_sup);

igluexmat_ind = find(igluexmat_ind); %find location of values for all finite values
iglusupmat_ind = find(iglusupmat_ind);
jrgecoexmat_ind = find(jregcoexmat_ind);
jrgecosupmat_ind = find(jrgecosupmat_ind);

[Cexsect,iexsecta,iexsectb] = intersect(igluexmat_ind,jrgecoexmat_ind); %do the locations overlap?
[Csupsect,isupsecta,isupsectb] = intersect(iglusupmat_ind,jrgecosupmat_ind);
[Cexdiffa,iexdiffa] = setdiff(igluexmat_ind,jrgecoexmat_ind);
[Csupdiffa,isupdiffa] = setdiff(iglusupmat_ind,jrgecosupmat_ind);
[Cexdiffb,iexdiffb] = setdiff(jrgecoexmat_ind,igluexmat_ind);
[Csupdiffb,isupdiffb] = setdiff(jrgecosupmat_ind,iglusupmat_ind);

iglusnfr_onlyex =iglumat_ex(Cexdiffa); %iglusnfr only
iglusnfr_onlysup =iglumat_sup(Csupdiffa);

iglucatex =  iglumat_ex(Cexsect);
jrgecocatex = jrgecomat_ex(Cexsect);
iglusnfr_jrgecoex = horzcat(iglucatex,jrgecocatex); %both
igluexcatsup =  iglumat_sup(Csupsect);
jrgecocatsup = jrgecomat_sup(Csupsect);
iglusnfr_jrgecosup = horzcat(igluexcatsup,jrgecocatsup);

jrgeco_onlyex = jrgecomat_ex(Cexdiffb);%jrgeco only
jrgeco_onlysup = jrgecomat_sup(Csupdiffb);

%clearvars -except iglusnfr_onlyex iglusnfr_onlysup iglusnfr_jrgecoex iglusnfr_jrgecosup jrgeco_onlyex jrgeco_onlysup % use for clearing but don't keep code here
 
%% now to give unity plots of coexcited reseponses
figure;scatter(iglusnfr_jrgecoex(:,1),iglusnfr_jrgecoex(:,2));title('scatter of significant excited z-scores across Ch1 and Ch2');xlabel('Ch1:iGluSnFR');ylabel('Ch2:jRGeCO');
 h = gca;
% h.YLim = h.XLim;
x = h.XLim;
y = h.YLim;
line(x,y,'Color','red','LineStyle','--');

%% now to give unity plots of cosuppressed responses
for i=1
    if ~isempty(iglusnfr_jrgecosup)
        figure;scatter(iglusnfr_jrgecosup(:,1),iglusnfr_jrgecosup(:,2));title('scatter of significant suppressed z-scores across Ch1 and Ch2');xlabel('Ch1:iGluSnFR');ylabel('Ch2:jRGeCO');
        h = gca;
%         h.XLim = h.YLim;
        x = h.XLim;
        y = h.YLim;
        line(x,y,'Color','red','LineStyle','--');
    else
        break
    end
end

%%now to give tables of solely excited, suppressed, and coactivated gloms
 
a = numel(iglusnfr_onlyex);
b = numel(iglusnfr_onlysup);
c = numel(jrgeco_onlyex);
d = numel(jrgeco_onlysup);
e = numel(iglusnfr_jrgecoex);
for i =1
    if size(iglusnfr_jrgecoex,2)~=0
        e = numel(iglusnfr_jrgecoex(:,1));
    end
end
f = numel(iglusnfr_jrgecosup);
for i =1
    if size(iglusnfr_jrgecosup,2)~=0
        f = numel(iglusnfr_jrgecosup(:,1));
    end
end
rois =  (size(results,1)-1);
odors = 24;
sumallgo = (rois*odors);

%% start table
iglusnfr_percexcited = [a/sumallgo];
iglusnfr_percsuppressed = [b/sumallgo];
iglusnfrgeco_percexcited = [e/sumallgo];
iglusnfrgeco_percsuppressed = [f/sumallgo];
rgeco_percexcited = [c/sumallgo];
rgeco_percsuppressed = [d/sumallgo];

% exc = [iglusnfr_percexcited,rgeco_percexcited,iglusnfrgeco_percexcited];
% sup = [iglusnfr_percsuppressed,rgeco_percsuppressed,iglusnfrgeco_percsuppressed];
% 
% T = table(exc,sup);
T = table(iglusnfr_percexcited,iglusnfr_percsuppressed,rgeco_percexcited,rgeco_percsuppressed,iglusnfrgeco_percexcited,iglusnfrgeco_percsuppressed)

%% now to get rows and col sigvars
vals=[];
ch1ex_fin = isfinite(results(:,:,1));
[row,col] = find(ch1ex_fin==1);
ch1ex_cat = horzcat(row,col);
for r = 1:size(ch1ex_cat,1)
    vals(r) = results(ch1ex_cat(r,1),ch1ex_cat(r,2),1);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch1ex_cat = horzcat(row,col,vals);
ch1ex_sortcat = sortrows(ch1ex_cat);
clearvars row col vals
    else
    end
end

vals = [];
ch1sup_fin = isfinite(results(:,:,2));
[row,col] = find(ch1sup_fin==1);
ch1sup_cat = horzcat(row,col);
for r = 1:size(ch1sup_cat,1)
    vals(r) = results(ch1sup_cat(r,1),ch1sup_cat(r,2),2);
end
vals = vals';
for i=1
    if ~isempty(vals)
        ch1sup_cat = horzcat(row,col,vals);
        ch1sup_sortcat = sortrows(ch1sup_cat);
        clearvars row col vals
    else 
    end
end


vals = [];
ch2ex_fin = isfinite(results(:,:,3));
[row,col] = find(ch2ex_fin==1);
ch2ex_cat = horzcat(row,col);
for r = 1:size(ch2ex_cat,1)
    vals(r) = results(ch2ex_cat(r,1),ch2ex_cat(r,2),3);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch2ex_cat = horzcat(row,col,vals);
ch2ex_sortcat = sortrows(ch2ex_cat);
clearvars row col vals
    else
    end
end

vals = [];
ch2sup_fin = isfinite(results(:,:,4));
[row,col] = find(ch2sup_fin==1);
ch2sup_cat = horzcat(row,col);
for r = 1:size(ch2sup_cat,1)
    vals(r) = results(ch2sup_cat(r,1),ch2sup_cat(r,2),4);
end
vals = vals';
for i=1
    if ~isempty(vals)
ch2sup_cat = horzcat(row,col,vals);
ch2sup_sortcat = sortrows(ch2sup_cat);
clearvars row col vals
    else
    end
end
 
openvar('ch1ex_cat')
openvar('ch2sup_cat')
openvar('ch2ex_cat')
openvar('ch2sup_cat')

clearvars -except ch1ex_cat ch1sup_cat ch2ex_cat ch2sup_cat ch1ex_sortcat ch1sup_sortcat ch2ex_sortcat ch2sup_sortcat results T

for i=1
if ~isempty(ch1ex_cat) && ~isempty(ch1sup_cat) && ~isempty(ch2ex_cat) && ~isempty(ch2sup_cat)
figure;scatter3(ch1ex_cat(:,1),ch1ex_cat(:,2),ch1ex_cat(:,3),'MarkerFaceColor','g')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch1sup_cat(:,1),ch1sup_cat(:,2),ch1sup_cat(:,3),'MarkerFaceColor','b')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Top Intensity (zscored dF values)');

hold 'on';scatter3(ch2ex_cat(:,1),ch2ex_cat(:,2),ch2ex_cat(:,3),'MarkerFaceColor','r')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch2sup_cat(:,1),ch2sup_cat(:,2),ch2sup_cat(:,3),'MarkerFaceColor','m')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');
legend('iGluSnFR +','iGluSnFR -','jRGECO +','jRECO -');
h = gca;
end
end

for i=1
    if isempty(ch2sup_cat)
figure;scatter3(ch1ex_cat(:,1),ch1ex_cat(:,2),ch1ex_cat(:,3),'MarkerFaceColor','g')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch1sup_cat(:,1),ch1sup_cat(:,2),ch1sup_cat(:,3),'MarkerFaceColor','b')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Top Intensity (zscored dF values)');

hold 'on';scatter3(ch2ex_cat(:,1),ch2ex_cat(:,2),ch2ex_cat(:,3),'MarkerFaceColor','r')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');
    end
end

for i=1
    if isempty(ch1sup_cat)
figure;scatter3(ch1ex_cat(:,1),ch1ex_cat(:,2),ch1ex_cat(:,3),'MarkerFaceColor','g')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch2ex_cat(:,1),ch2ex_cat(:,2),ch2ex_cat(:,3),'MarkerFaceColor','r')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch2sup_cat(:,1),ch2sup_cat(:,2),ch2sup_cat(:,3),'MarkerFaceColor','m')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');
    end
end

for i=1
    if isempty(ch1sup_cat) && isempty(ch2sup_cat)
figure;scatter3(ch1ex_cat(:,1),ch1ex_cat(:,2),ch1ex_cat(:,3),'MarkerFaceColor','g')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

hold 'on';scatter3(ch2ex_cat(:,1),ch2ex_cat(:,2),ch2ex_cat(:,3),'MarkerFaceColor','r')
xlabel('Glomerulus');
ylabel('Odor');
zlabel('Intensity (zscored dF values)');

    end
end

clear i 
clearvars

%%load in complete dataset
figure;histogram(ch1ex_sortcat(:,3),'FaceAlpha',0.25,'FaceColor','g');
hold on;
histogram(ch2ex_sortcat(:,3),'FaceAlpha',0.25,'FaceColor','r');
histogram(ch1sup_sortcat(:,3),'FaceAlpha',0.25,'FaceColor','b');
histogram(ch2sup_sortcat(:,3),'FaceAlpha',0.25,'FaceColor','m');
ylabel('Counts')
xlabel('Z-scored peak and trough dF values >= 7');
legend('iGluSnFR +','jRGECO +','iGluSnFR -','jRECO -');


figure; cdfplot(ch1ex_sortcat(:,3))
hold on;cdfplot(ch1sup_sortcat(:,3))
hold on;cdfplot(ch2ex_sortcat(:,3))
hold on;cdfplot(ch2sup_sortcat(:,3))
legend('iGluSnFR +','iGluSnFR -','jRGECO +','jRECO -');
xlabel('Z-scored peak and trough dF values >= 7');
ylabel('Counts');


%load timeseries analysis
[twocfile,twocpath] = uigetfile('*.mat','Choose Results file for twocolor imaging','MultiSelect','off');
if ~twocfile; return; end
load(twocfile); %clear twocfile twocpath;




% % h.YTicks = [];
% % h.YTickLabel = [];
% % h.YTick = [];
%T_sigresp = table(ch1ex_cat,ch1sup_cat,ch2ex_cat,ch2sup_cat)';



