clc
close all
%%load
cd('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\pcd iGluSnFR 3x odors\Total Dataset Analysis\Analysis .mat files\Updated Analysis')
load('pcd_conc_results_t2-t1');

%% pcd300
figure(1)
for r=1:size(pcd300_ET_T,1)
    plot([1 2 3],[pcd300_ET_T(r,1) pcd300_ET_T(r,2) pcd300_ET_T(r,3)],'--x');hold on;title('ethyl tiglate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(2)
for r=1:length(pcd300_MV_T)
plot([1 2 3],[pcd300_MV_T(r,1) pcd300_MV_T(r,2) pcd300_MV_T(r,3)],'--x');hold on;title('methyl valerate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(3)
for r=1:length(pcd300_EB_T)
plot([1 2 3],[pcd300_EB_T(r,1) pcd300_EB_T(r,2) pcd300_EB_T(r,3)],'--x');hold on;title('ethyl butyrate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(4)
for r=1:size(pcd300_MB_T,1)
    plot([1 2 3],[pcd300_MB_T(r,1) pcd300_MB_T(r,2) pcd300_MB_T(r,3)],'--x');hold on;title('methyl benzoate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

%% pcd310
figure(1)
for r=1:size(pcd310_ET_T,1)
    plot([1 2 3],[pcd310_ET_T(r,1) pcd310_ET_T(r,2) pcd310_ET_T(r,3)],'--x');hold on;title('ethyl tiglate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(2)
for r=1:length(pcd310_MV_T)
plot([1 2 3],[pcd310_MV_T(r,1) pcd310_MV_T(r,2) pcd310_MV_T(r,3)],'--x');hold on;title('methyl valerate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(3)
for r=1:length(pcd310_2hex_T)
plot([1 2 3],[pcd310_2hex_T(r,1) pcd310_2hex_T(r,2) pcd310_2hex_T(r,3)],'--x');hold on;title('2hexanone (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

figure(4)
for r=1:size(pcd310_MB_T,1)
    plot([1 2 3],[pcd310_MB_T(r,1) pcd310_MB_T(r,2) pcd310_MB_T(r,3)],'--x');hold on;title('methyl benzoate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end
%% pcd314
figure(1)
for r=1:size(pcd314_ET_T,1)
    plot([1 2 3],[pcd314_ET_T(r,1) pcd314_ET_T(r,2) pcd314_ET_T(r,3)],'--x');hold on;title('ethyl tiglate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

% figure(2);subplot(1,3,1);imagesc(pcd314_ET_T(:,1));xlabel('1x Conc.'),ylabel(' Ethyl Tiglate T2-T1 Value')
% subplot(1,3,2);imagesc(pcd314_ET_T(:,2));xlabel('3x Conc.'),ylabel(' Ethyl Tiglate T2-T1 Value')
% subplot(1,3,3);imagesc(pcd314_ET_T(:,3));xlabel('10x Conc.'),ylabel(' Ethyl Tiglate T2-T1 Value')

figure(3)
for r=1:length(pcd314_MV_T)
plot([1 2 3],[pcd314_MV_T(r,1) pcd314_MV_T(r,2) pcd314_MV_T(r,3)],'--x');hold on;title('methyl valerate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

% figure(4);subplot(1,3,1);imagesc(pcd314_MV_T(:,1));xlabel('1x Conc.'),ylabel(' methyl valerate T2-T1 Value')
% subplot(1,3,2);imagesc(pcd314_MV_T(:,2));xlabel('3x Conc.'),ylabel(' methyl valerate T2-T1 Value')
% subplot(1,3,3);imagesc(pcd314_MV_T(:,3));xlabel('10x Conc.'),ylabel(' methyl valerate T2-T1 Value')

figure(5)
for r=1:length(pcd314_2hex_T)
plot([1 2 3],[pcd314_2hex_T(r,1) pcd314_2hex_T(r,2) pcd314_2hex_T(r,3)],'--x');hold on;title('2hexanone (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

% figure(6);subplot(1,3,1);imagesc(pcd314_2hex_T(:,1));xlabel('1x Conc.'),ylabel('2hexanone T2-T1 Value')
% subplot(1,3,2);imagesc(pcd314_2hex_T(:,2));xlabel('3x Conc.'),ylabel('2hexanone T2-T1 Value')
% subplot(1,3,3);imagesc(pcd314_2hex_T(:,3));xlabel('10x Conc.'),ylabel(' 2hexanone T2-T1 Value')

figure(7)
for r=1:size(pcd314_EB_T,1)
    plot([1 2 3],[pcd314_EB_T(r,1) pcd314_EB_T(r,2) pcd314_EB_T(r,3)],'--x');hold on;title('ethyl butyrate (t2-t1)');xlabel('conc');ylabel('T2-T1 Value');
end

% figure(8);subplot(1,3,1);imagesc(pcd314_EB_T(:,1));xlabel('1x Conc.'),ylabel('Ethyl Butyrate T2-T1 Value')
% subplot(1,3,2);imagesc(pcd314_EB_T(:,2));xlabel('3x Conc.'),ylabel('Ethyl Butyrate T2-T1 Value')
% subplot(1,3,3);imagesc(pcd314_EB_T(:,3));xlabel('10x Conc.'),ylabel('Ethyl Butyrate T2-T1 Value')
% subplot(1,3,3);imagesc(pcd314_EB_T(:,3));xlabel('Conc.'),ylabel('T2-T1 Value')

%% now concantenated
clearvars
clc
cd('C:\Users\matt\Dropbox\mwlab_andrew\MATLAB_home\MATLAB_home_sub\MATLAB\Current Datasets and Scripts\Datasets and Scripts\main_datasets\2P\pcd iGluSnFR 3x odors\Total Dataset Analysis\Analysis .mat files\Updated Analysis');
load('pcd_conc_results_vertcat');

%% pcd273
load('pcd_conc_results_t2-t1.mat')
pcd273_2hex_1x = pcd273_2hex_T(:,1);
pcd273_2hex_10x = pcd273_2hex_T(:,3);
pcd273_EB_1x = pcd273_EB_T(:,1);
pcd273_EB_10x = pcd273_EB_T(:,3);
pcd273_MV_1x = pcd273_MV_T(:,1);
pcd273_MV_10x = pcd273_MV_T(:,3);
figure;scatter(pcd273_2hex_1x,pcd273_2hex_10x);hold on;scatter(pcd273_EB_1x,pcd273_EB_10x);scatter(pcd273_MV_1x,pcd273_MV_10x);
refline
pcd273_1x = vertcat(pcd273_2hex_1x,pcd273_MV_1x,pcd273_EB_1x);
pcd273_10x = vertcat(pcd273_2hex_10x,pcd273_MV_10x,pcd273_EB_10x);
figure;scatter(pcd273_1x,pcd273_10x);xlabel('1x T2-T1 all odorants');ylabel('10x T2-T1 all odorants');refline
mu = mean(pcd273_10x);
hline = refline([0 mu]);
hline.Color = 'r';

%% pcd300
pcd300_EB_1x = pcd300_EB_T(:,1);
pcd300_ET_1x = pcd300_ET_T(:,1);
pcd300_MB_1x = pcd300_MB_T(:,1);
pcd300_MV_1x = pcd300_MV_T(:,1);
pcd300_EB_10x = pcd300_EB_T(:,3);
pcd300_ET_10x = pcd300_ET_T(:,3);
pcd300_MB_10x = pcd300_MB_T(:,3);
pcd300_MV_10x = pcd300_MV_T(:,3);
pcd300_1x = vertcat(pcd300_EB_1x,pcd300_ET_1x,pcd300_MB_1x,pcd300_MV_1x);
pcd300_10x = vertcat(pcd300_EB_10x,pcd300_ET_10x,pcd300_MB_10x,pcd300_MV_10x);
figure;scatter(pcd300_1x,pcd300_10x);xlabel('1x T2-T1 all odorants');ylabel('10x T2-T1 all odorants');refline;mu = mean(pcd300_10x);
hline = refline([0 mu]);
hline.Color = 'r';
%% pcd310
pcd310_2hex_1x = pcd310_2hex_T(:,1);
pcd310_ET_1x = pcd310_ET_T(:,1);
pcd310_MB_1x = pcd310_MB_T(:,1);
pcd310_MV_1x = pcd310_MV_T(:,1);
pcd310_MV_10x = pcd310_MV_T(:,3);
pcd310_MB_10x = pcd310_MB_T(:,3);
pcd310_ET_10x = pcd310_ET_T(:,3);
pcd310_2hex_10x = pcd310_2hex_T(:,3);
pcd310_1x = vertcat(pcd310_2hex_1x,pcd310_ET_1x,pcd310_MB_1x,pcd310_MV_1x);
pcd310_10x = vertcat(pcd310_2hex_10x,pcd310_ET_10x,pcd310_MB_10x,pcd310_MV_10x);
figure;scatter(pcd310_1x,pcd310_10x);xlabel('1x T2-T1 all odorants');ylabel('10x T2-T1 all odorants');refline;mu = mean(pcd310_10x);
hline = refline([0 mu]);
hline.Color = 'r';
%% pcd314
%does not matter if all are included in the dataset as we are only using 1x and 10x concentrations
pcd314_2hex_1x = pcd314_2hex_T(:,1);
pcd314_ET_1x = pcd314_ET_T(:,1);
pcd314_MV_1x = pcd314_MV_T(:,1);
pcd314_MV_10x = pcd314_MV_T(:,3);
pcd314_ET_10x = pcd314_ET_T(:,3);
pcd314_2hex_10x = pcd314_2hex_T(:,3);
pcd314_1x = vertcat(pcd314_2hex_1x,pcd314_ET_1x,pcd314_MV_1x);
pcd314_10x = vertcat(pcd314_2hex_10x,pcd314_ET_10x,pcd314_MV_10x);
figure;scatter(pcd314_1x,pcd314_10x);xlabel('1x T2-T1 all odorants');ylabel('10x T2-T1 all odorants');refline;mu = mean(pcd314_10x);
hline = refline([0 mu]);
hline.Color = 'r';
%% now graph all together and make a refline from origin
allexp_1x = vertcat(pcd273_1x,pcd300_1x,pcd310_1x,pcd314_1x);
allexp_10x = vertcat(pcd273_10x,pcd300_10x,pcd310_10x,pcd314_10x);
update = find(allexp_1x>=-1);
allexp_1x_update = (allexp_1x(update));
allexp_10x_update = (allexp_10x(update));
figure;scatter(pcd273_1x,pcd273_10x,'x','r');hold on;
scatter(pcd300_1x,pcd300_10x,'x','k');hold on;
scatter(pcd310_1x,pcd310_10x,'x','b');hold on;
scatter(pcd314_1x,pcd314_10x,'x','g');hold on;
% gcf
% xlim = ylim;
% xtick = ticks;
figure;scatter(allexp_1x_update,allexp_10x_update);
line([1 1],[allexp_10x_update(1) allexp_10x_update(end)],'LineWidth',10)%double check if working
%%perform summary statistics 
[h,p,ci,stats]  = ttest(allexp_1x_update,allexp_10x_update);
%% now graph with colors for odorants, ignoring experiments
load('pcd_conc_results_t2-t1.mat')
hex_concant = vertcat(pcd273_2hex_T,pcd310_2hex_T,pcd314_2hex_T);
EB_concant = vertcat(pcd273_EB_T,pcd300_EB_T);
MV_concant = vertcat(pcd273_MV_T,pcd300_MV_T,pcd310_MV_T,pcd314_MV_T);
ET_concant = vertcat(pcd300_ET_T,pcd310_ET_T,pcd314_ET_T);
MB_concant = vertcat(pcd300_MB_T,pcd310_MB_T);
%now for concentrations
EB_1x = EB_concant(:,1);
EB_10x = EB_concant(:,3);

hex_1x = hex_concant(:,1);
hex_10x = hex_concant(:,3);

ET_1x = ET_concant(:,1);
ET_10x = ET_concant(:,3);

MV_1x = MV_concant(:,1);
MV_10x = MV_concant(:,3);

MB_1x = MB_concant(:,1);
MB_10x = MB_concant(:,3);

figure;scatter(EB_1x,EB_10x,'x','r');hold on;
scatter(ET_1x,ET_10x,'x','k');hold on;
scatter(hex_1x,hex_10x,'x','b');hold on;
scatter(MB_1x,MB_10x,'x','g');hold on;
scatter(MV_1x,MV_10x,'x','m');
legend('EB','ET','hex','MB','MV');hold off

% update = find(allexp_1x>=-1);
% allexp_1x_update = (allexp_1x(update));
% allexp_10x_update = (allexp_10x(update));
% figure;scatter(pcd273_1x,pcd273_10x,'x','r');hold on;
% scatter(pcd300_1x,pcd300_10x,'x','k');hold on;
% scatter(pcd310_1x,pcd310_10x,'x','b');hold on;
% scatter(pcd314_1x,pcd314_10x,'x','g');hold on;
%% now to plot the df/f T1 max by all t2-t1 values
clearvars
load('pcd_conc_results_t2-t1_t1max.mat')
pcd273_2hex_1xT = pcd273_2hex_T(:,1);
pcd273_2hex_10xT = pcd273_2hex_T(:,3);
pcd273_EB_1xT = pcd273_EB_T(:,1);
pcd273_EB_10xT = pcd273_EB_T(:,3);
pcd273_MV_1xT = pcd273_MV_T(:,1);
pcd273_MV_10xT = pcd273_MV_T(:,3);

pcd273_1xT = vertcat(pcd273_2hex_1xT,pcd273_EB_1xT,pcd273_MV_1xT);
pcd273_10xT = vertcat(pcd273_2hex_10xT,pcd273_EB_10xT,pcd273_MV_10xT);

pcd273_2hex_1xT1_max = pcd273_2hex_T1_max(:,1);
pcd273_2hex_10xT1_max = pcd273_2hex_T1_max(:,3);
pcd273_EB_1xT1_max = pcd273_EB_T1_max(:,1);
pcd273_EB_10xT1_max = pcd273_EB_T1_max(:,3);
pcd273_MV_1xT1_max = pcd273_MV_T1_max(:,1);
pcd273_MV_10xT1_max = pcd273_MV_T1_max(:,3);

pcd273_1xTmax = vertcat(pcd273_2hex_1xT1_max,pcd273_EB_1xT1_max,pcd273_MV_1xT1_max);
pcd273_10xTmax = vertcat(pcd273_2hex_10xT1_max,pcd273_EB_10xT1_max,pcd273_MV_10xT1_max);

pcd300_EB_1xT = pcd300_EB_T(:,1);
pcd300_ET_1xT = pcd300_ET_T(:,1);
pcd300_MB_1xT = pcd300_MB_T(:,1);
pcd300_MV_1xT = pcd300_MV_T(:,1);
pcd300_EB_10xT = pcd300_EB_T(:,3);
pcd300_ET_10xT = pcd300_ET_T(:,3);
pcd300_MB_10xT = pcd300_MB_T(:,3);
pcd300_MV_10xT = pcd300_MV_T(:,3);

pcd300_1xT = vertcat(pcd300_EB_1xT,pcd300_ET_1xT,pcd300_MB_1xT,pcd300_MV_1xT);
pcd300_10xT = vertcat(pcd300_EB_10xT,pcd300_ET_10xT,pcd300_MB_10xT,pcd300_MV_10xT);

pcd300_EB_1xT1_max = pcd300_EB_T1_max(:,1);
pcd300_ET_1xT1_max = pcd300_ET_T1_max(:,1);
pcd300_MB_1xT1_max = pcd300_MB_T1_max(:,1);
pcd300_MV_1xT1_max = pcd300_MV_T1_max(:,1);
pcd300_EB_10xT1_max = pcd300_EB_T1_max(:,3);
pcd300_ET_10xT1_max = pcd300_ET_T1_max(:,3);
pcd300_MB_10xT1_max = pcd300_MB_T1_max(:,3);
pcd300_MV_10xT1_max = pcd300_MV_T1_max(:,3);

pcd300_1xTmax = vertcat(pcd300_EB_1xT1_max,pcd300_ET_1xT1_max,pcd300_MB_1xT1_max,pcd300_MV_1xT1_max);
pcd300_10xTmax = vertcat(pcd300_EB_10xT1_max,pcd300_ET_10xT1_max,pcd300_MB_10xT1_max,pcd300_MV_10xT1_max);

pcd310_2hex_1xT = pcd310_2hex_T(:,1);
pcd310_ET_1xT = pcd310_ET_T(:,1);
pcd310_MB_1xT = pcd310_MB_T(:,1);
pcd310_MV_1xT = pcd310_MV_T(:,1);
pcd310_MV_10xT = pcd310_MV_T(:,3);
pcd310_MB_10xT = pcd310_MB_T(:,3);
pcd310_ET_10xT = pcd310_ET_T(:,3);
pcd310_2hex_10xT = pcd310_2hex_T(:,3);

pcd310_1xT = vertcat(pcd310_2hex_1xT,pcd310_ET_1xT,pcd310_MB_1xT,pcd310_MV_1xT);
pcd310_10xT = vertcat(pcd310_2hex_10xT,pcd310_ET_10xT,pcd310_MB_10xT,pcd310_MV_10xT);

pcd310_2hex_1xT1_max = pcd310_2hex_T1_max(:,1);
pcd310_ET_1xT1_max = pcd310_ET_T1_max(:,1);
pcd310_MB_1xT1_max = pcd310_MB_T1_max(:,1);
pcd310_MV_1xT1_max = pcd310_MV_T1_max(:,1);
pcd310_MV_10xT1_max = pcd310_MV_T1_max(:,3);
pcd310_MB_10xT1_max = pcd310_MB_T1_max(:,3);
pcd310_ET_10xT1_max = pcd310_ET_T1_max(:,3);
pcd310_2hex_10xT1_max = pcd310_2hex_T1_max(:,3);

pcd310_1xTmax = vertcat(pcd310_2hex_1xT1_max,pcd310_ET_1xT1_max,pcd310_MB_1xT1_max,pcd310_MV_1xT1_max);
pcd310_10xTmax = vertcat(pcd310_2hex_10xT1_max,pcd310_ET_10xT1_max,pcd310_MB_10xT1_max,pcd310_MV_10xT1_max);

pcd314_2hex_1xT = pcd314_2hex_T(:,1);
pcd314_ET_1xT = pcd314_ET_T(:,1);
pcd314_MV_1xT = pcd314_MV_T(:,1);
pcd314_MV_10xT = pcd314_MV_T(:,3);
pcd314_ET_10xT = pcd314_ET_T(:,3);
pcd314_2hex_10xT = pcd314_2hex_T(:,3);

pcd314_1xT = vertcat(pcd314_2hex_1xT,pcd314_ET_1xT,pcd314_MV_1xT);
pcd314_10xT = vertcat(pcd314_2hex_10xT,pcd314_ET_10xT,pcd314_MV_10xT);

pcd314_2hex_1xT1_max = pcd314_2hex_T1_max(:,1);
pcd314_ET_1xT1_max = pcd314_ET_T1_max(:,1);
pcd314_MV_1xT1_max = pcd314_MV_T1_max(:,1);
pcd314_MV_10xT1_max = pcd314_MV_T1_max(:,3);
pcd314_ET_10xT1_max = pcd314_ET_T1_max(:,3);
pcd314_2hex_10xT1_max = pcd314_2hex_T1_max(:,3);

pcd314_1xTmax = vertcat(pcd314_2hex_1xT1_max,pcd314_ET_1xT1_max,pcd314_MV_1xT1_max);
pcd314_10xTmax = vertcat(pcd314_2hex_10xT1_max,pcd314_ET_10xT1_max,pcd314_MV_10xT1_max);

allexp_1xT = vertcat(pcd273_1xT,pcd300_1xT,pcd310_1xT,pcd314_1xT);
allexp_10xT = vertcat(pcd273_10xT,pcd300_10xT,pcd310_10xT,pcd314_10xT);
allt2t1 = vertcat(allexp_1xT,allexp_10xT);
allt1max = vertcat(pcd273_1xTmax,pcd300_1xTmax,pcd310_1xTmax,pcd314_1xTmax);
allt10max = vertcat(pcd273_10xTmax,pcd300_10xTmax,pcd310_10xTmax,pcd314_10xTmax);
allt1and10max = vertcat(allt1max,allt10max);

figure;scatter(allt1max,allt10max);xlabel('{\Delta}F/F T1 1x');ylabel('{\Delta}F/F T1 10x');
figure;scatter(allt1and10max,allt2t1);xlabel('{\Delta}F/F T1 all conc');ylabel('T2-T1 all conc');

%% now to plot coherence?
clearvars
load('pcd_conc_results_t2-t1_t1coh.mat');
pcd273_2hex_1xT = pcd273_2hex_T(:,1);
pcd273_2hex_10xT = pcd273_2hex_T(:,3);
pcd273_EB_1xT = pcd273_EB_T(:,1);
pcd273_EB_10xT = pcd273_EB_T(:,3);
pcd273_MV_1xT = pcd273_MV_T(:,1);
pcd273_MV_10xT = pcd273_MV_T(:,3);

pcd273_1xT = vertcat(pcd273_2hex_1xT,pcd273_EB_1xT,pcd273_MV_1xT);
pcd273_10xT = vertcat(pcd273_2hex_10xT,pcd273_EB_10xT,pcd273_MV_10xT);

pcd273_2hex_1xT1_max = pcd273_2hex_T1_max(:,1);
pcd273_2hex_10xT1_max = pcd273_2hex_T1_max(:,3);
pcd273_EB_1xT1_max = pcd273_EB_T1_max(:,1);
pcd273_EB_10xT1_max = pcd273_EB_T1_max(:,3);
pcd273_MV_1xT1_max = pcd273_MV_T1_max(:,1);
pcd273_MV_10xT1_max = pcd273_MV_T1_max(:,3);

pcd273_1xTmax = vertcat(pcd273_2hex_1xT1_max,pcd273_EB_1xT1_max,pcd273_MV_1xT1_max);
pcd273_10xTmax = vertcat(pcd273_2hex_10xT1_max,pcd273_EB_10xT1_max,pcd273_MV_10xT1_max);

pcd273_2hex_1xC = pcd273_2hex_C(:,1);
pcd273_2hex_10xC = pcd273_2hex_C(:,3);
pcd273_EB_1xC = pcd273_EB_C(:,1);
pcd273_EB_10xC = pcd273_EB_C(:,3);
pcd273_MV_1xC = pcd273_MV_C(:,1);
pcd273_MV_10xC = pcd273_MV_C(:,3);

pcd273_1xC = vertcat(pcd273_2hex_1xC,pcd273_EB_1xC,pcd273_MV_1xC);
pcd273_10xC = vertcat(pcd273_2hex_10xC,pcd273_EB_10xC,pcd273_MV_10xC);

pcd300_EB_1xT = pcd300_EB_T(:,1);
pcd300_ET_1xT = pcd300_ET_T(:,1);
pcd300_MB_1xT = pcd300_MB_T(:,1);
pcd300_MV_1xT = pcd300_MV_T(:,1);
pcd300_EB_10xT = pcd300_EB_T(:,3);
pcd300_ET_10xT = pcd300_ET_T(:,3);
pcd300_MB_10xT = pcd300_MB_T(:,3);
pcd300_MV_10xT = pcd300_MV_T(:,3);

pcd300_1xT = vertcat(pcd300_EB_1xT,pcd300_ET_1xT,pcd300_MB_1xT,pcd300_MV_1xT);
pcd300_10xT = vertcat(pcd300_EB_10xT,pcd300_ET_10xT,pcd300_MB_10xT,pcd300_MV_10xT);

pcd300_EB_1xT1_max = pcd300_EB_T1_max(:,1);
pcd300_ET_1xT1_max = pcd300_ET_T1_max(:,1);
pcd300_MB_1xT1_max = pcd300_MB_T1_max(:,1);
pcd300_MV_1xT1_max = pcd300_MV_T1_max(:,1);
pcd300_EB_10xT1_max = pcd300_EB_T1_max(:,3);
pcd300_ET_10xT1_max = pcd300_ET_T1_max(:,3);
pcd300_MB_10xT1_max = pcd300_MB_T1_max(:,3);
pcd300_MV_10xT1_max = pcd300_MV_T1_max(:,3);

pcd300_1xTmax = vertcat(pcd300_EB_1xT1_max,pcd300_ET_1xT1_max,pcd300_MB_1xT1_max,pcd300_MV_1xT1_max);
pcd300_10xTmax = vertcat(pcd300_EB_10xT1_max,pcd300_ET_10xT1_max,pcd300_MB_10xT1_max,pcd300_MV_10xT1_max);

pcd300_EB_1xC = pcd300_EB_C(:,1);
pcd300_ET_1xC = pcd300_ET_C(:,1);
pcd300_MB_1xC = pcd300_MB_C(:,1);
pcd300_MV_1xC = pcd300_MV_C(:,1);
pcd300_EB_10xC = pcd300_EB_C(:,3);
pcd300_ET_10xC = pcd300_ET_C(:,3);
pcd300_MB_10xC = pcd300_MB_C(:,3);
pcd300_MV_10xC = pcd300_MV_C(:,3);

pcd300_1xC = vertcat(pcd300_EB_1xC,pcd300_ET_1xC,pcd300_MB_1xC,pcd300_MV_1xC);
pcd300_10xC = vertcat(pcd300_EB_10xC,pcd300_ET_10xC,pcd300_MB_10xC,pcd300_MV_10xC);

pcd310_2hex_1xT = pcd310_2hex_T(:,1);
pcd310_ET_1xT = pcd310_ET_T(:,1);
pcd310_MB_1xT = pcd310_MB_T(:,1);
pcd310_MV_1xT = pcd310_MV_T(:,1);
pcd310_MV_10xT = pcd310_MV_T(:,3);
pcd310_MB_10xT = pcd310_MB_T(:,3);
pcd310_ET_10xT = pcd310_ET_T(:,3);
pcd310_2hex_10xT = pcd310_2hex_T(:,3);

pcd310_1xT = vertcat(pcd310_2hex_1xT,pcd310_ET_1xT,pcd310_MB_1xT,pcd310_MV_1xT);
pcd310_10xT = vertcat(pcd310_2hex_10xT,pcd310_ET_10xT,pcd310_MB_10xT,pcd310_MV_10xT);

pcd310_2hex_1xT1_max = pcd310_2hex_T1_max(:,1);
pcd310_ET_1xT1_max = pcd310_ET_T1_max(:,1);
pcd310_MB_1xT1_max = pcd310_MB_T1_max(:,1);
pcd310_MV_1xT1_max = pcd310_MV_T1_max(:,1);
pcd310_MV_10xT1_max = pcd310_MV_T1_max(:,3);
pcd310_MB_10xT1_max = pcd310_MB_T1_max(:,3);
pcd310_ET_10xT1_max = pcd310_ET_T1_max(:,3);
pcd310_2hex_10xT1_max = pcd310_2hex_T1_max(:,3);

pcd310_1xTmax = vertcat(pcd310_2hex_1xT1_max,pcd310_ET_1xT1_max,pcd310_MB_1xT1_max,pcd310_MV_1xT1_max);
pcd310_10xTmax = vertcat(pcd310_2hex_10xT1_max,pcd310_ET_10xT1_max,pcd310_MB_10xT1_max,pcd310_MV_10xT1_max);

pcd310_2hex_1xC = pcd310_2hex_C(:,1);
pcd310_ET_1xC = pcd310_ET_C(:,1);
pcd310_MB_1xC = pcd310_MB_C(:,1);
pcd310_MV_1xC = pcd310_MV_C(:,1);
pcd310_MV_10xC = pcd310_MV_C(:,3);
pcd310_MB_10xC = pcd310_MB_C(:,3);
pcd310_ET_10xC = pcd310_ET_C(:,3);
pcd310_2hex_10xC = pcd310_2hex_C(:,3);

pcd310_1xC = vertcat(pcd310_2hex_1xC,pcd310_ET_1xC,pcd310_MB_1xC,pcd310_MV_1xC);
pcd310_10xC = vertcat(pcd310_2hex_10xC,pcd310_ET_10xC,pcd310_MB_10xC,pcd310_MV_10xC);

pcd314_2hex_1xT = pcd314_2hex_T(:,1);
pcd314_ET_1xT = pcd314_ET_T(:,1);
pcd314_MV_1xT = pcd314_MV_T(:,1);
pcd314_MV_10xT = pcd314_MV_T(:,3);
pcd314_ET_10xT = pcd314_ET_T(:,3);
pcd314_2hex_10xT = pcd314_2hex_T(:,3);

pcd314_1xT = vertcat(pcd314_2hex_1xT,pcd314_ET_1xT,pcd314_MV_1xT);
pcd314_10xT = vertcat(pcd314_2hex_10xT,pcd314_ET_10xT,pcd314_MV_10xT);

pcd314_2hex_1xT1_max = pcd314_2hex_T1_max(:,1);
pcd314_ET_1xT1_max = pcd314_ET_T1_max(:,1);
pcd314_MV_1xT1_max = pcd314_MV_T1_max(:,1);
pcd314_MV_10xT1_max = pcd314_MV_T1_max(:,3);
pcd314_ET_10xT1_max = pcd314_ET_T1_max(:,3);
pcd314_2hex_10xT1_max = pcd314_2hex_T1_max(:,3);

pcd314_1xTmax = vertcat(pcd314_2hex_1xT1_max,pcd314_ET_1xT1_max,pcd314_MV_1xT1_max);
pcd314_10xTmax = vertcat(pcd314_2hex_10xT1_max,pcd314_ET_10xT1_max,pcd314_MV_10xT1_max);

pcd314_2hex_1xC = pcd314_2hex_C(:,1);
pcd314_ET_1xC = pcd314_ET_C(:,1);
pcd314_MV_1xC = pcd314_MV_C(:,1);
pcd314_MV_10xC = pcd314_MV_C(:,3);
pcd314_ET_10xC = pcd314_ET_C(:,3);
pcd314_2hex_10xC = pcd314_2hex_C(:,3);

pcd314_1xC = vertcat(pcd314_2hex_1xC,pcd314_ET_1xC,pcd314_MV_1xC);
pcd314_10xC = vertcat(pcd314_2hex_10xC,pcd314_ET_10xC,pcd314_MV_10xC);

allexp_1xT = vertcat(pcd273_1xT,pcd300_1xT,pcd310_1xT,pcd314_1xT);
allexp_10xT = vertcat(pcd273_10xT,pcd300_10xT,pcd310_10xT,pcd314_10xT);
allt2t1 = vertcat(allexp_1xT,allexp_10xT);
allt1max = vertcat(pcd273_1xTmax,pcd300_1xTmax,pcd310_1xTmax,pcd314_1xTmax);
allt10max = vertcat(pcd273_10xTmax,pcd300_10xTmax,pcd310_10xTmax,pcd314_10xTmax);
allt1and10max = vertcat(allt1max,allt10max);
allC_1x = vertcat(pcd273_1xC,pcd300_1xC,pcd310_1xC,pcd314_1xC);
allC_10x = vertcat(pcd273_10xC,pcd300_10xC,pcd310_10xC,pcd314_10xC);
allC1and10 = vertcat(allC_1x,allC_10x);

figure;scatter(allt1max,allt10max);xlabel('{\Delta}F/F T1 1x');ylabel('{\Delta}F/F T1 10x');
figure;scatter(allt1and10max,allt2t1);xlabel('{\Delta}F/F T1 all conc');ylabel('T2-T1 all conc');
figure;scatter(allC_1x,allC_10x);xlabel('Coherence @ 1x');ylabel('Coherence @ 10x')
figure;scatter(allC_1x,allC_10x);xlabel('Coherence @ 1x');ylabel('Coherence @ 10x');xlim([0 1]),ylim([0 1])

%% now plot coherence for individual experiments/odorants
figure;scatter(pcd273_1xC,pcd273_10xC, 'o', 'MarkerFaceColor', 'r');
hold on;scatter(pcd300_1xC,pcd300_10xC, 'o', 'MarkerFaceColor', 'k');
hold on;scatter(pcd310_1xC,pcd310_10xC, 'o', 'MarkerFaceColor', 'b');
hold on;scatter(pcd314_1xC,pcd314_10xC, 'o', 'MarkerFaceColor', 'g');xlabel('2Hz Coherence @ 1x');ylabel('2Hz Coherence @ 10x');legend('pcd273','pcd300','pcd310','pcd314');hold off;
hex_concantC = vertcat(pcd273_2hex_C,pcd310_2hex_C,pcd314_2hex_C);
EB_concantC = vertcat(pcd273_EB_C,pcd300_EB_C);
MV_concantC = vertcat(pcd273_MV_C,pcd300_MV_C,pcd310_MV_C,pcd314_MV_C);
ET_concantC = vertcat(pcd300_ET_C,pcd310_ET_C,pcd314_ET_C);
MB_concantC = vertcat(pcd300_MB_C,pcd310_MB_C);
%now for concentrations
EB_1xC = EB_concantC(:,1);
EB_3xC = EB_concantC(:,2);
EB_10xC = EB_concantC(:,3);

hex_1xC = hex_concantC(:,1);
hex_3xC = hex_concantC(:,2);
hex_10xC = hex_concantC(:,3);

ET_1xC = ET_concantC(:,1);
ET_3xC = ET_concantC(:,2);
ET_10xC = ET_concantC(:,3);

MV_1xC = MV_concantC(:,1);
MV_3xC = MV_concantC(:,2);
MV_10xC = MV_concantC(:,3);

MB_1xC = MB_concantC(:,1);
MB_3xC = MB_concantC(:,2);
MB_10xC = MB_concantC(:,3);

figure;scatter(EB_1xC,EB_10xC,'o','r');hold on;
scatter(ET_1xC,ET_10xC,'o','k');hold on;
scatter(hex_1xC,hex_10xC,'o','b');hold on;
scatter(MB_1xC,MB_10xC,'o','g');hold on;
scatter(MV_1xC,MV_10xC,'o','m');xlabel('2Hz Coherence @ 1x');ylabel('2Hz Coherence @ 10x');
legend('EB','ET','hex','MB','MV');hold off

%%now for all odorants and experiments, plot across 3 concentrations:
figure;
for i=1:length(EB_concantC)
    plot([1 2 3],[EB_1xC(i), EB_3xC(i), EB_10xC(i)],'--o');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant ethyl-butyrate G-O pairs to 2Hz sniffing')    

figure;
for i=1:length(hex_concantC)
    plot([1 2 3],[hex_1xC(i), hex_3xC(i), hex_10xC(i)],'--o');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant 2-hexanone G-O pairs to 2Hz sniffing')

figure;
for i=1:length(ET_concantC)
    plot([1 2 3],[ET_1xC(i), ET_3xC(i), ET_10xC(i)],'--o');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant ethyl tiglate G-O pairs to 2Hz sniffing')

figure;
for i=1:length(MV_concantC)
    plot([1 2 3],[MV_1xC(i), MV_3xC(i), MV_10xC(i)],'--o');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant methyl valerate G-O pairs to 2Hz sniffing')

figure;
for i=1:length(MB_concantC)
    plot([1 2 3],[MB_1xC(i), MB_3xC(i), MB_10xC(i)],'--o');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant methyl benzoate G-O pairs to 2Hz sniffing')

odorconcant = vertcat(EB_concantC,ET_concantC,hex_concantC,MV_concantC,MB_concantC);
% graph all responses across concentrations
figure;for i=1:length(EB_concantC)
    plot([1 2 3],[EB_1xC(i), EB_3xC(i), EB_10xC(i)],'--o','Color','r','MarkerFaceColor','r');hold on;
end
for i=1:length(ET_concantC)
    plot([1 2 3],[ET_1xC(i), ET_3xC(i), ET_10xC(i)],'--o','Color','k','MarkerFaceColor','k');hold on;
end
for i=1:length(hex_concantC)
    plot([1 2 3],[hex_1xC(i), hex_3xC(i), hex_10xC(i)],'--o','Color','b','MarkerFaceColor','b');hold on;
end
for i=1:length(MB_concantC)
    plot([1 2 3],[MB_1xC(i), MB_3xC(i), MB_10xC(i)],'--o','Color','g','MarkerFaceColor','g');hold on;
end
for i=1:length(MV_concantC)
    plot([1 2 3],[MV_1xC(i), MV_3xC(i), MV_10xC(i)],'--o','Color','m','MarkerFaceColor','m');hold on;
end
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant methyl benzoate G-O pairs to 2Hz sniffing')

meanEBconcant = mean(EB_concantC);
meanETconcant = mean(ET_concantC);
meanhexconcant = mean(hex_concantC);
meanMBconcant = mean(MB_concantC);
meanMVconcant = mean(MV_concantC);

figure;plot([1 2 3],[meanEBconcant(:,1), meanEBconcant(:,2), meanEBconcant(:,3)],'-o','Color','r','MarkerFaceColor','r');hold on;
plot([1 2 3],[meanETconcant(:,1), meanETconcant(:,2), meanETconcant(:,3)],'-o','Color','k','MarkerFaceColor','k');hold on;
plot([1 2 3],[meanhexconcant(:,1), meanhexconcant(:,2), meanhexconcant(:,3)],'-o','Color','b','MarkerFaceColor','b');hold on;
plot([1 2 3],[meanMBconcant(:,1), meanMBconcant(:,2), meanMBconcant(:,3)],'-o','Color','g','MarkerFaceColor','g');hold on;
plot([1 2 3],[meanMVconcant(:,1), meanMVconcant(:,2), meanMVconcant(:,3)],'-o','Color','m','MarkerFaceColor','m');hold on;
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Mean Coherence to 2Hz Sniff');xlabel('Concentration');title('Coherence of Signficant methyl benzoate G-O pairs to 2Hz sniffing')
legend('EB','ET','hex','MB','MV');hold off

figure
subplot(5,1,1)
boxplot([EB_1xC, EB_3xC, EB_10xC],'PlotStyle','compact','Colors','r');
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');
subplot(5,1,2)
boxplot([ET_1xC,ET_3xC,ET_10xC],'PlotStyle','compact','Colors','k');
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');
subplot(5,1,3)
boxplot([hex_1xC,hex_3xC,hex_10xC],'PlotStyle','compact','Colors','b');
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');
subplot(5,1,4)
boxplot([MB_1xC,MB_3xC,MB_10xC],'PlotStyle','compact','Colors','g');
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');
subplot(5,1,5)
boxplot([MV_1xC,MV_3xC,MV_10xC],'PlotStyle','compact','Colors','m');
gca;xticks([1 2 3]);xticklabels({'1x','3x','10x'});ylabel('Coherence to 2Hz Sniff');xlabel('Concentration');

%% now to the text calculations in the paper
[h,p,ci,stats] = ttest(EB_concantC(:,1),EB_concantC(:,3))
[h,p,ci,stats] = ttest(MV_concantC(:,1),MV_concantC(:,3))
[h,p,ci,stats] = ttest(hex_concantC(:,1),hex_concantC(:,3))
[h,p,ci,stats] = ttest(ET_concantC(:,1),ET_concantC(:,3))% from 1 to 10x
[h,p,ci,stats] = ttest(ET_concantC(:,1),ET_concantC(:,2))%from 1 to 3x
mean_ET1x =mean(ET_concantC(:,1));
mean_ET3x =mean(ET_concantC(:,2));
mean_ET10x =mean(ET_concantC(:,3));
SEM_ET1x = std(ET_concantC(:,1))/sqrt(length(ET_concantC(:,1)));
SEM_ET3x = std(ET_concantC(:,2))/sqrt(length(ET_concantC(:,2)));
SEM_ET10x = std(ET_concantC(:,3))/sqrt(length(ET_concantC(:,3)));

figure;cdfplot(allexp_1xT)
hold on; cdfplot(allexp_10xT)
[h, p, k]=kstest2(allexp_1xT,allexp_10xT) % may not include correct sample size

figure;cdfplot(allt1max);
hold 'on';cdfplot(allt10max);
[h, p, k]=kstest2(allt1max,allt10max)

figure;cdfplot(allC_1x)
hold on;cdfplot(allC_10x)
[h, p, k]=kstest2(allC_1x,allC_10x)


