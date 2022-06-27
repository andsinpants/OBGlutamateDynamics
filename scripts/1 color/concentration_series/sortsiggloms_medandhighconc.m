%%load DON'T CLEAR!
%% set indexes of signficiant responses at highest concentration
goodspans_exc1x = goodspans_exc1x';
goodspans_sup1x= goodspans_sup1x';
sigexcroi_num1x = goodspans_exc1x(:,1:end);
sigsuproi_num1x = goodspans_sup1x(:,1:end);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now the same but for low concentrations
%% load odor file
disp('Select Medium Conc. Odor Plotdata File')
uiopen('load')

for i = 1:length(myplotdata.file.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.file.roi(i);
end

for j = 1:length(myplot_r) 
allplot_odor3x(j,:) = myplot_r(j).odor.avgtrial.series';
end

myplot_times = myplot_r(1).odor.avgtrial.time;
% allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented
%% chose same ROIs as in medium concentration series
for i=1
    if isempty(sigexcroi_num1x)
        break
    end
excrois_3x = allplot_odor3x(sigexcroi_num1x,:);
end

for i=1
    if isempty(sigsuproi_num1x)
        break
    end
suprois_3x = allplot_odor3x(sigsuproi_num1x,:);
end
%% smooth responses at medium concentration
excrois_3xsmth= [];
for i=1
if ~exist('excrois_3x','var') || isempty(excrois_3x)
    break
end
for j=1:length(excrois_3x(:,1))
    ys=excrois_3x(j,:);
    smoothb=fastsmooth(ys,15,3,1);
    excrois_3xsmth=[excrois_3xsmth;smoothb];
end
end

suprois_3xsmth= [];
for i=1
if ~exist('suprois_3x','var') || isempty(suprois_3x)
    break
end
for j=1:length(suprois_3x(:,1))
    ys=suprois_3x(j,:);
    smoothb=fastsmooth(ys,15,3,1);
    suprois_3xsmth=[suprois_3xsmth;smoothb];
end
end


%% norm -1 to 1 for lowest concentration
for i=1
    if ~exist('excrois_3xsmth','var') || isempty(excrois_3xsmth)
        break
    end
for j=1:size(excrois_3xsmth,1) %% in progress
    normexc3x(j,:) = normalised_diff(excrois_3xsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end


normsup3x=[];
for i=1
    if ~exist('suprois_3xsmth','var') || isempty(suprois_3xsmth)
        break
    end
for j=1:size(suprois_3xsmth,1) %% in progress
    normsup3x(j,:) = normalised_diff(suprois_3xsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normall3x = vertcat(normexc3x,normsup3x);

allplot_sigodorexc3x = [];
allplot_sigodorsup3x = [];
for i=1:length(goodspans_exc1x)
allplot_sigodorexc3x(i,:) = allplot_odor3x(goodspans_exc1x(i),:);
end
for i=1:length(goodspans_sup1x)
allplot_sigodorsup3x(i,:) = allplot_odor3x(goodspans_sup1x(i),:);
end
allplot_sigodor3x = vertcat(allplot_sigodorexc3x,allplot_sigodorsup3x);
%%%%%%%%%%%%%%%%%%%%%%
%% create imagesc fig
figure;imagesc(normexc3x);colorbar;title('Significantly Excited ROIs Normalized (3x)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
figure;imagesc(normsup3x);colorbar;title('Significantly Suppressed ROIs Normalized  (3x)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
figure;imagesc(normall3x);colorbar;title('Significantly Excited and Suppressed ROIs Normalized (3x)');
colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});


%% plot sorted, sig rois
figure;plot(normexc3x','r');hold on;plot(normsup3x','b'),title('Significantly Excited and Suppressed ROIs Normalized and Sorted (Traces) (3x)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now the same but for medium concentrations
%% load odor file
disp('Select Highest Conc. Odor Plotdata File')
uiopen('load')

for i = 1:length(myplotdata.file.roi) %move from struct into matrix
    myplot_r(i) = myplotdata.file.roi(i);
end

for j = 1:length(myplot_r) 
allplot_odor10x(j,:) = myplot_r(j).odor.avgtrial.series';
end

myplot_times = myplot_r(1).odor.avgtrial.time;
% allplot_odor = allplot_odor';
timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

%% chose same ROIs as in lowest concentration series
for i=1
    if isempty(sigexcroi_num1x)
        break
    end
excrois_10x = allplot_odor10x(sigexcroi_num1x,:);
end

for i=1
    if isempty(sigsuproi_num1x)
        break
    end
suprois_10x = allplot_odor10x(sigsuproi_num1x,:);
end
%% smooth responses at highest concentration
excrois_10xsmth= [];
for i=1
if ~exist('excrois_10x','var') || isempty(excrois_10x)
    break
end
for j=1:length(excrois_10x(:,1))
    ys=excrois_10x(j,:);
    smoothb=fastsmooth(ys,15,3,1);
    excrois_10xsmth=[excrois_10xsmth;smoothb];
end
end

suprois_10xsmth= [];
for i=1
if ~exist('suprois_10x','var') || isempty(suprois_10x)
    break
end
for j=1:length(suprois_10x(:,1))
    ys=suprois_10x(j,:);
    smoothb=fastsmooth(ys,15,3,1);
    suprois_10xsmth=[suprois_10xsmth;smoothb];
end
end


%% norm -1 to 1 for highest concentration
for i=1
    if ~exist('excrois_10xsmth','var') || isempty(excrois_10xsmth)
        break
    end
for j=1:size(excrois_10xsmth,1) %% in progress
    normexc10x(j,:) = normalised_diff(excrois_10xsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end


normsup10x=[];
for i=1
    if ~exist('suprois_10xsmth','var') || isempty(suprois_10xsmth)
        break
    end
for j=1:size(suprois_10xsmth,1) %% in progress
    normsup10x(j,:) = normalised_diff(suprois_10xsmth(j,1:todoron:todoroff:end)); %normalize to odor only
end
end

normall10x = vertcat(normexc10x,normsup10x);

allplot_sigodorexc10x = [];
allplot_sigodorsup10x = [];
for i=1:length(goodspans_exc1x)
allplot_sigodorexc10x(i,:) = allplot_odor10x(goodspans_exc1x(i),:);
end
for i=1:length(goodspans_sup1x)
allplot_sigodorsup10x(i,:) = allplot_odor10x(goodspans_sup1x(i),:);
end
allplot_sigodor10x = vertcat(allplot_sigodorexc10x,allplot_sigodorsup10x);
%%%%%%%%%%%%%%%%%%%%%%
%% create imagesc fig
figure;imagesc(normexc10x);colorbar;title('Significantly Excited ROIs Normalized (10x)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
figure;imagesc(normsup10x);colorbar;title('Significantly Suppressed ROIs Normalized  (10x)');colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
figure;imagesc(normall10x);colorbar;title('Significantly Excited and Suppressed ROIs Normalized (10x)');
colormap(gca,nawhimar_auto);
c = colorbar;
c.Label.String = 'Norm. \DeltaF/F';
gca;xlabel({'Frames'});ylabel({'ROIs'});


%% plot sorted, sig rois
figure;plot(normexc10x','r');hold on;plot(normsup10x','b'),title('Significantly Excited and Suppressed ROIs Normalized and Sorted (Traces) (10x)')
%% save variables
prompt = ('Do you want to save? Y/N?')
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if ~exist('suprois_3x','var') || isempty(suprois_3x)
        save ('excandsuproiunsorted_3x.mat', 'excrois_3xsmth','allplot_sigodor3x')
     elseif tf==1 && ~isempty(suprois_3x)
        save ('excandsuproisunsorted_3x.mat', 'excrois_3xsmth', 'suprois_3xsmth','allplot_sigodor3x')
    end
end

for i=1
    if ~exist('suprois_10x','var') || isempty(suprois_10x)
        save ('excandsuproiunsorted_10x.mat', 'excrois_10xsmth','allplot_sigodor10x')
     elseif tf==1 && ~isempty(suprois_10x)
        save ('excandsuproisunsorted_10x.mat', 'excrois_10xsmth', 'suprois_10xsmth','allplot_sigodor10x')
    end
end

for i=1
    if ~exist('suprois_3x','var') || isempty(suprois_3x)
        save ('excandsuprois_sortandnorm3x.mat', 'normexc3x')
     elseif tf==1 && ~isempty(suprois_3x)
        save ('excandsuprois_sortandnorm3x.mat', 'normexc3x','normsup3x')
    end
end

for i=1
    if ~exist('suprois_10x','var') || isempty(suprois_10x)
        save ('excandsuprois_sortandnorm10x.mat', 'normexc10x')
     elseif tf==1 && ~isempty(suprois_10x)
        save ('excandsuprois_sortandnorm10x.mat', 'normexc10x','normsup10x')
    end
end


%% close all figures
prompt = ('Do you want to close all figures? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf == 1
    close all
else
end
%clearvars -except excrois10x suprois10x 