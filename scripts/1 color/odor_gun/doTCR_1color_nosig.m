close all
clearvars

[atfile,atpath] = uigetfile('*.mat','Choose Channel #1 File(s) (2Hz filter?)','MultiSelect','off');
if ~atfile; return; end
ch1_at = load(fullfile(atpath,atfile)); clear atpath atfile;
% [supfile,suppath] = uigetfile('*.mat','Choose Channel #1 Suppressive File(s) (0.5Hz filter?)','MultiSelect','off');
% if ~supfile; return; end
% ch1_sup = load(fullfile(suppath,supfile)); clear suppath supfile;

numrois = length(ch1_at.myplotdata.avgfile.roi(:));
numodors = length(ch1_at.myplotdata.avgfile.roi(1).odor(:));
for i=1:numrois
for j=1:numodors
alltraces(i,j,:) = ch1_at.myplotdata.avgfile.roi(i).odor(j).avgtrial.series;
end
end

%% save matrix
prompt = ('Do you want to save traces.mat? Y/N?');
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
if tf == 1
        save ('alltraces.mat', 'alltraces')
end
