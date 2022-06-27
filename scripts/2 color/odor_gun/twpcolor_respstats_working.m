close all; 
clear;
clc;

[exfile,expath] = uigetfile('*.mat','Choose Channel #1 Excitatory File(s) (2Hz filter?)','MultiSelect','off');
if ~exfile; return; end
ch1_ex = load(fullfile(expath,exfile)); clear expath exfile;
[supfile,suppath] = uigetfile('*.mat','Choose Channel #1 Suppressive File(s) (0.5Hz filter?)','MultiSelect','off');
if ~supfile; return; end
ch1_sup = load(fullfile(suppath,supfile)); clear suppath supfile;

bplotbyodor = 1; %1 means figure# same as odor#; 0 means figure# is roi#
ibase = 1:275; %preodor frames
iodor = 301:825; %postodor frames
minthresh = -7; maxthresh = 7;

odornums = [];
for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
    for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornums = union(odornums,ch1_ex.myplotdata.avgfile.roi(r).odor(o).number,'sorted');
    end
end
results = nan(length(ch1_ex.myplotdata.avgfile.roi),odornums(end),2);

figure(256); hold on; %odor number vs significant min/max values
title('Channel #1');
exbasemeanall = [];
tmpmaxmat = [];
for r = 1:length(ch1_ex.myplotdata.avgfile.roi)
    for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exbasemeanall = vertcat(exbasemeansub,exbasemeanall);
    end
    exbasestd = std(exbasemeanall);
     for o = 1:length(ch1_ex.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_ex.myplotdata.avgfile.roi(r).odor(o).number;
        exbasemean = mean(ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        exbasemeansub = ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-exbasemean;
        exzts = (ch1_ex.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-exbasemean)./exbasestd;
        tmpmax = prctile(exzts(ibase+25),95);
        tmpmaxmat(r,o) = tmpmax;
      end
end

figure;histogram(tmpmaxmat(1:42,:));

supbasemeanall = [];
tmpminmat = [];
for r = 1:length(ch1_sup.myplotdata.avgfile.roi)
    for o = 1:length(ch1_sup.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_sup.myplotdata.avgfile.roi(r).odor(o).number;
        supbasemean = mean(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supbasemeanall = vertcat(supbasemeansub,supbasemeanall);
    end
    supbasestd = std(supbasemeanall);
     for o = 1:length(ch1_sup.myplotdata.avgfile.roi(r).odor)
        odornum = ch1_sup.myplotdata.avgfile.roi(r).odor(o).number;
        supbasemean = mean(ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase));
        supbasemeansub = ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series(ibase)-supbasemean;
        supzts = (ch1_sup.myplotdata.avgfile.roi(r).odor(o).avgtrial.series-supbasemean)./supbasestd;
        tmpmin = prctile(supzts(ibase),10);
        tmpminmat(r,o) = tmpmin;
      end
end

figure;histogram(tmpminmat(1:42,:));
