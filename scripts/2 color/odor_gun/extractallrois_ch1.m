clc
clear all
close all
%% clear console and vars
clc
clearvars
%% load odor file 1Hz filtered
% disp('Select Avg Ch1 Odor Gun File (1Hz filtered)')
% uiopen('load')
% 
% for i = 1:length(myplotdata.avgfile.roi) %move from struct into matrix
%     myplot_r1hz(i) = myplotdata.avgfile.roi(i);
% end
% 
% myplot_onum = [myplot_r1hz(1).odor.number];
% myplot_times = [myplot_r1hz(1).odor(1).avgtrial.time];
% 
% disp(myplot_onum)
% disp(1:12)
% 
% for i=1:length(myplot_onum)
%     OdorNumID(i) = myplot_onum(i);
% end
% OdorNumID = OdorNumID';
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat1=repelem(OdorNumID(1),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat2=repelem(OdorNumID(2),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat3=repelem(OdorNumID(3),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat4=repelem(OdorNumID(4),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat5=repelem(OdorNumID(5),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat6=repelem(OdorNumID(6),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat7=repelem(OdorNumID(7),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat8=repelem(OdorNumID(8),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat9=repelem(OdorNumID(9),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat10=repelem(OdorNumID(10),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat11=repelem(OdorNumID(11),r)';
% end
% 
% for r = length(myplot_r1hz)
%     OdorNumID_mat12=repelem(OdorNumID(12),r)';
% end
% 
% OdorNumID = vertcat(OdorNumID_mat1,OdorNumID_mat2,OdorNumID_mat3,OdorNumID_mat4,OdorNumID_mat5,OdorNumID_mat6,OdorNumID_mat7,OdorNumID_mat8,OdorNumID_mat9,OdorNumID_mat10,OdorNumID_mat11,OdorNumID_mat12);
% totplot_odor = [];
% allplot_odor = [];
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor1(r,:) = myplot_r1hz(r).odor(1).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor2(r,:) = myplot_r1hz(r).odor(2).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor3(r,:) = myplot_r1hz(r).odor(3).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor4(r,:) = myplot_r1hz(r).odor(4).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor5(r,:) = myplot_r1hz(r).odor(5).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor6(r,:) = myplot_r1hz(r).odor(6).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor7(r,:) = myplot_r1hz(r).odor(7).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor8(r,:) = myplot_r1hz(r).odor(8).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor9(r,:) = myplot_r1hz(r).odor(9).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor10(r,:) = myplot_r1hz(r).odor(10).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor11(r,:) = myplot_r1hz(r).odor(11).avgtrial.series';
% end
% 
% for r = 1:length(myplot_r1hz)
%     totplot_odor12(r,:) = myplot_r1hz(r).odor(12).avgtrial.series';
% end
% 
% allplot_odor1hz = vertcat(totplot_odor1,totplot_odor2,totplot_odor3,totplot_odor4,totplot_odor5,totplot_odor6,totplot_odor7,totplot_odor8,totplot_odor9,totplot_odor10,totplot_odor11,totplot_odor12);
% totplot_odorch1hz = horzcat(OdorNumID,allplot_odor1hz);
% 
% timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented
% 
% figure;imagesc(allplot_odor1hz);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
% xlabel('Frames')
% ylabel('ROI Number')%df/f allplot of all ROIs
% 
% %clearvars -except allplot_odor timemat_odor myplot_onum onum
% 
% figure;plot(allplot_odor1hz');
% % find times of stimulus onset and offset
% % prompt = 'Odor on?';
% % odoron = input(prompt);
% % prompt = 'Odor off?';
% % odoroff = input(prompt);
% % timemat_odoron = find(timemat_odor(1,:)>=odoron);
% % timemat_odoroff = find(timemat_odor(1,:)<=odoroff);
% % todoron = timemat_odoron(1);
% % todoroff = timemat_odoroff(end);
%% now extract unfiltered
%% load odor file unfiltered
disp('Select Avg Ch1 Odor Gun File (unfiltered (MW 2Hz filtered)')
uiopen('load')

for i = 1:length(myplotdata.avgfile.roi) %move from struct into matrix
    myplot_runfilt(i) = myplotdata.avgfile.roi(i);
end

myplot_onum = [myplot_runfilt(1).odor.number];
myplot_times = [myplot_runfilt(1).odor(1).avgtrial.time];

disp(myplot_onum)
disp(1:12)

for i=1:length(myplot_onum)
    OdorNumID(i) = myplot_onum(i);
end
OdorNumID = OdorNumID';

for r = length(myplot_runfilt)
    OdorNumID_mat1=repelem(OdorNumID(1),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat2=repelem(OdorNumID(2),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat3=repelem(OdorNumID(3),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat4=repelem(OdorNumID(4),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat5=repelem(OdorNumID(5),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat6=repelem(OdorNumID(6),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat7=repelem(OdorNumID(7),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat8=repelem(OdorNumID(8),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat9=repelem(OdorNumID(9),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat10=repelem(OdorNumID(10),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat11=repelem(OdorNumID(11),r)';
end

for r = length(myplot_runfilt)
    OdorNumID_mat12=repelem(OdorNumID(12),r)';
end

OdorNumID = vertcat(OdorNumID_mat1,OdorNumID_mat2,OdorNumID_mat3,OdorNumID_mat4,OdorNumID_mat5,OdorNumID_mat6,OdorNumID_mat7,OdorNumID_mat8,OdorNumID_mat9,OdorNumID_mat10,OdorNumID_mat11,OdorNumID_mat12);
totplot_odor = [];
allplot_odor = [];

for r = 1:length(myplot_runfilt)
    totplot_odor1(r,:) = myplot_runfilt(r).odor(1).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor2(r,:) = myplot_runfilt(r).odor(2).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor3(r,:) = myplot_runfilt(r).odor(3).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor4(r,:) = myplot_runfilt(r).odor(4).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor5(r,:) = myplot_runfilt(r).odor(5).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor6(r,:) = myplot_runfilt(r).odor(6).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor7(r,:) = myplot_runfilt(r).odor(7).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor8(r,:) = myplot_runfilt(r).odor(8).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor9(r,:) = myplot_runfilt(r).odor(9).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor10(r,:) = myplot_runfilt(r).odor(10).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor11(r,:) = myplot_runfilt(r).odor(11).avgtrial.series';
end

for r = 1:length(myplot_runfilt)
    totplot_odor12(r,:) = myplot_runfilt(r).odor(12).avgtrial.series';
end

allplot_odoruf = vertcat(totplot_odor1,totplot_odor2,totplot_odor3,totplot_odor4,totplot_odor5,totplot_odor6,totplot_odor7,totplot_odor8,totplot_odor9,totplot_odor10,totplot_odor11,totplot_odor12);
totplot_odorch1uf = horzcat(OdorNumID,allplot_odoruf);

timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

figure;imagesc(allplot_odoruf);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

