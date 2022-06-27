%%extract all responses to g
disp('Select Avg or Nonavg Odor Gun File')
[file,path] = uigetfile('*.mat')
load(file)

tfplotdata = regexp(file,'and');
tfpd = isempty(tfplotdata);

for i=1
    if tfpd==0
        for i = 1:length(myplotdata.avgfile.roi) %move from struct into matrix
        myplot_r(i) = myplotdata.avgfile.roi(i);
        end
    end
    if tfpd==1
        for i = 1:length(myplotdata.file.roi) %move from struct into matrix
        myplot_r(i) = myplotdata.file.roi(i);
        end
    end
end

myplot_onum = [myplot_r(1).odor.number];
myplot_times = [myplot_r(1).odor(1).avgtrial.time];

disp(myplot_onum)
disp(1:12)

for i=1:length(myplot_onum)
    OdorNumID(i) = myplot_onum(i);
end
OdorNumID = OdorNumID';

for r = length(myplot_r)
    OdorNumID_mat1=repelem(OdorNumID(1),r)';
end

for r = length(myplot_r)
    OdorNumID_mat2=repelem(OdorNumID(2),r)';
end

for r = length(myplot_r)
    OdorNumID_mat3=repelem(OdorNumID(3),r)';
end

for r = length(myplot_r)
    OdorNumID_mat4=repelem(OdorNumID(4),r)';
end

for r = length(myplot_r)
    OdorNumID_mat5=repelem(OdorNumID(5),r)';
end

for r = length(myplot_r)
    OdorNumID_mat6=repelem(OdorNumID(6),r)';
end

for r = length(myplot_r)
    OdorNumID_mat7=repelem(OdorNumID(7),r)';
end

for r = length(myplot_r)
    OdorNumID_mat8=repelem(OdorNumID(8),r)';
end

for r = length(myplot_r)
    OdorNumID_mat9=repelem(OdorNumID(9),r)';
end

for r = length(myplot_r)
    OdorNumID_mat10=repelem(OdorNumID(10),r)';
end

for r = length(myplot_r)
    OdorNumID_mat11=repelem(OdorNumID(11),r)';
end

for r = length(myplot_r)
    OdorNumID_mat12=repelem(OdorNumID(12),r)';
end

OdorNumID = vertcat(OdorNumID_mat1,OdorNumID_mat2,OdorNumID_mat3,OdorNumID_mat4,OdorNumID_mat5,OdorNumID_mat6,OdorNumID_mat7,OdorNumID_mat8,OdorNumID_mat9,OdorNumID_mat10,OdorNumID_mat11,OdorNumID_mat12);
totplot_odor = [];
allplot_odor = [];

for r = 1:length(myplot_r)
    totplot_odor1(r,:) = myplot_r(r).odor(1).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor2(r,:) = myplot_r(r).odor(2).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor3(r,:) = myplot_r(r).odor(3).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor4(r,:) = myplot_r(r).odor(4).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor5(r,:) = myplot_r(r).odor(5).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor6(r,:) = myplot_r(r).odor(6).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor7(r,:) = myplot_r(r).odor(7).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor8(r,:) = myplot_r(r).odor(8).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor9(r,:) = myplot_r(r).odor(9).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor10(r,:) = myplot_r(r).odor(10).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor11(r,:) = myplot_r(r).odor(11).avgtrial.series';
end

for r = 1:length(myplot_r)
    totplot_odor12(r,:) = myplot_r(r).odor(12).avgtrial.series';
end

allplot_odor = vertcat(totplot_odor1,totplot_odor2,totplot_odor3,totplot_odor4,totplot_odor5,totplot_odor6,totplot_odor7,totplot_odor8,totplot_odor9,totplot_odor10,totplot_odor11,totplot_odor12);
totplot_odor1ID = horzcat(OdorNumID_mat1,totplot_odor1);
totplot_odor2ID = horzcat(OdorNumID_mat2,totplot_odor2);
totplot_odor3ID = horzcat(OdorNumID_mat3,totplot_odor3);
totplot_odor4ID = horzcat(OdorNumID_mat4,totplot_odor4);
totplot_odor5ID = horzcat(OdorNumID_mat5,totplot_odor5);
totplot_odor6ID = horzcat(OdorNumID_mat6,totplot_odor6);
totplot_odor7ID = horzcat(OdorNumID_mat7,totplot_odor7);
totplot_odor8ID = horzcat(OdorNumID_mat8,totplot_odor8);
totplot_odor9ID = horzcat(OdorNumID_mat9,totplot_odor9);
totplot_odor10ID = horzcat(OdorNumID_mat10,totplot_odor10);
totplot_odor11ID = horzcat(OdorNumID_mat11,totplot_odor11);
totplot_odor12ID = horzcat(OdorNumID_mat12,totplot_odor12);
totplot_odor = horzcat(OdorNumID,allplot_odor);

timemat_odor = myplot_times; %use for now, should be the same timeseries across odorants presented

figure;imagesc(allplot_odor);colorbar;title('\DeltaF/F heatmap');rb = nawhimar_auto();colormap(rb)
xlabel('Frames')
ylabel('ROI Number')%df/f allplot of all ROIs

%% now save vars
prompt = ('Do you want to save vars? Y/N?');
savequestion = input(prompt,'s');
tfsave=strcmpi(savequestion,'Y');
if tfsave==1
    save( 'allodorsandgloms.mat', 'totplot_odor1ID', 'totplot_odor2ID', 'totplot_odor3ID', 'totplot_odor4ID', 'totplot_odor5ID', 'totplot_odor6ID', 'totplot_odor7ID', 'totplot_odor8ID','totplot_odor9ID','totplot_odor10ID','totplot_odor11ID','totplot_odor12ID');
else
end