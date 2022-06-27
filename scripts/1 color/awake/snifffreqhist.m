% clc
% clearvars
%% multiselect trials
[file,path] = uigetfile({'*.mat'},'Pick Files to Load','MultiSelect','on');

cd(path)
load(file)

prompt = 'Is this a Behavioral data file? Y/N?';
savequestion = input(prompt,'s');
tf=strcmpi(savequestion,'Y');
for i=1
    if tf==1
        for i=1:length(behaviordata.trials)
        tdatanames{i} = behaviordata.trials(i).name';
        end
   
        tdatanames = tdatanames';
%% get frequency of sniffing bouts
% convert to cell
        for i=1:length(behaviordata.trials)
            tdatadetsniffs = behaviordata.trials(i).det_sniffs.inhalations_freq;
            tdatasniff{i} = {tdatadetsniffs(1:end)};
        end

        for i=1:length(behaviordata.trials)
            tdatamat = horzcat(tdatasniff{1,1:i});
        %     tdatamat_all = horzcat(tdatamat{1,1:i});
        end

        for i=1:length(behaviordata.trials)
        tdatamat_all = horzcat(tdatamat{1,1:i});
        end

        %% plot frequency of sniffing bouts
        figure;histogram(tdatamat_all);xlabel('Sniffing Frequencies (Hz)');ylabel('Counts');
    end
    if tf==0
        for i=1:length(trialsdata.trials)
        tdatanames{i} = trialsdata.trials(i).name';
        end
   
        tdatanames = tdatanames';
%% get frequency of sniffing bouts
% convert to cell
        for i=1:length(trialsdata.trials)
            tdatadetsniffs = trialsdata.trials(i).det_sniffs.troughs_freq;
            tdatasniff{i} = {tdatadetsniffs(1:end)};
        end

        for i=1:length(trialsdata.trials)
            tdatamat = horzcat(tdatasniff{1,1:i});
        %     tdatamat_all = horzcat(tdatamat{1,1:i});
        end

        for i=1:length(trialsdata.trials)
        tdatamat_all = horzcat(tdatamat{1,1:i});
        end

        %% plot frequency of sniffing bouts
        figure;histogram(tdatamat_all);xlabel('Sniffing Frequencies (Hz)');ylabel('Counts');
    end
end