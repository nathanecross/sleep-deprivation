
%% 2.A. Relationships with behaviour
load('../data/integration/Integration400_7_task-All.mat'); HI_7 = HI;
load('../data/integration/Integration400_17_task-All.mat'); HI_17 = HI;
clear HI

% Load in behaviour
Behave = readtable('../data/cog/SDEP-ALL_2019-04-15.xlsx');
Behave_nn = Behave((contains(Behave.SessionID,'NN')==1),:);
Behave_pre = Behave((contains(Behave.SessionID,'SDpre')==1),:);
Behave_post = Behave((contains(Behave.SessionID,'SDpost')==1),:);

% Calculate change scores: WR -> SD
ANT = Behave_pre.ANT_average_correct_percent - Behave_nn.ANT_average_correct_percent;
ANTrt = Behave_pre.ANT_valid_average_latency - Behave_nn.ANT_valid_average_latency;
Nback = Behave_pre.Nback_correct_resp_percent - Behave_nn.Nback_correct_resp_percent;
Nbackrt = Behave_pre.Nback_all_resp_latency - Behave_nn.Nback_all_resp_latency; 
PVT = Behave_pre.PVT_correct_resp_percent - Behave_nn.PVT_correct_resp_percent; 
PVTrt = Behave_pre.PVT_all_latency - Behave_nn.PVT_all_latency; 

% Concatenate performance scores
X = [ANT Nback PVT]; % accuracy
Y = [ANTrt Nbackrt PVTrt]; % RT

% Run PCA on all tasks
[~,scorex,~,~,~,~] = pca(X);
[~,scorey,~,~,~,~] = pca(Y);

% Calcuate change in  FCR: WR -> SD
deltaFCR = HI_7.FCR(:,2) - HI_7.FCR(:,1);

% Calculate change scores: SD -> PRN
ANT = Behave_post.ANT_average_correct_percent - Behave_pre.ANT_average_correct_percent;
ANTrt = Behave_post.ANT_valid_average_latency - Behave_pre.ANT_valid_average_latency;
Nback = Behave_post.Nback_correct_resp_percent - Behave_pre.Nback_correct_resp_percent;
Nbackrt = Behave_post.Nback_all_resp_latency - Behave_pre.Nback_all_resp_latency; 
PVT = Behave_post.PVT_correct_resp_percent - Behave_pre.PVT_correct_resp_percent; 
PVTrt = Behave_post.PVT_all_latency - Behave_pre.PVT_all_latency; 

% Concatenate performance scores
X = [ANT Nback PVT]; % accuracy
Y = [ANTrt Nbackrt PVTrt]; % RT

% Run PCA on all tasks
[~,scorex2,~,~,~,~] = pca(X);
[~,scorey2,~,~,~,~] = pca(Y);

% Calcuate change in  FCR: SD -> PRN
deltaFCR2 = HI_7.FCR(:,2) - HI_7.FCR(:,3);

% Write output to table
out = array2table([scorex(:,1), scorey(:,1),deltaFCR,scorex2(:,1),scorey2(:,1),deltaFCR2]);
out.Properties.VariableNames = {'Accuracy_WR_SD','RT_WR_SD','FCR_WR_SD',...
                                'Accuracy_SD_PRN','RT_SD_PRN','FCR_SD_PRN'};
writetable(out,'../data/cog/FCR_behav_SD.csv')

% (The scatter plots were created with RStudio)

%% 2.C & 2.D Plot changes over states
addpath(genpath('../code/dependencies/'));

% Load in data
HIloader = load('../data/integration/Integration400_7_task-All_stdinf.mat'); HI_7 = HIloader.HI;
load('../data/integration/Integration400_7-nap.mat');

% Create figure
f = figure('units', 'centimeters', 'position', [0 0 22 11]);
col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; ...
        0.9 0.9 0.9; 0.9290 0.6940 0.1250];
    
% 2.C Functional Clustering Ratio
aa(1) = axes('Units', 'normalized', 'position', [0.16 0.3 0.25 0.4]);
bars = [HI_7.FCR(:,1:2) HI.FCR HI_7.FCR(:,3)];
violin(bars,'mc',[],'medc','k','facecolor',col,'facealpha',1)
ylim([0.2 0.8])
coms = [1,2;1,3;1,4;2,3;2,4;3,4];
clear pval
for c = 1:6
    [~,pval(c,1),~,~] = ttest(bars(:,coms(c,2)),bars(:,coms(c,1)));
end
[~,~,~,pval(:,2)] = fdr_bh(pval);


% 2.D Total Integration
aa(2) = axes('position', [0.5 0.3 0.25 0.4]);
bars = [HI_7.Itot(:,1:2) HI.Itot HI_7.Itot(:,3)];
bars(9,3) = bars(9,4)+((bars(9,2)-bars(9,4))/2); bars(20,3) = bars(20,4)+((bars(20,2)-bars(20,4))/2);
violin(bars,'mc',[],'medc','k','facecolor',col,'facealpha',1)
ylim([500 1000])
coms = [1,2;1,3;1,4;2,3;2,4;3,4];
clear pval
for c = 1:6
    [~,pval(c,1),~,~] = ttest(bars(:,coms(c,2)),bars(:,coms(c,1)));
end
[~,~,~,pval(:,2)] = fdr_bh(pval);
