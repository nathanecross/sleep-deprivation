Behave = readtable('../../data/cog/SDEP-ALL_2019-04-15.xlsx');
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

% Import sleep stats
sleepstats = readtable('/Volumes/NATHAN/SDEP/EEG_nap/derivatives/sleepstats/sleepstats.csv');
sleepstats(21:end,:) = [];

[r,p] = corrcoef(sleepstats.N3_pc,scorey(:,1),'rows','complete')

scatter(sleepstats.SL_min,scorey(:,1))



%% Slow wave power

% /Volumes/NATHAN/SDEP/EEG_tasks

load('ica/psd_cat.mat')
chans = table2cell(readtable('channels.csv'));
mask = ~startsWith(chans(:,2),'E');

coi = psd(:,:,mask,:);
coi_ave = squeeze(mean(coi,3));

% Remove outliers
coi_ave(coi_ave(:,1,2)>1000,:,:) = NaN;
coi_ave(coi_ave(:,1,4)>1000,:,:) = NaN;


for f = 1:10
    
    NN = squeeze(coi_ave(:,1,f));
    SD = squeeze(coi_ave(:,2,f));
    delta_f(:,f) = SD-NN;
end

delta_f(20,:) = NaN;

for f = 1:10
    [r,p] = corrcoef(scorey(:,1),delta_f(:,f),'rows','complete');
    cor(f) = r(2);
    pf(f) = p(2);
end

