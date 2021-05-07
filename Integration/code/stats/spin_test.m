%% Spin test for significance
addpath(genpath('../dependencies'));




%% 7 Networks
keyword = 'task-All';
nparc = '400'; NumberOfParcels = str2num(nparc);
nnetw = '7'; NumberOfNetworks = str2num(nnetw);
load(sprintf('../../data/integration/Integration%s_%s_%s_stdinf.mat',nparc,nnetw,keyword)); HI_7 = HI;
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
load(sprintf('../../labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw))
Yeo_Clusters_ref = load(sprintf('../../labels/fsaverage5/YeoNetworks/1000subjects_clusters%s_ref.mat',nnetw)); % 17 Networks
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

% Calculate change in Integration per network
for n = 1:NumberOfNetworks
     [~,~,~,stats] = ttest(HI_7.Networks.Itot(n,:,2) - HI_7.Networks.Itot(n,:,1));
     I7_change(n) = stats.tstat;
end
for p = 1:NumberOfParcels
    I7_changevec(p) = I7_change(Yeo_Shf(p));   
end

% Project onto surface
SP = SurfStatAvSurf({'../../labels/fsaverage5/surf/lh.pial', ...
                     '../../labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = I7_changevec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);

% Run Permutations and save data
csvwrite(sprintf('../../data/integration/I%s_changeveclh.csv',nnetw),OnSurf(1:10242)')
csvwrite(sprintf('../../data/integration/I%s_changevecrh.csv',nnetw),OnSurf(10243:end)')
readleft = sprintf('../../data/integration/I%s_changeveclh.csv',nnetw);
readright = sprintf('../../data/integration/I%s_changevecrh.csv',nnetw);
SpinPermuFS(readleft,readright,100,sprintf('../../data/stats/I%s_change_perm.mat',nnetw))

%% 17 Networks
keyword = 'task-All';
nparc = '400'; NumberOfParcels = str2num(nparc);
nnetw = '17'; NumberOfNetworks = str2num(nnetw);
load(sprintf('../../data/integration/Integration%s_%s_%s_stdinf.mat',nparc,nnetw,keyword)); HI_17 = HI;
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
load(sprintf('../../labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw))
Yeo_Clusters_ref = load(sprintf('../../labels/fsaverage5/YeoNetworks/1000subjects_clusters%s_ref.mat',nnetw)); % 17 Networks
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

% Calculate change in Integration per network
for n = 1:NumberOfNetworks
    [~,~,~,stats] = ttest(HI_17.Networks.Itot(n,:,2) - HI_17.Networks.Itot(n,:,1));
    I17_change(n) = stats.tstat;
end
for p = 1:NumberOfParcels
    I17_changevec(p) = I17_change(Yeo_Shf(p));   
end

% Project onto surface
SP = SurfStatAvSurf({'../../labels/fsaverage5/surf/lh.pial', ...
                     '../../labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = I17_changevec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);

% Run Permutations and save data
csvwrite(sprintf('../../data/integration/I%s_changeveclh.csv',nnetw),OnSurf(1:10242)')
csvwrite(sprintf('../../data/integration/I%s_changevecrh.csv',nnetw),OnSurf(10243:end)')
readleft = sprintf('../../data/integration/I%s_changeveclh.csv',nnetw);
readright = sprintf('../../data/integration/I%s_changevecrh.csv',nnetw);
SpinPermuFS(readleft,readright,100,sprintf('../../data/stats/I%s_change_perm.mat',nnetw))
