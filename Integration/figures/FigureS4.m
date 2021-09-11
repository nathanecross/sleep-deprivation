%% FIGURE S7
%  Plot total, within-system and between-system integration across subjects

%% Setup
addpath(genpath('dependencies/'));
fs = 'fsaverage5'; %<- (Choose number of vertices as desired)
nparc = '400'; %<- (Choose number of parcellations as desired)
nnetw = '7'; %<- (Choose number of networks as desired)
NumberOfParcels = str2num(nparc);
NumberOfNets = str2num(nnetw);
network_ids = 1:NumberOfNets;
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in Yeo network mask
mask_rois = Yeo_Shf; 
load('SubjectName.mat');
NumberOfSubjects = length(SubjectName);
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];

load(sprintf('labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw)) 
Yeo_Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters007_ref.mat'); % 7 Networks
%Yeo_Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters017_ref.mat'); % 17 Networks

% Map Schaeffer parcellations to the Yeo Networks
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

%% Calculate Hierachical Integration within the whole brain during wake

load(sprintf('data/Parcels%s_%s_task-All.mat', string(NumberOfParcels),nnetw));

for t = 1:length(times(:,1))
    for z = 1:NumberOfSubjects
        for p = 1:NumberOfParcels
            mat(:,p) = Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfMeanclean.(char("p"+p))';
        end
    fmri{z,1} = mat;
    end
    opt_inference = struct('method', 'standard-inference', 'n_samplings', 1000, 'nu_max', 800);
    sprintf("Running Session -> %s",(char(times(t,3))))
    HI.(char(times(t,3))) = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);
end
save('data/Integration400_7_stdinf.mat', 'HI');

%% Calculate summaries of total, within and between integration during wake

load('data/Integration400_7_stdinf.mat'); 
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}]; 

for t = 1:3
    for z = 1:20
        Iwsi(z,t) = sum(diag(HI.(char(times(t,3))){z, 1}.int_intra.mean));
        Ibsi(z,t) = HI.(char(times(t,3))){z, 1}.int_inter.mean;
        Itoti(z,t) = HI.(char(times(t,3))){z, 1}.int_total.mean;
    end
end

%% Calculate Hierachical Integration within the whole brain during nap

load(sprintf('data/Parcels%s_17_nap-Allchunks.mat', string(NumberOfParcels)));
load('SubjectName_n18.mat');
NumberOfSubjects = length(SubjectName);

for z = 1:NumberOfSubjects
    mat = Parcels.nap.(char(SubjectName{z})).ShfNet';
    fmri{z,1} = mat;
end
opt_inference = struct('method', 'standard-inference', 'n_samplings', 1000, 'nu_max', 800);
HI.Nap = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);

save('data/Integration400_7-nap_stdinf.mat','HI');

%% Calculate summaries of total, within and between integration during nap

HIloader = load('data/Integration400_7-nap_stdinf.mat'); HI = HIloader.HI;

for z = 1:18
        Iwsin(z,1) = sum(diag(HI.Nap{z, 1}.int_intra.mean));
        Ibsin(z,1) = HI.Nap{z, 1}.int_inter.mean;
        Itotin(z,1) = HI.Nap{z, 1}.int_total.mean;
end

Itoti(:,4) = [Itotin(1:8);NaN;Itotin(9:end);NaN];
Iwsi(:,4) = [Iwsin(1:8);NaN;Iwsin(9:end);NaN];
Ibsi(:,4) = [Ibsin(1:8);NaN;Ibsin(9:end);NaN];

% Rearrange so that nap is 3rd
Itoti = [Itoti(:,1:2), Itoti(:,4), Itoti(:,3)];
Ibsi = [Ibsi(:,1:2), Ibsi(:,4), Ibsi(:,3)];
Iwsi = [Iwsi(:,1:2), Iwsi(:,4), Iwsi(:,3)];

%% Plot 

col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; ...
        0.9 0.9 0.9; 0.9290 0.6940 0.1250];
    
%Individual matrices
figure;
data = [Itoti, zeros(20,1),Iwsi,zeros(20,1),Ibsi];
color = [col;[0,0,0]; col;[0,0,0]; col];
violin(data,'mc',[],'medc','k','facecolor',color,'facealpha',1)
ylim([10 1100])
legend('off')    