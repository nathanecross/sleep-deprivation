% 
%  Script to run integration analyses on cortical functional data using  
%  hierarchical (Bayesian) group inference, using communites estimated using
%  the louvain community detection algorithm (taken from Conn Toolbox)
%
%  This script is written with the expectation that the BOLD timeseries
%  have been preprocessed and clustered into parcels or ROIs. Each ROI
%  is also required to be assigned to a higher-level/lower-order group or
%  "network" (ie. the hierarchical aspect of these integration analyses).

%% Setup
sub_dir = '../../data/Parcels_cat/';
fs = 'fsaverage5';                  %<-- (Choose number of vertices as desired)
nparc = '400';                      %<-- (Choose number of parcellations as desired)
NumberOfParcels = str2num(nparc);
nnetw = '7';                        %<-- (Choose number of networks as desired)
NumberOfNets = str2num(nnetw);
network_ids = 1:NumberOfNets;
load('../../labels/louv_M.mat'); %loads in louvain clusters mask
mask_rois = louv_M; 
load(sprintf('../../labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw))  %loads in parcel mask

% Map Schaeffer parcellations to the Yeo Networks     
Yeo_Clusters_ref = load(sprintf('../../labels/fsaverage5/YeoNetworks/1000subjects_clusters%s_ref.mat',nnetw));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

% Load in subject and study information
load('../../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = length(SubjectName);
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];

%% Run Hierachical Integration within the whole brain

sprintf("Running Hierachical Integration within the whole brain")
addpath(genpath('../dependencies/'));
for t = 1:3
    for z = 1:NumberOfSubjects
        load( [sub_dir sprintf('%s/%s/%s_task-All.mat',  ...
            (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})))] );
    fmri{z,1} = Parcels;
    end
    opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 800);
    sprintf("Running Session -> %s",(char(times(t,3))))
    HI.(char(times(t,3))) = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);
end

Itot = zeros(NumberOfSubjects,3); I_inter = zeros(NumberOfSubjects,3); I_intra = zeros(NumberOfSubjects,7,3);

for t=1:3
    for z = 1:NumberOfSubjects
        Itot(z,t) = HI.(char(times(t,3))){z, 1}.int_total.mean;
        I_inter(z,t) = HI.(char(times(t,3))){z, 1}.int_inter.mean;
        for n = 1:NumberOfNets
            I_intra(z,n,t) = HI.(char(times(t,3))){z, 1}.int_intra.mean(n,n);
        end
        Iws(z,t) = squeeze(sum(I_intra(z,:,t)));
        Ibs(z,t) = I_inter(z,t);
        FCR(z,t) = Iws(z,1)/Ibs(z,t);
    end
    I_intra_mean(:,t) = mean(squeeze(I_intra(:,:,t)),1)';    
end

I_intra_diff(:,1) = I_intra_mean(:,2) - I_intra_mean(:,1);
I_intra_diff(:,2) = I_intra_mean(:,3) - I_intra_mean(:,2);
I_intra_diff(:,3) = I_intra_mean(:,3) - I_intra_mean(:,1);

Yeo_labels = [lh_labels_Shf_Yeo; rh_labels_Shf_Yeo];
for n = 1:17
    for y = 1:20484
        if Yeo_labels(y,2) == n
            Yeo_labels(y,3) = I_intra_diff(n,1);
        end
    end
end

HI.Itot = Itot; HI.Iws = Iws; HI.Ibs = Ibs; HI.FCR = FCR;

% Stats
stats(1,1) = HI.(char(times(1,3))){21,1}.int_total.mean;
stats(1,2) = HI.(char(times(1,3))){21, 1}.int_total.var;
stats(1,3) = sum(HI.(char(times(1,3))){21,1}.int_total.samples>...
                   HI.(char(times(2,3))){21,1}.int_total.samples)...
                   /length(HI.(char(times(1,3))){21,1}.int_total.samples);
stats(1,4) = HI.(char(times(2,3))){21,1}.int_total.mean;
stats(1,5) = HI.(char(times(2,3))){21, 1}.int_total.var;
stats(1,6) = sum(HI.(char(times(3,3))){21,1}.int_total.samples>...
                   HI.(char(times(2,3))){21,1}.int_total.samples)...
                   /length(HI.(char(times(2,3))){21,1}.int_total.samples);
stats(1,7) = HI.(char(times(3,3))){21,1}.int_total.mean;
stats(1,8) = HI.(char(times(3,3))){21, 1}.int_total.var;  

stats = array2table(stats);
stats.Properties.VariableNames = {'WRmean','WRstd','sig','SDmean','SDstd','sig2','PRNmean','PRNstd'}


HI.stats = stats;

fprintf("Saving Data For Whole Cortex -> 7 Networks");
save(sprintf('../../data/Integration%s_7_task-All_louv.mat',nparc), 'HI');


%% Bayesian inference for FCR
load('../../data/Integration400_7_task-All_louv.mat')

% Load samples from Gibbs sampling for each measure
for t = 1:3
    for s = 1:1000
        Iws(s,t) = sum(diag(HI.(char(times(t,3))){21, 1}.int_intra.samples(:,:,s)));
    end
    Itot(:,t) = HI.(char(times(t,3))){21, 1}.int_total.samples;
    Ibs(:,t) = HI.(char(times(t,3))){21, 1}.int_inter.samples;
    FCR(:,t) = Iws(:,t)./Ibs(:,t);
end

% Create table & compare in a Bayesian framework
sig = NaN(6,2);
sig(1,1) = sum(squeeze(FCR(:,1))>squeeze(FCR(:,2)))/1000;
sig(1,2) = sum(squeeze(FCR(:,3))>squeeze(FCR(:,2)))/1000;
table = [mean(FCR(:,1)), std(FCR(:,1)), sig(1,1), mean(FCR(:,2)), std(FCR(:,2)),...
         sig(1,2),mean(FCR(:,3)), std(FCR(:,3))]

