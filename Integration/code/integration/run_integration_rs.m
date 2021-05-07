%  (Formatted for resting-state) 
%
%  Script to run integration analyses on cortical functional data using  
%  hierarchical (Bayesian) group inference.
%
%  This script is written with the expectation that the BOLD timeseries
%  have been preprocessed and clustered into parcels or ROIs. Each ROI
%  is also required to be assigned to a higher-level/lower-order group or
%  "network" (ie. the hierarchical aspect of these integration analyses).

% cd('/home/ncross/projects/def-ttdangvu/SDEP/') % (main directory containing scripts and data)



%% Setup
addpath(genpath('../dependencies/'));
addpath(genpath(pwd));
fs = 'fsaverage5'; %<-------------------------------------- (Choose number of vertices as desired)
nparc = '100'; %<------------------------------------------ (Choose number of parcellations as desired)
nnetw = '7'; %<------------------------------------------- (Choose number of networks as desired)  
NumberOfParcels = str2num(nparc);
NumberOfNets = str2num(nnetw);
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in Yeo network mask
network_ids = 1:NumberOfNets;
mask_rois = Yeo_Shf;   
load(sprintf('../../labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw)) 

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

% 1. Arrange data and calculate integration for subjects 
for t = 1:3
    % a. Organise BOLD timeseries data into matrix (time x parcels)
    for z = 1:NumberOfSubjects
        load(sprintf('../../data/Parcels/%s/%s/%s_task-CROSS.mat',  ...
            (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z}))));
    fmri{z,1} = double(Parcels(:,1:120)'); % put all subjects' matrices into cell array
    end
    % b. Specify options and calculate integration on data matrix
    opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 400);
    HI.(char(times(t,3))) = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);
end

% 2. Create summary table
stats(1,1) = HI.(char(times(1,3))){21,1}.int_total.mean;
stats(1,2) = HI.(char(times(1,3))){21, 1}.int_total.var;
stats(1,3) = sum(HI.(char(times(1,3))){21,1}.int_total.samples>...
                   HI.(char(times(2,3))){21,1}.int_total.samples)...
                   /length(HI.(char(times(1,3))){21,1}.int_total.samples);  % Bayesian inference SD>WR
stats(1,4) = HI.(char(times(2,3))){21,1}.int_total.mean;
stats(1,5) = HI.(char(times(2,3))){21, 1}.int_total.var;
stats(1,6) = sum(HI.(char(times(3,3))){21,1}.int_total.samples>...
                   HI.(char(times(2,3))){21,1}.int_total.samples)...
                   /length(HI.(char(times(2,3))){21,1}.int_total.samples);  % Bayesian inference PRN>SD
stats(1,7) = HI.(char(times(3,3))){21,1}.int_total.mean;
stats(1,8) = HI.(char(times(3,3))){21, 1}.int_total.var; 

stats = array2table(stats);
stats.Properties.VariableNames = {'WRmean','WRstd','sig_SD_WR','SDmean','SDstd','sig_PN_SD','PRNmean','PRNstd'};

stats




