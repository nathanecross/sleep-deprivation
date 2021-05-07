% 
% Script to run integration analyses on cortical functional data using  
%  hierarchical (Bayesian) group inference.
%
%  This script is written with the expectation that the BOLD timeseries
%  have been preprocessed and clustered into parcels or ROIs. Each ROI
%  is also required to be assigned to a higher-level/lower-order group or
%  "network" (ie. the hierarchical aspect of these integration analyses).

%% Setup
sub_dir = '../../data/Parcels_regr_cat/';
keyword = 'task-All';
fs = 'fsaverage5';                  %<-- (Choose number of vertices as desired)
nparc = '400';                      %<-- (Choose number of parcellations as desired)
NumberOfParcels = str2num(nparc);
nnetw = '7';                        %<-- (Choose number of networks as desired)
NumberOfNets = str2num(nnetw);
network_ids = 1:NumberOfNets;
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in network mask
mask_rois = Yeo_Shf; 
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
     
%% 1. Run Hierachical Integration within the whole brain

% a. Arrange data and calculate integration for subjects 
sprintf("Running Hierachical Integration within the whole brain")
addpath(genpath('../../../dependencies/'));
for t = 1:3
    for z = 1:NumberOfSubjects
        load( [sub_dir sprintf('%s/%s/%s_%s.mat',  ...
            (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword)] );
    fmri{z,1} = Parcels';
    end
    opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 800);
    sprintf("Running Session -> %s",(char(times(t,3))))
    HI.(char(times(t,3))) = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);
end

save(sprintf('../../data/integration/Integration%s_%s_%s',nparc,nnetw,kewyord),'HI');


% b. Create summary table
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

               
%% 2. Run Hierachical Integration within 17 -> 7 networks
nnetw = '17';
% a. Assign each of the 17 networks to the 7 networks
new_nets = {};
for n = 1:str2num(nnetw)
    new_nets(n,1) = {n};
end
new_nets(1:2,2) = {1}; new_nets(1:2,3) = {'Visual'}; 
new_nets(3:4,2) = {2}; new_nets(3:4,3) = {'SomMot'}; 
new_nets(5:6,2) = {3}; new_nets(5:6,3) = {'DorsAttn'}; 
new_nets(7:8,2) = {4}; new_nets(7:8,3) = {'SalVentAttn'}; 
new_nets(9:10,2) = {5}; new_nets(9:10,3) = {'Limbic'}; 
new_nets(11:13,2) = {6}; new_nets(11:13,3) = {'Cont'}; 
new_nets(17,2) = {6}; new_nets(17,3) = {'Cont'}; 
new_nets(14:16,2) = {7}; new_nets(14:16,3) = {'Default'}; 
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in Yeo network mask
network_ids = 1:7;
mask_rois_net = [];
for i=1:length(network_ids)
    for ii=1:NumberOfParcels
        mask_rois_net(ii,1) = cell2mat(new_nets(Yeo_Shf(ii),2));
    end
end

ShfLabel = [(1:NumberOfParcels)' Yeo_Shf mask_rois_net];

% b. Arrange data and calculate integration for subjects 
sprintf('Running integration within each 7 network')
for n = 1:7
    SN = ShfLabel(ShfLabel(:,3)==n,:);
    mask_rois_net = SN(:,2);
    network_ids_net = unique(SN(:,2));
    for t = 1:3
        sprintf("Running -> Network %s, Session %s",string(n),string(t))
        clear fmri
        for z = 1:NumberOfSubjects
            load( [sub_dir sprintf('%s/%s/%s_%s.mat',  ...
                (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword)] );
            fmri{z,1} = Parcels';
        end
        opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 800);
        HI.Networks.(char(times(t,3))).(char("N"+string(n))) = hierarchical_integration(fmri, network_ids_net, mask_rois_net, opt_inference, true);
    end
end

HI.Networks.Itot = Ii; HI.Networks.Ibs = Iibs; HI.Networks.Iws = Iiws; HI.Networks.FCR = FCRi;

save(sprintf('../../data/integration/Integration%s_7_17_%s',nparc,kewyord),'HI');

% c. Continue summary table
for n = 1:7
    stats(n+1,1) = HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.mean;
    stats(n+1,2) = HI.Networks.(char(times(1,3))).(char("N"+string(n))){21, 1}.int_total.var;
    stats(n+1,3) = sum(HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.samples>...
                       HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples)...
                       /length(HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.samples);
    stats(n+1,4) = HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.mean;
    stats(n+1,5) = HI.Networks.(char(times(2,3))).(char("N"+string(n))){21, 1}.int_total.var;
    stats(n+1,6) = sum(HI.Networks.(char(times(3,3))).(char("N"+string(n))){21,1}.int_total.samples>...
                       HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples)...
                       /length(HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples);
    stats(n+1,7) = HI.Networks.(char(times(3,3))).(char("N"+string(n))){21,1}.int_total.mean;
    stats(n+1,8) = HI.Networks.(char(times(3,3))).(char("N"+string(n))){21, 1}.int_total.var;    
end

%% 3. Run Hierachical Integration within each 17 network
sprintf("Running Hierachical Integration within each 17 network")

% a. Arrange data and calculate integration for subjects
nnetw = '17'; NumberOfNets=str2num(nnetw);
load(sprintf('../../data/ConMat%s_%s_task-All.mat',nparc,nnetw));
matrix = ConMat.Session_mean.Control; 
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in Yeo network mask
mask_rois = Yeo_Shf; 
matrix(ConMat.Session_mean_p.Control>0.05)=0; % threshold group average correlation matrix
ShfLabel = intraclass_clustering(matrix,nparc,nnetw); % NOTE: Best results with 400 parcels, definitely not 100!!
load(sprintf('../../data/Parcels%s_%s_task-All.mat',nparc,nnetw));
ShfLabel = [ShfLabel mask_rois];
for n = 1:NumberOfNets
    SN = ShfLabel(ShfLabel(:,3)==n,:);
    mask_rois_net = SN(:,2);
    network_ids_net = unique(SN(:,2));
    for t = 1:3
        sprintf("Running -> Network %s, Session %s",string(n),string(t))
        clear fmri
        for z = 1:NumberOfSubjects
            load( [sub_dir sprintf('%s/%s/%s_%s.mat',  ...
                (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword)] );
            fmri{z,1} = Parcels';
        end
        opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 800);
        try
            HI.Networks.(char(times(t,3))).(char("N"+string(n))) = hierarchical_integration(fmri, network_ids_net, mask_rois_net, opt_inference, true);
        catch
            opt_inference = struct('method', 'standard-inference', 'n_samplings', 1000, 'nu_max', 800);
            HI.Networks.(char(times(t,3))).(char("N"+string(n))) = hierarchical_integration(fmri, network_ids_net, mask_rois_net, opt_inference, true);
        end
    end

end
Ii = zeros(NumberOfNets,20,3); Ii_inter = zeros(NumberOfNets,20,3); Ii_intra = zeros(NumberOfNets,20,12,3); 

for n = 1:NumberOfNets
    for z = 1:20
        Ii(n,z,1) = HI.Networks.Control.(char("N"+string(n))){z, 1}.int_total.mean;
        Ii(n,z,2) = HI.Networks.PreNap.(char("N"+string(n))){z, 1}.int_total.mean;
        Ii(n,z,3) = HI.Networks.PostNap.(char("N"+string(n))){z, 1}.int_total.mean;
        Ii_inter(n,z,1) = HI.Networks.Control.(char("N"+string(n))){z, 1}.int_inter.mean;
        Ii_inter(n,z,2) = HI.Networks.PreNap.(char("N"+string(n))){z, 1}.int_inter.mean;
        Ii_inter(n,z,3) = HI.Networks.PostNap.(char("N"+string(n))){z, 1}.int_inter.mean;
        for m = 1:length(HI.Networks.Control.(char("N"+string(n))){z, 1}.int_intra.mean)
            Ii_intra(n,z,m,1) = HI.Networks.Control.(char("N"+string(n))){z, 1}.int_intra.mean(m,m);
            Ii_intra(n,z,m,2) = HI.Networks.PreNap.(char("N"+string(n))){z, 1}.int_intra.mean(m,m);
            Ii_intra(n,z,m,3) = HI.Networks.PostNap.(char("N"+string(n))){z, 1}.int_intra.mean(m,m);
        end
    end
end

for z = 1:20
    Iiws(:,z,1) = squeeze(sum(Ii_intra(:,z,:,1),3));
    Iiws(:,z,2) = squeeze(sum(Ii_intra(:,z,:,2),3));
    Iiws(:,z,3) = squeeze(sum(Ii_intra(:,z,:,3),3));
    Iibs(:,z,1) = Ii_inter(:,z,1);
    Iibs(:,z,2) = Ii_inter(:,z,2);
    Iibs(:,z,3) = Ii_inter(:,z,3);
    FCRi(:,z,1) = Iiws(:,z,1)./Iibs(:,z,1);
    FCRi(:,z,2) = Iiws(:,z,2)./Iibs(:,z,2);
    FCRi(:,z,3) = Iiws(:,z,3)./Iibs(:,z,3);
end
HI.Networks.Itot = Ii;
HI.Networks.Iws = Iiws;
HI.Networks.Ibs = Iibs;
HI.Networks.FCR = FCRi;

save(sprintf('../../data/integration/Integration%s_%s_%s',nparc,nnetw,keyword),'HI');

% b. Continue summary table
for n = 1:17
    try
        stats(n+8,1) = HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.mean;
        stats(n+8,2) = HI.Networks.(char(times(1,3))).(char("N"+string(n))){21, 1}.int_total.var;
        stats(n+8,3) = sum(HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.samples>...
                           HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples)...
                           /length(HI.Networks.(char(times(1,3))).(char("N"+string(n))){21,1}.int_total.samples);
        stats(n+8,4) = HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.mean;
        stats(n+8,5) = HI.Networks.(char(times(2,3))).(char("N"+string(n))){21, 1}.int_total.var;
        stats(n+8,6) = sum(HI.Networks.(char(times(3,3))).(char("N"+string(n))){21,1}.int_total.samples>...
                           HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples)...
                           /length(HI.Networks.(char(times(2,3))).(char("N"+string(n))){21,1}.int_total.samples);
        stats(n+8,7) = HI.Networks.(char(times(3,3))).(char("N"+string(n))){21,1}.int_total.mean;
        stats(n+8,8) = HI.Networks.(char(times(3,3))).(char("N"+string(n))){21, 1}.int_total.var;
    catch
        stats(n+8,:) = NaN;
    end  
end
s = array2table(stats);
s.Properties.VariableNames = {'WRmean','WRstd','sig','SDmean','SDstd','sig2','PRNmean','PRNstd'};

HI_stats = s;

save('../../data/stats/HI_stats.mat', 'HI_stats');


%% Stats for assemblies
ass_stat = NaN(17,5,3); ass_dfff = NaN(17,5,3); ass_diffpc = NaN(17,5,3);
for n = 1:17
    for l = 1:size(HI.Networks.Control.(char("N"+string(n))){21, 1}.int_intra.samples,1)
        Inws(l,:,1) = squeeze(HI.Networks.Control.(char("N"+string(n))){21, 1}.int_intra.samples(l,l,:));
        Inws(l,:,2) = squeeze(HI.Networks.PreNap.(char("N"+string(n))){21, 1}.int_intra.samples(l,l,:));
        Inws(l,:,3) = squeeze(HI.Networks.PostNap.(char("N"+string(n))){21, 1}.int_intra.samples(l,l,:));
        
        ass_diff(n,l,1) = mean(Inws(l,:,2)-Inws(l,:,1));
        ass_diff(n,l,2) = mean(Inws(l,:,3)-Inws(l,:,2));
        ass_dfff(n,l,3) = mean(Inws(l,:,3)-Inws(l,:,1));
        
        ass_diffpc(n,l,1) = mean((Inws(l,:,2)-Inws(l,:,1))./Inws(l,:,1))*100;
        ass_diffpc(n,l,2) = mean((Inws(l,:,3)-Inws(l,:,2))./Inws(l,:,2))*100;
        ass_dfffpc(n,l,3) = mean((Inws(l,:,3)-Inws(l,:,1))./Inws(l,:,1))*100;
        
        ass_stat(n,l,1) = sum(Inws(l,:,1)>Inws(l,:,2))/1000;
        ass_stat(n,l,2) = sum(Inws(l,:,2)>Inws(l,:,3))/1000;
        ass_stat(n,l,3) = sum(Inws(l,:,1)>Inws(l,:,3))/1000;
    end
end
%SD
HI.Assemblies.WR_SD_stats = ass_stat(:,:,1);
HI.Assemblies.WR_SD_diff = ass_diff(:,:,1);
HI.Assemblies.WR_SD_diffpc = ass_diffpc(:,:,1);
%PRN
HI.Assemblies.SD_PRN_stats = ass_stat(:,:,2);
HI.Assemblies.SD_PRN_diff = ass_diff(:,:,2);
HI.Assemblies.SD_PRN_diffpc = ass_diffpc(:,:,2);

save(sprintf('../../data/integration/Integration%s_%s_%s',nparc,nnetw,keyword),'HI');
