%  (Nap)
%  Script to run integration analyses on cortical functional data using  
%  standard inference to create subject-level measures of integration.
%
%  This script is written with the expectation that the BOLD timeseries
%  have been preprocessed and clustered into parcels or ROIs. Each ROI
%  is also required to be assigned to a higher-level/lower-order group or
%  "network" (ie. the hierarchical aspect of these integration analyses).

%% Setup
sub_dir = '../../data/Parcels100/';
fs = 'fsaverage5';                  %<-- (Choose number of vertices as desired)
nparc = '100';                      %<-- (Choose number of parcellations as desired)
NumberOfParcels = str2num(nparc);
nnetw = '7';                        %<-- (Choose number of networks as desired)
NumberOfNets = str2num(nnetw);
network_ids = 1:NumberOfNets;
load(sprintf('../../labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in network mask
mask_rois = Yeo_Shf; 
load(sprintf('../../labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw))  %loads in parcel mask
times = [{'ses-01','run-01','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}];
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
     
%% 1. Run Hierachical Integration within the whole brain

% a. Arrange data and calculate integration for subjects 
sprintf("Running Hierachical Integration within the whole brain")
addpath(genpath('../../../dependencies/'));
clear fmri
for z = 1:NumberOfSubjects
    filepath = [sub_dir sprintf('%s/%s/', (char(SubjectName{z})),'nap_trunc')];
    files = dir(filepath); files = struct2table(files);
    files = files.name; 
    if isempty(files)
        sprintf('Skipping %s',char(SubjectName{z}))
    else
        filename = char(files(contains(files,'N2')));
        filename = [filepath filename];
        load( filename );
        fmri{z,1} = Parcels';
    end
end
mask = cellfun(@isempty,fmri);
fmri = fmri(~mask);
opt_inference = struct('method', 'standard-inference', 'n_samplings', 1000, 'nu_max', 800);
fprintf("Running Stage2 \n" )
HI.Stage2 = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);

clear fmri

for z = 1:NumberOfSubjects
    filepath = [sub_dir sprintf('%s/%s/', (char(SubjectName{z})),'nap_trunc')];
    files = dir(filepath); files = struct2table(files);
    files = files.name; 
    if isempty(files)
        sprintf('Skipping %s',char(SubjectName{z}))
    else
        filename = char(files(contains(files,'N3')));
        if isempty(filename)
            sprintf('Skipping %s',char(SubjectName{z}))
        else
            filename = [filepath filename];
            load( filename );
            fmri{z,1} = Parcels';
        end
    end
end
mask = cellfun(@isempty,fmri);
fmri = fmri(~mask);
opt_inference = struct('method', 'standard-inference', 'n_samplings', 1000, 'nu_max', 800);
fprintf("Running Stage3 \n" )
HI.Stage3 = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);



sum((HI.Stage2{17,1}.int_total.samples>...
                   HI.Stage3{12,1}.int_total.samples)...
                   /1000)
for z = 1:18
    if z<13
    FCR_N3(z) = sum(diag(HI.Stage3{z, 1}.int_intra.mean))/(HI.Stage3{z, 1}.int_inter.mean);
    end
    FCR_N2(z) = sum(diag(HI.Stage2{z, 1}.int_intra.mean))/(HI.Stage2{z, 1}.int_inter.mean);
end

FCR_N3 = [NaN(1,1) FCR_N3(1:2) NaN(1,1) FCR_N3(3:6) NaN(1,1) FCR_N3(7:9) ...
          FCR_N3(10) NaN(1,1) FCR_N3(11) NaN(1,1) FCR_N3(12) NaN(1,1)];
plot([FCR_N2' FCR_N3'])
% fprintf("Saving Data For Whole Cortex -> 7 Networks");
% save(sprintf('../../data/integration/Integration%s_7-nap.mat',nparc), 'HI');  

               
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
        for z = 1:20
        mat  = zeros(624,length(SN(:,1)));
            for p = 1:length(SN(:,1))
            mat(:,p) = Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfMeanclean.(char("p"+SN(p,1)))';
            end
        fmri{z,1} = mat;
        end
        opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 800);
        HI.Networks.(char(times(t,3))).(char("N"+string(n))) = hierarchical_integration(fmri, network_ids_net, mask_rois_net, opt_inference, true);
    end

end

HI.Networks.Itot = Ii; HI.Networks.Ibs = Iibs; HI.Networks.Iws = Iiws; HI.Networks.FCR = FCRi;
save(sprintf('../../data/Integration%s_7_17_nap.mat',nparc),'HI');

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
load(sprintf('../../data/Parcels%s_%s_nap.mat',nparc,nnetw));
ShfLabel = [ShfLabel mask_rois];
for n = 1:NumberOfNets
    SN = ShfLabel(ShfLabel(:,3)==n,:);
    mask_rois_net = SN(:,2);
    network_ids_net = unique(SN(:,2));
    for t = 1:3
        sprintf("Running -> Network %s, Session %s",string(n),string(t))
        for z = 1:20
            mat  = zeros(624,length(SN(:,1)));
            for p = 1:length(SN(:,1))
                mat(:,p) = Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfMeanclean.(char("p"+SN(p,1)))';
            end
            fmri{z,1} = mat;
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

save(sprintf('../../data/Integration%s_17_nap.mat',nparc), 'HI');