addpath(genpath('dependencies/'));

%% Setup
fs = 'fsaverage5'; %<-------------------------------------- (Choose number of vertices as desired)
nparc = '100'; %<------------------------------------------ (Choose number of parcellations as desired)
nnetw = '7'; %<------------------------------------------- (Choose number of networks as desired)
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];
     
NumberOfParcels = str2num(nparc);
NumberOfNets = str2num(nnetw);
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc)); %loads in Yeo network mask
network_ids = 1:NumberOfNets;
mask_rois = Yeo_Shf; 
load(sprintf('data/ConMat%s_%s_CROSS.mat',nparc,nnetw));
load(sprintf('data/Parcels%s_%s_CROSS.mat',nparc,nnetw));
lendata = length(Parcels.Control.sub01.ShfNet);
load('SubjectName.mat');  
NumberOfSubjects=length(SubjectName);
load(sprintf('labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw)) 

% SELECT ONE
Yeo_Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters007_ref.mat'); % 7 Networks
%Yeo_Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters017_ref.mat'); % 17 Networks

%	v. Map Schaeffer parcellations to the Yeo Networks (temporary workaround)
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

%% Figure S3.A
holdmat=NaN(NumberOfSubjects, 2);
for i=1:NumberOfParcels
    for j = 1:NumberOfParcels
            for z = 1:NumberOfSubjects
              if isfield(ConMat.conditions.(char(times(1,3))),char((SubjectName{z}))) == 1
                  if isnan(ConMat.conditions.(char(times(1,3))).(char(SubjectName{z})).ShfCorr_z(i,j)) ~= 1
                        holdmat(z,1) = ConMat.conditions.(char(times(1,3))).(char(SubjectName{z})).ShfCorr_z(i,j);
                        holdmat(z,2) = ConMat.conditions.(char(times(2,3))).(char(SubjectName{z})).ShfCorr_z(i,j);
                  else
                    continue
                  end
              else
                  continue
              end
            end
            a = array2table(holdmat);
            a.Properties.VariableNames = {'Session1','Session2'}; % add in the co-variates here
            s1 = holdmat(:,1);
            s2 = holdmat(:,2);
            [~,p,~,stats] = ttest(s2, s1);
            SDEP_Tmat(i,j) = stats.tstat;
            SDEP_Pmat(i,j) = p;
    end 
end

for i=1:NumberOfParcels
    for j = 1:NumberOfParcels
        SDEP_Tmat(j,i) = SDEP_Tmat(i,j);
    end
end

matrix = SDEP_Tmat;
[matrix_ordered, matrix_bordered, L] = reorder_matrices(matrix, nparc, nnetw); 

figure;   
imagesc(matrix_bordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([-7, 7]);
colormap(cmap); % <--- go to town and play around on this 

matrix(SDEP_Pmat>0.05)=0;
[matrix_ordered, matrix_bordered, L] = reorder_matrices(matrix, nparc, nnetw); %'#ad2f00'
cmap = customcolormap([0 0.2 0.35 0.5 0.8 0.9 1], {'#6b1200','#E23603','#FF8D33',  '#000000', '#336EFF', '#336EFF', '#033c9e'});
figure;   
imagesc(matrix_bordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([-7, 7]);
colormap(cmap); % <--- go to town and play around on this 


%% Run Hierachical Integration within the whole brain

% 1. Arrange data and calculate integration for subjects 
for t = 1:3
    % a. Organise BOLD timeseries data into matrix (time x parcels)
    for z = 1:20
    mat  = zeros(lendata,NumberOfParcels);
        for p = 1:NumberOfParcels
        mat(:,p) = Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfNet(p,:)';
        end
    fmri{z,1} = mat; % put all subjects' matrices into cell array
    end
    % b. Specify options and calculate integration on data matrix
    opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 200);
    HI.(char(times(t,3))) = hierarchical_integration(fmri, network_ids, mask_rois, opt_inference, true);
end

% 2. Reorganise output into summaries
Itot = zeros(20,3); I_inter = zeros(20,3); I_intra = zeros(20,NumberOfNets,3);
for z = 1:20
    % Total integration
    Itot(z,1) = HI.Control{z, 1}.int_total.mean;
    Itot(z,2) = HI.PreNap{z, 1}.int_total.mean;
    Itot(z,3) = HI.PostNap{z, 1}.int_total.mean;
    % Between-systems integration
    I_inter(z,1) = HI.Control{z, 1}.int_inter.mean;
    I_inter(z,2) = HI.PreNap{z, 1}.int_inter.mean;
    I_inter(z,3) = HI.PostNap{z, 1}.int_inter.mean;
    % Within-systems integration
    for n = 1:NumberOfNets
        I_intra(z,n,1) = HI.Control{z, 1}.int_intra.mean(n,n);
        I_intra(z,n,2) = HI.PreNap{z, 1}.int_intra.mean(n,n);
        I_intra(z,n,3) = HI.PostNap{z, 1}.int_intra.mean(n,n);
    end
end

I_intra_mean(:,1) = mean(squeeze(I_intra(:,:,1)),1)';
I_intra_mean(:,2) = mean(squeeze(I_intra(:,:,2)),1)';
I_intra_mean(:,3) = mean(squeeze(I_intra(:,:,3)),1)';
I_intra_diff(:,1) = I_intra_mean(:,2) - I_intra_mean(:,1);
I_intra_diff(:,2) = I_intra_mean(:,3) - I_intra_mean(:,2);
I_intra_diff(:,3) = I_intra_mean(:,3) - I_intra_mean(:,1);

Yeo_labels = [lh_labels_Shf_Yeo; rh_labels_Shf_Yeo];
for n = 1:NumberOfNets
    for y = 1:20484
        if Yeo_labels(y,2) == n
            Yeo_labels(y,3) = I_intra_diff(n,1);
        end
    end
end

for z = 1:20
    Iws(z,1) = squeeze(sum(I_intra(z,:,1)));
    Iws(z,2) = squeeze(sum(I_intra(z,:,2)));
    Iws(z,3) = squeeze(sum(I_intra(z,:,3)));
    Ibs(z,1) = I_inter(z,1);
    Ibs(z,2) = I_inter(z,2);
    Ibs(z,3) = I_inter(z,3);
    FCR(z,1) = Iws(z,1)/Ibs(z,1);
    FCR(z,2) = Iws(z,2)/Ibs(z,2);
    FCR(z,3) = Iws(z,3)/Ibs(z,3);
end

HI.Itot = Itot; HI.Iws = Iws; HI.Ibs = Ibs; HI.FCR = FCR;
save(sprintf('data/Integration%s_%s_CROSS.mat',nparc,nnetw), 'HI');


%% For 17 -> 7 networks
nnetw = '17';
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
% Assign each of the 17 networks to the 7 networks
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

network_ids = 1:7;
mask_rois_net = [];
for i=1:length(network_ids)
    for ii=1:NumberOfParcels
        mask_rois_net(ii,1) = cell2mat(new_nets(Yeo_Shf(ii),2));
    end
end

ShfLabel = [(1:NumberOfParcels)' Yeo_Shf mask_rois_net];

% Run integration within each 7 network
for n = 1:7
    SN = ShfLabel(ShfLabel(:,3)==n,:);
    mask_rois_net = SN(:,2);
    network_ids_net = unique(SN(:,2));
    for t = 1:3
        for z = 1:20
        mat  = zeros(lendata,length(SN(:,1)));
            for p = 1:length(SN(:,1))
            mat(:,p) = Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfNet(SN(p,1),:)';
            end
        fmri{z,1} = mat;
        end
        opt_inference = struct('method', 'hierarchical-inference', 'n_samplings', 1000, 'nu_max', 200);
        HI.Networks.(char(times(t,3))).(char("N"+string(n))) = hierarchical_integration(fmri, network_ids_net, mask_rois_net, opt_inference, true);
    end

end

% Obtain results
Ii = zeros(7,20,3); Ii_inter = zeros(7,20,3); Ii_intra = zeros(7,20,6,3); Ii_intra_mean = zeros(7,6,3);

for n = 1:7
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
    Ii_intra_mean(n,:,1) = mean(squeeze(Ii_intra(n,:,:,1)),1)';
    Ii_intra_mean(n,:,2) = mean(squeeze(Ii_intra(n,:,:,2)),1)';
    Ii_intra_mean(n,:,3) = mean(squeeze(Ii_intra(n,:,:,3)),1)';
    Ii_intra_diff(n,:,1) = Ii_intra_mean(n,:,2) - Ii_intra_mean(n,:,1);
    Ii_intra_diff(n,:,2) = Ii_intra_mean(n,:,3) - Ii_intra_mean(n,:,2);
    Ii_intra_diff(n,:,3) = Ii_intra_mean(n,:,3) - Ii_intra_mean(n,:,1);
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


for n = 1:7
    for t = 1:3
        I_tab(n,t*2-1) = mean(Ii(n,:,t));
        I_tab(n,t*2) = std(Ii(n,:,t));
        Ibs_tab(n,t*2-1) = mean(Iibs(n,:,t));
        Ibs_tab(n,t*2) = std(Iibs(n,:,t));
        Iws_tab(n,t*2-1) = mean(Iiws(n,:,t));
        Iws_tab(n,t*2) = std(Iiws(n,:,t));
        FCR_tab(n,t*2-1) = mean(FCRi(n,:,t));
        FCR_tab(n,t*2) = std(FCRi(n,:,t));
    end
end
nnetw='7';
HI.Networks.Itot = Ii; HI.Networks.Ibs = Iibs; HI.Networks.Iws = Iiws; HI.Networks.FCR = FCRi;

save(sprintf('data/Integration%s_%s_CROSS.mat',nparc,nnetw), 'HI');

%% Figure S3.B 
% Hierarchical integration
clear a
load('data/Integration100_7_CROSS.mat');
cmap = customcolormap([0 0.2 0.4 0.5 0.8 0.9 1], {'#6b1200','#E23603','#FF8D33',  '#bababa', '#336EFF', '#336EFF', '#033c9e'});

% [~,~,~,stats] = ttest(HI.Itot(:,2),HI.Itot(:,1));
% data(1,1) = stats.tstat;
% order7 = [6,7,3,5,4,1,2]; %reorder 7 networks
% order17 = [11,12,13,17,14,15,16,5,6,9,10,7,8,3,4,1,2];
% for n = 1:7
%     [~,~,~,stats] = ttest(HI.Networks.Itot(n,:,2),HI.Networks.Itot(n,:,1));
%     data(order7(n),2) = stats.tstat;
% end
% for n =1:17
%     [~,~,~,stats] = ttest(HI_17.Networks.Itot(n,:,2),HI_17.Networks.Itot(n,:,1));
%     data(order17(n),3) = stats.tstat;
% end
data = NaN(7,2);
d = (HI.PreNap{21, 1}.int_total.mean-HI.Control{21, 1}.int_total.mean)/...
    HI.Control{21, 1}.int_total.mean*100;
data(1:7,1) = [d; NaN(6,1)];

for n=1:7
    data(n,2) = ((HI.Networks.PreNap.(char("N" + n)){21, 1}.int_total.mean - ...
        HI.Networks.Control.(char("N" + n)){21, 1}.int_total.mean))/...
        HI.Networks.Control.(char("N" + n)){21, 1}.int_total.mean*100;
end

data(1:7,2) = [data(6,2);data(7,2);data(3,2);data(5,2);data(4,2);data(1,2);data(2,2)];
i = 0;
for c = 1:size(data,2)
    for r = 1:length(data(:,c))
        if isnan(data(r,c)) ~= 1
            i = i+1;
            dat(i,1) = data(r,c);
        end 
    end
end

addpath('dependencies/customcolormap')
addpath('dependencies/fireice')
clim = [-50 50];


f = figure;

% cortex
a(1) = axes('position', [0.1 0.05 0.1 0.9]);
imagesc(dat(1)); axis off;
colormap(a(1), cmap)
caxis(a(1), clim)

% networks 1:7
for i = 2:8
    if i == 2
        a(i) = axes('position', [0.3 1.01-(0.13*i) 0.1 0.2]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    elseif i == 3
        a(i) = axes('position', [0.3 0.95-(0.13*i) 0.1 0.15]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    else
        a(i) = axes('position', [0.3 0.85-(0.1*i) 0.1 0.08]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    end
end


% S3.C Plot on surface
SP = SurfStatAvSurf({'../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/lh.pial', '../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/rh.pial'});

% 7 networks
oldorder = [data(7,2);data(6,2);data(3,2);data(5,2);data(4,2);data(1,2);data(2,2)];
nnetw = '7'; nparc ='400'; NumberOfParcels = str2num(nparc);
load(sprintf('labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end
for i = 1:NumberOfParcels
    vec(i,1) = oldorder(Yeo_Shf(i),1);   
end

parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
%nonwall = [2:50 51:101];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.5 0.70 0.3 0.3], [0.5 0.47 0.3 0.3], ...
            1, 2, clim, cmap);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.5 0.24 0.3 0.3], [0.5 0.01 0.3 0.3], ...
            1, 2, clim, cmap);

a(10) = axes('position', [0.9 0.05 0.05 0.9]);
imagesc((flipud([clim(1):0.01:clim(end)]'))); axis off
colormap(a(10), cmap)
caxis(a(10), clim);  
