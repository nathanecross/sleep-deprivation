% Script to generate Figure 5 in Cross et al. 2021 (PLOS Biology)
% 
% NOTE: In order for this script to run, the repository must first be cloned,
%       and the working directory must be set to the same directory in which 
%       this script is located.


%% Setup
load('../data/subject.mat');
SubjectName = subject.SubjectList;
keyword = 'task-All'; NumberOfSubjects=20;
NumberOfParcels = 400; nparc = '400'; nnetw= '17';
vols = {'amygdala_lh','amygdala_rh',...
        'thalamus_lh','thalamus_rh',...
        'caudate_lh','caudate_rh',...
        'hippocampus_lh','hippocampus_rh',...
        'pallidum_lh','pallidum_rh',...
        'putamen_lh','putamen_rh'};
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];
addpath(genpath('../'));

%% (Pre) Create cortical-subcortical connectivty matrices
%
%   The parcelwise timeseries for each individual were generated using the following scripts:
%
%   The timeseries of BOLD data for each individual per state: 
%       e.g. '../data/Subjectdata_nii_regr_cat/sub-01/Control/sub-01_task-All.mat'         
%       and
%       e.g. '../data/Parcels_regr_cat/sub-01/Control/sub-01_task-All.mat'         
%       were generated with the script: '../code/import_setup.m'
%
 
 % WR, SD and PRN
 for t = 1:length(times)
    for z = 1:NumberOfSubjects
        load(sprintf('../data/Subjectdata_nii_regr_cat/%s/%s/%s_%s.mat', ...
             (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword));
        load(sprintf('../data/Parcels_regr_cat/%s/%s/%s_%s.mat', ...
             (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword)); 
        for p=1:NumberOfParcels
            for v = 1:length(vols) 
                 [r,~] = corrcoef(Parcels(p,:),subjectdata.(char(vols(v))));
                 submat(t,z,p,v) = r(2);
            end
        end
    end
 end
 
% Nap
for z = 1:NumberOfSubjects
    if contains((char(SubjectName{z})),'17') == 1 
        submat(4,z,:,:) = NaN;
    elseif contains((char(SubjectName{z})),'34')
        submat(4,z,:,:) = NaN;
    else
        load(sprintf('../data/Subjectdata_nii_cat/%s/nap_trunc/%s_nap.mat', ...
             (char(SubjectName{z})),(char(SubjectName{z}))));
        load(sprintf('../data/Parcels_cat/%s/nap_trunc/%s_nap.mat', ...
             (char(SubjectName{z})),(char(SubjectName{z})))); 
        for p=1:NumberOfParcels
            for v = 1:length(vols) 
                 [r,~] = corrcoef(Parcels(p,:),subjectdata.(char(vols(v))));
                 submat(4,z,p,v) = r(2);
            end
        end
    end
end

%Resting state
 keyword = 'task-CROSS';
 NumberOfParcels = 100;
 for t = 1:length(times)
    for z = 1:NumberOfSubjects
        load(sprintf('../data/Subjectdata_nii/%s/%s/%s_%s.mat', ...
             (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword));
        load(sprintf('../data/Parcels/%s/%s/%s_%s.mat', ...
             (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword)); 
        for p=1:NumberOfParcels
            for v = 1:length(vols) 
                 [r,~] = corrcoef(Parcels(p,:),subjectdata.(char(vols(v))));
                 submat_rs(t,z,p,v) = r(2);
            end
        end
    end
 end

 
%% 5.A. Ttests for differences across states
NumberOfParcels = 400;
for p=1:NumberOfParcels
    for v = 1:length(vols)
            for z = 1:NumberOfSubjects
                holdmat(z,1) = submat(1,z,p,v);
                holdmat(z,2) = submat(2,z,p,v);
                holdmat(z,3) = submat(3,z,p,v);
                holdmat(z,4) = submat(4,z,p,v);
            end
            a = array2table(holdmat);
            a.Properties.VariableNames = {'WR','SD','PRN','NREM'}; % add in the co-variates here
            s1 = holdmat(:,1);
            s2 = holdmat(:,2);
            s3 = holdmat(:,3);
            s4 = holdmat(:,4);
            [~,pval,~,stats] = ttest(s2, s1);
            SDEP_vols_Tmat(p,v) = stats.tstat;
            SDEP_vols_Pmat(p,v) = pval;
            [~,pval,~,stats] = ttest(s4, s2);
            Recovery_vols_Tmat(p,v) = stats.tstat;
            Recovery_vols_Pmat(p,v) = pval;
            [~,pval,~,stats] = ttest(s4, s1);
            ConRec_vols_Tmat(p,v) = stats.tstat;
            ConRec_vols_Pmat(p,v) = pval;
            [~,pval,~,stats] = ttest(s3, s1);
            NREM_vols_Tmat(p,v) = stats.tstat;
            NREM_vols_Pmat(p,v) = pval;
    end 
end

% Plot the correlation matrix
I = customcolormap([0 0.2 0.35 0.5 0.8 0.9 1], ...
    {'#6b1200','#E23603','#FF8D33',  '#000000', '#336EFF', '#336EFF', '#033c9e'});
lut = transpose(table2cell(readtable(sprintf('../../labels/Schaefer2018_%sParcels_%sNetworks_order.txt', ...
                                             nparc,nnetw))));

% Reorder vectors for visualisation so that lh and rh regions are together                                         
matrix = SDEP_vols_Tmat';
for l = 1:length(matrix(:,1))
         new = reorder_vectors(matrix(l,:)',nparc,nnetw,lut);
         matrix(l,:) = new';
end

% Plot
figure;   
imagesc(matrix');
set(gca, 'XTickLabel', vols, 'FontSize', 5);
caxis([-8, 8]);
colormap(I); % <--- go to town and play around on this 


%% Figure 5.B - Parcelwise change in thalamocortical connectivity
load(sprintf('../labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
Yeo_17Clusters_ref = load(sprintf('../labels/fsaverage5/YeoNetworks/1000subjects_clusters%s_ref.mat',nnetw)); % 17 Networks
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.rh_labels(i);
end
I = customcolormap([0 0.2 0.4 0.5 0.6 0.8 1], ...
[1,1,0.75; 1,1,0.046875; 1,0.03125,0; 0.729,0.729,0.729; 0,0,0.9375; ...
0,0.96875,1; 0.75,1,1]);
I = customcolormap([0 0.2 0.35 0.5 0.8 0.9 1], {'#6b1200','#E23603','#FF8D33',  '#bababa', '#336EFF', '#336EFF', '#033c9e'});
figure;
vec = mean(SDEP_vols_Tmat(:,3:4),2);
sig = SDEP_vols_Pmat(:,3:4);
sig = sig<0.05;
sig = sum(sig,2);
vec(sig<2)=0;
SP = SurfStatAvSurf({'../labels/fsaverage5/surf/lh.pial', ...
                     '../labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 6);
clim = [-9 9];
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.2 0.75 0.2 0.2], [0.39 0.75 0.2 0.2], ...
            1, 2, clim, I);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.2 0.59 0.2 0.2], [0.39 0.59 0.2 0.2], ...
            1, 2, clim, I); 

%% Figure 5.C  Thalamocortical connectivity and integration
%
%   Hierarchical models of integration were generated using the following scripts:
%
%   Bayesian inference statistics of integration changes: 
%       e.g. '../data/stats/Integration400_7_task-All.mat'           
%       was generated by: '../code/integration/run_integration_h.m'

clear I17_changevec
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}];
load(sprintf('../data/integration/Integration%s_%s_task-All.mat',nparc,nnetw)); HI_17 = HI;
load(sprintf('../labels/Yeo%s_Shf%s.mat', nnetw, nparc))
nparc = '400'; NumberOfParcels = str2num(nparc);
nnetw = '17'; NumberOfNetworks = str2num(nnetw);
NumberOfSubjects=20;

% Difference in thalamocortical connectivity: WR -> SD
for p=1:NumberOfParcels
    for v = 1:length(vols)
        for z = 1:NumberOfSubjects
                    diffmat(z,v,p) = submat(2,z,p,v) - submat(1,z,p,v);
        end   
    end
end
thaldiffmat = squeeze(mean(diffmat(:,3:4,:),2));


% Difference in integraton within 17 networks: WR -> SD
for n = 1:NumberOfNetworks
    for z = 1:20
        for t = 1:2
            holdmat(z,t) = HI_17.Networks.(char(times(t,3))).(char("N" + string(n))){z, 1}.int_total.mean;
        end
    end
    I17_change(n,:) = holdmat(:,2) - holdmat(:,1);
end
for p = 1:NumberOfParcels
    I17_changevec(:,p) = I17_change(Yeo_Shf(p),:);   
end
for p = 1:NumberOfParcels
    [rval, pval] = corrcoef(I17_changevec(:,p), thaldiffmat(:,p));
    rvec(p,1) = rval(2);
    pvec(p,1) = pval(2);
end

vec = rvec;
vec(pvec>0.05)=0;

% Plot the results
hotred = customcolormap([0 0.25 0.5 0.75 1],{'#660320','#DD7661','#F8F3F0','#4692C0','#093262'});
clim = [-1 1];
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.2 0.75 0.2 0.2], [0.39 0.75 0.2 0.2], ...
            1, 2, clim, hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.2 0.59 0.2 0.2], [0.39 0.59 0.2 0.2], ...
            1, 2, clim, hotred);  

%% Figure 5.D Thalamocortical connectivity
thal_all_task = mean(submat(:,:,:,3:4),4);
thal_all_task = mean(thal_all_task,3)';
thal_all_task = [thal_all_task(:,1:2),thal_all_task(:,4),thal_all_task(:,3)]; %reorder for plotting
thal_all_rs = mean(submat_rs(:,:,:,3:4),4);
thal_all_rs = mean(thal_all_rs,3)';

bars = [thal_all_task(:,1) thal_all_rs(:,1) thal_all_task(:,2) thal_all_rs(:,2) ...
    thal_all_task(:,3) thal_all_task(:,4) thal_all_rs(:,3)];
col = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.8078 0.1098 0.1882; ...
           0.9 0.9 0.9; 0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250];
colflip = flip(col);
%boxplots
violin(bars,'mc',[],'medc','k','facecolor',col,'facealpha',1);
h = findobj(gca,{'Type','Patch'}) %,'-and',{'Tag','Box'});
for j=1:length(h)
    hp = patch(get(h(j),'XData'),get(h(j),'YData'),colflip(j,:),'FaceAlpha',0.2);
    if j == 1 ||  j == 4 || j == 6
        hatchfill2(hp,'single','HatchAngle',45,'HatchDensity',100,'LineWidth',1,'HatchLineWidth',1.2)
    end
end
hold on; 
plot([0 8],[0,0],'k--')

for c = 1:size(bars,2)
    [~,p,~,~] = ttest(bars(:,c))
end
