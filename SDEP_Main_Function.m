%% Main analysis pipeline for sleep deprivation study %%
% 
%             z             Nathan Cross
%            z              2019:2020**
%           z               
%      _ _ z
%   '( -l- )'  ___
%  /\ \ O /   /   \     |
% /__\  ^----/-----\----|
%        \-----|
%                                              **2020 doesn't really exist
%                                   
%% 1. Setup
%   i. Set variables 
subject_dir = [{'../nifti_bids/derivatives/xcpEngine/output_bids/clean_nogsr_surfsmooth6/'}]; %[{'../nifti_bids/derivatives/clean_nogsr_surfsmooth6/'};
     % <--- MUST have forward slash / at end
keyword = 'task-All'; %<-------------------------------- Set to specific task - this is case sensitive (for all tasks: keyword1 = 'task-All')
format = ['.gii'];  %<--- Set to format in which you would like to import (e.g. volumetric = '.nii'; surface = '.gii'; both = ['.nii','.gii'])
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap']; % <-------------- Set up study specific session & run combinations (LEAVE AS IS FOR SDEP)

% Needed ONLY IF keyword1 = 'task-All'
file = {'task-ANT', 'task-Nback', 'task-PVT'}; % list + order of tasks to concatenate 

%   ii. Specify Schaefer Parcellation and Yeo Network templates 
fs = 'fsaverage5'; %<-------------------------------------- (Choose number of vertices as desired)
nparc = '400'; %<------------------------------------------ (Choose number of parcellations as desired)
nnetw = '17'; %<------------------------------------------- (Choose number of networks as desired)

%  iii. Load in the Schaeffer parcellations
load(sprintf('labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw)) 
NumberOfParcels = max(rh_labels_Shf); % Number of Parcels used (saved for future)

%  iv. Load in the Yeo network templates
Yeo_7Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters007_ref.mat'); % 7 Networks
Yeo_17Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters017_ref.mat'); % 17 Networks

%	v. Map Schaeffer parcellations to the Yeo Networks (temporary workaround)
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.rh_labels(i);
end

%% 2. Import your data
%
% NOTE: If the data has been imported and Schaeffer parcellations created
% previously, you can skip ahead to Section 5.
%
% Option A. Import one task only.
%
Subject = import_data(subject_dir, keyword, format, times); % Imports data 
save(sprintf('data/Subject_%s.mat',keyword), 'Subject','-v7.3');

%% 
% (Optional). Concatenate all tasks.
%surfaces
keyword = 'task-All'; % Note: change keyword to 'task-All' 
Subject = concatenate_timeseries(file);
save(sprintf('data/Subject_%s.mat',keyword), 'Subject','-v7.3');
%% volumes
keyword = 'task-All'; % Note: change keyword to 'task-All' 
Subject = concatenate_timeseries_vols(file,times);
save(sprintf('data/Subject_vols_%s.mat',keyword), 'Subject','-v7.3');

%% Combine surfs and vols
keyword = 'task-All'; 
Global_timeseries = combine_surf_vol(keyword);
save(sprintf('data/Global_timeseries_%s.mat',keyword), 'Global_timeseries','-v7.3');
%%
Tasktimes = concatenate_taskonsets(file,times);
save('times/Tasktimes.mat', 'Tasktimes');

%% 3.1. Map data to Schaeffer parcellations
filename = sprintf('Subject_%s',keyword);
Schaeffer = map_Schaeffer_parcellations(filename,nparc,lh_labels_Shf,rh_labels_Shf);
save(sprintf('data/Schaeffer%s_%s.mat',nparc,keyword), 'Schaeffer','-v7.3');

%% 3.2. Regress out confounds of interest
filename = sprintf('Schaeffer%s_%s',nparc,keyword);
Schaeffer = task_onset_regression(filename,nparc);
save(sprintf('data/Schaeffer%s_%s.mat',nparc,keyword), 'Schaeffer','-v7.3');

%% 4.1. Create connectivity matrices between Schaeffer parcellations
filename = sprintf('Schaeffer%s_%s',nparc,keyword);
ConMat = create_connectivity_matrices(filename,nparc,lh_labels_Shf,rh_labels_Shf);
save(sprintf('data/ConMat%s_%s.mat',nparc,keyword), 'ConMat','-v7.3');
%%
CovMat = create_coraviance_matrices(filename,nparc,lh_labels_Shf,rh_labels_Shf);
save(sprintf('data/CovMat%s_%s.mat',nparc,keyword), 'CovMat','-v7.3');

%% (Optional) Load back in previously imported data
%
% You may run this if the data has previously been imported and setup
% with the Schaeffer parcellations to save time.
%
load(sprintf('data/ConMat%s_%s.mat',nparc,keyword));

%% 4.2. Reorder connectivity matrices based on Yeo & Krienen (2011) networks
%
for t = 1:length(times)
    matrix = ConMat.Session_mean_z.(char(times(t,3)));
    [matrix_ordered, matrix_bordered] = reorder_matrices(matrix, nparc, nnetw); 

    L = matrix_bordered(2,2:end);
    figure;   
    imagesc(matrix_bordered); 
    set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
    set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
    set(gca, 'YTickLabel', L, 'FontSize', 5);
    caxis([-1.5, 1.5]);
    colormap fireice; % <--- got to town and play around on this 

end

%% 5.1. Model of changes in z-scores - Fmat, Pmat, Pmat_fdr

% Apply GLM to connectivity matrices
[Fmat, Pmat] = model_rm(ConMat,nparc);
    
% Apply FDR corrections
Pvec = Pmat(:);
Pvec_fdr = fdr_bh(Pvec);
Pmat_fdr = reshape(Pvec_fdr,[NumberOfParcels, NumberOfParcels]);

%Plot F map
[matrix_ordered, matrix_bordered] = reorder_matrices(Fmat, nparc, nnetw); 
L = matrix_bordered(2,2:end);
figure;   
imagesc(matrix_bordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([-8, 8]);
colormap fireice; % <--- got to town and play around on this 

%Plot P map
[matrix_ordered, matrix_bordered] = reorder_matrices(Pmat_fdr, nparc, nnetw); 
L = matrix_bordered(2,2:end);
figure;   
imagesc(matrix_bordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([0, 0.05]);
colormap hot; % <--- got to town and play around on this 

% Save variables
GLM.Fmat = Fmat;
GLM.Pmat = Pmat;
GLM.Pmat_fdr = Pmat_fdr;

%% 5.2 Conduct post-hoc T-tests 
SubjectName = fieldnames(ConMat.times.Control);
NumberOfSubjects = length(fieldnames(ConMat.times.Control));
holdmat=NaN(NumberOfSubjects, length(times(:,3)));

for i=1:NumberOfParcels
    for j = 1:NumberOfParcels
        %if Pmat_fdr(i, j) <0.05 == 1 %check whether significant from ANOVA
            for z = 1:NumberOfSubjects
              if isfield(ConMat.times.(char(times(1,3))),(SubjectName{z})) == 1
                  if isnan(ConMat.times.(char(times(1,3))).(SubjectName{z}).ShfCorr_z(i,j)) ~= 1
                        holdmat(z,1) = ConMat.times.(char(times(1,3))).(SubjectName{z}).ShfCorr_z(i,j);
                        holdmat(z,2) = ConMat.times.(char(times(2,3))).(SubjectName{z}).ShfCorr_z(i,j);
                        holdmat(z,3) = ConMat.times.(char(times(3,3))).(SubjectName{z}).ShfCorr_z(i,j);
                  else
                    continue
                  end
              else
                  continue
              end
            end
%         else
%             continue
%         end
        a = array2table(holdmat);
        a.Properties.VariableNames = {'Session1','Session2','Session3'}; % add in the co-variates here
        s1 = holdmat(:,1);
        s2 = holdmat(:,2);
        s3 = holdmat(:,3);
        [h,p,ci,stats] = ttest(s2, s1);
        SDEP_Tmat(i,j) = stats.tstat;
        SDEP_Pmat(i,j) = p;
        [h,p,ci,stats] = ttest(s3, s2);
        Recovery_Tmat(i,j) = stats.tstat;
        Recovery_Pmat(i,j) = p;
        [h,p,ci,stats] = ttest(s3, s1);
        ConRec_Tmat(i,j) = stats.tstat;
        ConRec_Pmat(i,j) = p;
    end 
end

GLM.SDEP_Tmat = SDEP_Tmat;
GLM.SDEP_Pmat = SDEP_Pmat;
GLM.Recovery_Tmat = Recovery_Tmat;
GLM.Recovery_Pmat = Recovery_Pmat;
GLM.ConRec_Tmat = ConRec_Tmat;
GLM.ConRec_Pmat = ConRec_Pmat;

save(sprintf('data/GLM%s_%s.mat',nparc,keyword), 'GLM','-v7.3');

%Apply FDR corrections
SDEP_Pmat_l = tril(SDEP_Pmat);
SDEP_Pvec = SDEP_Pmat_l(:);
[h, crit_p, adj_ci_cvrg, SDEP_Pvec_fdr] = fdr_bh(SDEP_Pvec);
SDEP_Pmat_fdr = reshape(SDEP_Pvec_fdr,[NumberOfParcels, NumberOfParcels]);
SDEP_Pmat_fdr(isnan(SDEP_Pmat_fdr))=2; SDEP_Pmat_fdr(SDEP_Pmat_fdr==0)=2;

matrix = SDEP_Tmat;
[matrix_ordered, matrix_bordered] = reorder_matrices(matrix, nparc, nnetw);
L = matrix_bordered(2,2:end);
figure;   
imagesc(matrix_ordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([-5 5]);
colormap fireice;
