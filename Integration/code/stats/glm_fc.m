%% Setup 
NumberOfParcels = 400;
NumberOfNetworks = 17;
keyword = 'task-All';
sub_dir = '../../data/fc/';
load('../../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = length(SubjectName);
nparc = char(string(NumberOfParcels)); nnetw = char(string(NumberOfNetworks));
holdmat1=NaN(NumberOfSubjects, NumberOfParcels,NumberOfParcels);
holdmat2=NaN(NumberOfSubjects, NumberOfParcels,NumberOfParcels);
SDEP_Tmat=NaN(NumberOfParcels,NumberOfParcels); SDEP_Pmat=NaN(NumberOfParcels,NumberOfParcels); 
addpath(genpath(pwd));
%% 1. Conduct T-test between WR and SD

% 1a. Load in corr mats
for z = 1:NumberOfSubjects
    load( [sub_dir char(SubjectName{z}) '/Control/' char(SubjectName{z}) '_' keyword '.mat']) 
    holdmat1(z,:,:) = ConMat.corr_z;
    load( [sub_dir char(SubjectName{z}) '/PreNap/' char(SubjectName{z}) '_' keyword '.mat']) 
    holdmat2(z,:,:) = ConMat.corr_z;            
end

% 1b. Perform t-tests on edges
for i=1:NumberOfParcels
    for j = 1:NumberOfParcels
            v1 = holdmat1(:,i,j);
            v2 = holdmat2(:,i,j);
            [~,p,~,stats] = ttest(v2, v1);
            SDEP_Tmat(i,j) = stats.tstat;
            SDEP_Pmat(i,j) = p;
    end 
end

for i=1:NumberOfParcels
    for j = 1:NumberOfParcels
        SDEP_Tmat(j,i) = SDEP_Tmat(i,j);
    end
end

GLM.SDEP_Tmat = SDEP_Tmat;
GLM.SDEP_Pmat = SDEP_Pmat;

% 1c. Apply FDR corrections
lut = transpose(table2cell(readtable(sprintf('../../labels/Schaefer2018_%sParcels_%sNetworks_order.txt', nparc,nnetw))));
[SDEP_Pmat_ordered, ~, ~] = reorder_matrices(SDEP_Pmat, lut, nparc, nnetw);
SDEP_Pmat_l = triu(SDEP_Pmat_ordered,1);
SDEP_Pmat_l(SDEP_Pmat_l==0)=NaN;
[~, ~, ~, SDEP_Pmat_fdr] = fdr_bh(SDEP_Pmat_l);
SDEP_Pmat_fdr(isnan(SDEP_Pmat_fdr))=2; 
GLM.SDEP_Pmat_fdr_ordered = SDEP_Pmat_fdr;

%% 2. ANOVA model of changes between SD and PRN (controlling for TST)

% 2a. Set Contrasts to apply and covariates to include
contrasts = {'PreNap','PostNap'};
covariates = [{'TST'},{[150;220;100;11;201;230;260;111;289;240;206;260;233;252;184;179;205;226;140;190]}]; 

% 2b. Apply rm-GLM to connectivity matrices
[GLM.PRN_Fmat, GLM.PRN_Pmat, GLM.PRN_Diff, GLM.PRN_Tmat, GLM.PRN_PostPmat] = model_rm(sub_dir, SubjectName,contrasts, covariates);

% 2c. Apply FDR corrections
PRN_Pmat = GLM.PRN_Pmat;
[PRN_Pmat_ordered, ~, ~] = reorder_matrices(PRN_Pmat, lut, nparc, nnetw);
PRN_Pmat_l = triu(PRN_Pmat_ordered,1);
PRN_Pmat_l(PRN_Pmat_l==0)=NaN;
[~, ~, ~, PRN_Pmat_fdr] = fdr_bh(PRN_Pmat_l);
PRN_Pmat_fdr(isnan(PRN_Pmat_fdr))=2; 
GLM.PRN_Pmat_fdr_ordered = PRN_Pmat_fdr;
figure; imagesc(GLM.PRN_Pmat_fdr_ordered)


%% 3. Save output
save(sprintf('../../data/stats/GLM%s_%s_%s.mat',nparc,nnetw,keyword), 'GLM','-v7.3');

