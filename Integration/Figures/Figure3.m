%% Setup
keyword = 'task-All';
nparc = '400';
NumberOfParcels = str2num(nparc);
nnetw = '17';
NumberOfNets = str2num(nnetw); 
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];
load('../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = subject.NumberOfSubjects;


% Map Schaeffer parcellations to the Yeo Networks
load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
load(sprintf('../labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw))
Yeo_Clusters_ref = load(sprintf('../labels/fsaverage5/YeoNetworks/1000subjects_clusters%s_ref.mat',nnetw)); % 17 Networks
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
end

addpath(genpath('../code/dependencies'));

%% (Pre) Calculate amplitude fluctuation 
for t = 1:3
    sprintf('Timepoint: %s',(char(times(t,3))))
    for z = 1:NumberOfSubjects
        sprintf('Subject: %s',(char(SubjectName{z})))
        load(sprintf('../data/Parcels_regr_cat/%s/%s/%s_%s.mat', ...
                     (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword));
        for p=1:NumberOfParcels
                STDDEV_surf_Parcels(t,z,p,:) = std(Parcels(p,:));  
        end
        
        % 7 Networks 
        nnetw = '7';
        NumberOfNets = str2num(nnetw); 
        load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
        for n = 1:NumberOfNets
            ii = 0;
            for p=1:NumberOfParcels
                if Yeo_Shf(p) == n
                        ii = ii + 1;
                        Anet7(n,p,:) = Parcels(p,:);
                end
            end
            STDDEV_surf_Nets7(t,z,n) = std(mean(squeeze(Anet7(n,:,:)),1));
        end
        
         % 7 Networks 
        nnetw = '17';
        NumberOfNets = str2num(nnetw); 
        load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
        for n = 1:NumberOfNets
            ii = 0;
            for p=1:NumberOfParcels
                if Yeo_Shf(p) == n
                        ii = ii + 1;
                        Anet17(n,p,:) = Parcels(p,:);
                end
            end
            STDDEV_surf_Nets17(t,z,n) = std(mean(squeeze(Anet17(n,:,:)),1));
        end
        
    end
end


%% Figure 3.A Changes in parcelwise signal fluctuation WR -> SD 

% Setup variables
holdmat=NaN(NumberOfSubjects, 2);
SDEP_Tmat=NaN(NumberOfParcels,1); SDEP_Pmat=NaN(NumberOfParcels,1); 
Recovery_Tmat=NaN(NumberOfParcels,1); Recovery_Pmat=NaN(NumberOfParcels,1);
ConRec_Tmat=NaN(NumberOfParcels,1); ConRec_Pmat=NaN(NumberOfParcels,1);

% Statistical tests for change in amplitude fluctuation WR -> SD 
for i=1:NumberOfParcels
    for z = 1:NumberOfSubjects
      holdmat(z,1) = STDDEV_surf_Parcels(1,z,i);
      holdmat(z,2) = STDDEV_surf_Parcels(2,z,i);
    end
    a = array2table(holdmat);
    a.Properties.VariableNames = {'WR','SD'};
    s1 = holdmat(:,1);
    s2 = holdmat(:,2);
    [~,p,~,stats] = ttest(s2, s1);
    SDEP_Tmat(i,1) = stats.tstat;
    SDEP_Pmat(i,1) = p;
end


% Plot labels to surface
hotgrey = customcolormap([0 0.1 0.2 0.5 0.8 1], {'#FFFEE6','#FFF000', '#FEBC00', '#807F7F', '#0004FE','#00EBFE' });
cmap_grad = cbrewer('div', 'Spectral', 256, 'linear');
clim = [-3 3];
SP = SurfStatAvSurf({'../labels/fsaverage5/surf/lh.pial', ...
                     '../labels/fsaverage5/surf/rh.pial'});
vec = SDEP_Tmat;
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 8);
f = figure('units', 'centimeters', 'position', [0 0 20 20]); 
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.2 0.75 0.2 0.2], [0.4 0.75 0.2 0.2], ...
            1, 2, clim, hotgrey);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.2 0.59 0.2 0.2], [0.4 0.59 0.2 0.2], ...
            1, 2, clim, hotgrey);

%% Figure 3.B - Scatterplots of Integration Vs Amplitude of fluctuations

% Load in total integration stats for 7 and 17 networks
load('../data/stats/HI_stats.mat');
SP = SurfStatAvSurf({'../../labels/fsaverage5/surf/lh.pial', ...
                     '../../labels/fsaverage5/surf/rh.pial'});


% i. Scatterplot: ? Integration vs ? Amplitude fluctuation, 7 Networks
nnetw = '7';
load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
% Arrange data vector
data(1:7,1) = (HI_stats.SDmean(2:8) - HI_stats.WRmean(2:8)); data(8:17,1:2) = NaN;
for i = 1:NumberOfParcels
    vec(i,1) = data(Yeo_Shf(i),1);   
end
figure;
% Set colormap
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
Yeo7 = customcolormap([0 0.166 0.333 0.5 0.666 0.832 1], ...
    {Seven,Six,Five,Four,Three,Two,One});
% Plot
scatter(vec,SDEP_Tmat, 300, Yeo_Shf, 'filled', 'MarkerEdgeColor', 'k','LineWidth',0.1); lsline
colormap(Yeo7)
[r,~] = corrcoef(vec,SDEP_Tmat);
rsq = r(2)^2
% Spin test for permutation testing of significance
toMap(nonwall) = SDEP_Tmat;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
readleft = sprintf('../data/integration/I%s_changeveclh.csv',nnetw);
readright = sprintf('../data/integration/I%s_changevecrh.csv',nnetw);
pval = pvalvsNull(readleft,readright,OnSurf(1:10242)',OnSurf(10243:end)',...
                  100,sprintf('../data/stats/I%s_change_perm.mat',nnetw))

% ii. Scatterplot: ? Integration vs ? Amplitude fluctuation, 17 Networks
clear data vec
nnetw = '17';
load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));
% Arrange data vector
data(1:17,1) = (HI_stats.SDmean(9:25) - HI_stats.WRmean(9:25));
for i = 1:NumberOfParcels
    vec(i,1) = data(Yeo_Shf(i),1);   
end
figure;
% Set colormap
One = '#3D0042'; Two = '#801982'; Three = '#5DACDB'; Four = '#1A43DB'; Five = '#459738';
Six = '#45D738'; Seven = '#FF00BF'; Eight = '#FFA5FC'; Nine = '#F5FCCE'; Ten = '#FFFC84';
Eleven = '#F4B750'; Twelve = '#F49E00'; Thirteen = '#FF8832'; Fourteen = '#AE3C43';
Fifteen = '#FF3C43'; Sixteen = '#FF0000'; Seventeen = '#000000';
Yeo17 = customcolormap([0 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5 ...
    0.5625 0.625 0.6875 0.75 0.8125 0.875 0.9375 1], ...
    {Seventeen, Sixteen, Fifteen, Fourteen,Thirteen,Twelve,Eleven,Ten, Nine, Eight, ...
    Seven,Six,Five,Four,Three,Two,One});
% Plot
scatter(vec,SDEP_Tmat, 200, Yeo_Shf, 'filled', 'MarkerEdgeColor', 'k','LineWidth',0.1); lsline
colormap(Yeo17)
[r,p] = corrcoef(vec,SDEP_Tmat);
rsq = r(2)^2
% Spin test for permutation testing of significance
toMap(nonwall) = SDEP_Tmat;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
readleft = sprintf('../data/integration/I%s_changeveclh.csv',nnetw);
readright = sprintf('../data/integration/I%s_changevecrh.csv',nnetw);
pval = pvalvsNull(readleft,readright,OnSurf(1:10242)',OnSurf(10243:end)',...
                  100,sprintf('../data/stats/I%s_change_perm.mat',nnetw))
              
              
% iii. Scatterplot: ? Integration vs ? Amplitude fluctuation, 57 clusters
clear data vec
nnetw = '17';
load(sprintf('../labels/Yeo%s_Shf%s_Ass.mat',nnetw,nparc));
% Arrange data vector
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end
for i = 1:NumberOfParcels
    x = ShfLabel(i,3);
    y = num2str(ShfLabel(i,2)); y=str2num(y(end));
    vec(i,1) = HI.Assemblies.WR_SD_diff(x,y);
end
figure;
scatter(vec,SDEP_Tmat, 150, [95/255, 97/255, 99/255], 'filled', 'MarkerEdgeColor', 'k','LineWidth',0.1); lsline
[r,p] = corrcoef(vec,SDEP_Tmat);
rsq = r(2)^2
% Spin test for permutation testing of significance
toMap(nonwall) = SDEP_Tmat;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
readleft = sprintf('../data/integration/I%s_changeveclh.csv',nnetw);
readright = sprintf('../data/integration/I%s_changevecrh.csv',nnetw);
pval = pvalvsNull(readleft,readright,OnSurf(1:10242)',OnSurf(10243:end)',...
                  100,sprintf('../data/stats/I%s_change_perm.mat',nnetw))
              

%% Figure 3.C

% 17 Networks
nparc = '400'; NumberOfParcels = str2num(nparc);
nnetw = '17'; NumberOfNetworks = str2num(nnetw);
load(sprintf('../data/integration/Integration%s_%s_task-All.mat',nparc,nnetw));
load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));

% Change in Integration per network
for n = 1:NumberOfNetworks
    I17_change(n,:) = HI.Networks.Itot(n,:,2) - HI.Networks.Itot(n,:,1);
end
for p = 1:NumberOfParcels
    I17_changevec(p,:) = I17_change(Yeo_Shf(p),:);   
end

% Change in amplitude fluctuation per parcel
for z = 1:20
    for p = 1:NumberOfParcels
        STD_changevec(p,z) = STDDEV_surf_Parcels(2,z,p) - STDDEV_surf_Parcels(1,z,p);
    end
end

% Average per network
for n = 1:NumberOfNetworks
    STD_Netchange(n,:) = mean(STD_changevec(Yeo_Shf==n,:),1);   
end
for p = 1:NumberOfParcels
    STD_Netchagevec(p,:) = STD_Netchange(Yeo_Shf(p),:);   
end

% Correlate the two
for p = 1:NumberOfParcels
    [rval, pval] = corrcoef(I17_changevec(p,:),STD_Netchagevec(p,:));
    amp_int(p) = rval(2);
    amp_int_sig(p) = pval(2);
end

[~,~,~,amp_int_sig_corr] = fdr_bh(amp_int_sig);


parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
%nonwall = [2:50 51:101];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = amp_int;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
clim = [-1,1];
figure;
hotred = customcolormap([0 0.25 0.5 0.75 1],{'#660320','#DD7661','#F8F3F0','#4692C0','#093262'});
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.2 0.75 0.2 0.2], [0.2 0.59 0.2 0.2], ...
            1, 2, clim, hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.39 0.75 0.2 0.2], [0.39 0.59 0.2 0.2], ...
            1, 2, clim, hotred);
        
        
%% Figure 3.D Global Signal

load('../regressors/GS.mat')

% violin plots
figure;
col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.9290 0.6940 0.1250];
violin(GS.gsf(:,1:3),'mc',[],'medc','k','facecolor',col,'facealpha',1)
ylim([0 6])


%% Figure 3.E Global Signal vs integration
nnetw='7';
load(sprintf('../data/integration/Integration%s_%s_task-All.mat',nparc,nnetw));
out = array2table([GS.gsf(:,1),HI.Itot(:,1),GS.gsf(:,2),HI.Itot(:,2),GS.gsf(:,3),HI.Itot(:,3)]);
out.Properties.VariableNames = {'gsf_WR','I_WR','gsf_SD','I_SD','gsf_PRN','I_PRN'};
writetable(out,'../data/integration/gsf_Int.csv')

% (scatterplots were created in R studio)

%% Figure 3.F Global signal contributors

nparc = '400'; NumberOfParcels = str2num(nparc);
keyword = 'task-All';
nnetw = '17'; NumberOfNets = str2num(nnetw); 
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];
NumberOfSubjects = 20; 
load(sprintf('../labels/Yeo%s_Shf%s.mat',nnetw,nparc));



for t = 1:length(times)
    for z = 1:NumberOfSubjects
        sub = erase(SubjectName{z},'-');
        load(sprintf('../data/Parcels_regr_cat/%s/%s/%s_%s.mat', ...
                     (char(SubjectName{z})),(char(times(t,3))),(char(SubjectName{z})),keyword));
        for p = 1:NumberOfParcels
           [rr,pp] = corrcoef(Parcels(p,:), GS.gs.(char(times(t,3))).(sub));
           rval(t,z,p) = rr(2);
           pval(t,z,p) = pp(2);
        end
    end
end

global_contr = squeeze(mean(rval,2));


% Plot labels to surface

SP = SurfStatAvSurf({'../../labels/fsaverage5/surf/lh.pial', ...
                     '../../labels/fsaverage5/surf/rh.pial'});
% create colormaps
I = customcolormap([0 0.1 0.2 0.5 0.8 0.9 1], {'#E23603','#FF8D33', '#FFE633', '#bababa', '#33C7FF', '#336EFF', '#336EFF'});

% plot results
for t = 1:3
    vec = global_contr(t,:)';
    clim = [-0.5 0.5];
    lim = prctile(vec,[0 90]);
    parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
    nonwall = [2:201 203:402];
    toMap = zeros(1,length(unique(parc)));
    toMap(nonwall) = vec;
    OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
    OnSurfS = SurfStatSmooth(OnSurf, SP, 8);
    BoSurfStat_calibrate2Views(OnSurfS, SP, ...
                [0.25 1-(0.3*t) 0.3 0.3], [0.55 1-(0.3*t) 0.3 0.3], ...
                1, 2, clim, I); 
end