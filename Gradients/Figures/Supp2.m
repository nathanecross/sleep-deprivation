%
% Supplementary Figure 2
% 
% Change in task activations for A. PVT, B. 2-back, C. ANT

%% Load dependencies and setup
addpath('dependencies/surfstat');
addpath('dependencies/cbrewer');
addpath('dependencies/customcolormap');
addpath('gradients');
NumberOfParcels = 400;
load('labels/fsaverage5/ShfParcels/ShfLabels400_17.mat'); 
FS = SurfStatAvSurf({'../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/lh.pial', '../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf(:,1)' rh_labels_Shf(:,1)'];
nonwall = [2:201 203:402];
load('gradients/embedding.mat');
tasks = {'PVT','Nback', 'ANT'};
u(:,1) = 1:400;
u(:,2:4) = embedding(:,1:3,1)*-1;

%%  Obtain active areas in response to task onset
% all tasks
clear BL PRE POST
for t = 1:length(tasks)
    load(sprintf('data/activations_%s.mat',tasks{t}));
    BL.(char(sprintf('%s',tasks{t}))) = Activations.BL; 
    PRE.(char(sprintf('%s',tasks{t}))) = Activations.PRE; 
    POST.(char(sprintf('%s',tasks{t}))) = Activations.POST; 
    [~,p,~,stats] = ttest(BL.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(1,:) = stats.tstat;
    actis.(char(sprintf('group%sp',tasks{t})))(1,:) = p;
    [~,~,~,stats] =ttest(PRE.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(2,:) = stats.tstat;
    [~,~,~,stats] = ttest(POST.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(3,:) = stats.tstat;
    for p = 1:NumberOfParcels
        [~,pv,~,stats] = ttest(PRE.(char(sprintf('%s',tasks{t})))(:,p), BL.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,1,:) = [0 0 stats.sd stats.tstat pv];
        [~,pv,~,stats] = ttest(POST.(char(sprintf('%s',tasks{t})))(:,p), PRE.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,2,:) = [0 0 stats.sd stats.tstat pv];
        [~,pv,~,stats] = ttest(POST.(char(sprintf('%s',tasks{t})))(:,p), BL.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,3,:) = [0 0 stats.sd stats.tstat pv];
    end
end

%% Supplementry Figure 2
% 2A. PVT
for prep_figureSupp2A = 1
    [sigBL,~,~,adjp] = fdr_bh(actis.groupPVTp);
    clear stats table
    for p = 1:400 % for each node
        in_table = table([1:size(BL.PVT,1)]', BL.PVT(:,p), ...
            PRE.PVT(:,p), POST.PVT(:,p), ...
            'VariableNames', {'subject', 'state1', 'state2', 'state3'});

        rm = fitrm(in_table, 'state1-state3~1');
        ranovatbl = ranova(rm);
        stats(p,:) = table2array(ranovatbl(1,:));
    end
    rois = fdr_bh(stats(:,5)); % FDR corrected across all cortex at BL
    table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
    parcels_stats = [num2cell(table.Var1) table.Var2 num2cell(stats(:,4:5))];
    uncorr_rois = parcels_stats(stats(:,5)<0.05,:); % Uncorrected changes in activations
    p = [1:400; stats(:,5)']';
    % FDR corrected across activations at BL
    nonsigBL = ~sigBL;
    all = [1:400; stats(:,5)'];
    keep = all(:,logical(sigBL));
    holdd = all(:,nonsigBL); holdd(2,:) = 1;
    [h, crit_p, adj_ci_cvrg, keep(2,:)]=fdr_bh(keep(2,:));
    fdrp = array2table([keep holdd]');
    fdrp = sortrows(fdrp);
    fdrp = table2array(fdrp);
    vec = stats(:,4);       
end

for figureSupp2A = 1
    f = figure('units', 'centimeters', 'position', [0 0 50 20]);
    vec(p(:,2)>0.05)=0;
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    %a(3) = axes('position', [0.15 0.83 0.1 0.01]);
    pos = [0.01 0.8 0.1 0.1; 0.01 0.68 0.1 0.1; 0.07 0.68 0.1 0.1; 0.07 0.8 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 

    vec(fdrp(:,2)>0.05)=0; 
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    pos = [0.01 0.5 0.1 0.1; 0.01 0.38 0.1 0.1; 0.07 0.38 0.1 0.1; 0.07 0.5 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 
end


% Supp2B. Nback
for prep_figureSupp2A = 1
    [sigBL,~,~,adjp] = fdr_bh(actis.groupNbackp);
    clear stats table
    for p = 1:400 % for each node
        in_table = table([1:size(BL.Nback,1)]', BL.Nback(:,p), ...
            PRE.Nback(:,p), POST.Nback(:,p), ...
            'VariableNames', {'subject', 'state1', 'state2', 'state3'});

        rm = fitrm(in_table, 'state1-state3~1');
        ranovatbl = ranova(rm);
        stats(p,:) = table2array(ranovatbl(1,:));
    end
    rois = fdr_bh(stats(:,5)); % FDR corrected across all cortex at BL
    table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
    parcels_stats = [num2cell(table.Var1) table.Var2 num2cell(stats(:,4:5))];
    uncorr_rois = parcels_stats(stats(:,5)<0.05,:); % Uncorrected changes in activations
    p = [1:400; stats(:,5)']';
    % FDR corrected across activations at BL
    nonsigBL = ~sigBL;
    all = [1:400; stats(:,5)'];
    keep = all(:,logical(sigBL));
    holdd = all(:,nonsigBL); holdd(2,:) = 1;
    [h, crit_p, adj_ci_cvrg, keep(2,:)]=fdr_bh(keep(2,:));
    fdrp = array2table([keep holdd]');
    fdrp = sortrows(fdrp);
    fdrp = table2array(fdrp);
    vec = stats(:,4);       
end

for figureSupp2A = 1
    vec(p(:,2)>0.05)=0;
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    %a(3) = axes('position', [0.15 0.83 0.1 0.01]);
    pos = [0.2 0.8 0.1 0.1; 0.2 0.68 0.1 0.1; 0.26 0.68 0.1 0.1; 0.26 0.8 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 

    vec(fdrp(:,2)>0.05)=0; 
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    pos = [0.2 0.5 0.1 0.1; 0.2 0.38 0.1 0.1; 0.26 0.38 0.1 0.1; 0.26 0.5 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 
end

% Supp2C. ANT
for prep_figureSupp2A = 1
    [sigBL,~,~,adjp] = fdr_bh(actis.groupANTp);
    clear stats table
    for p = 1:400 % for each node
        in_table = table([1:size(BL.ANT,1)]', BL.ANT(:,p), ...
            PRE.ANT(:,p), POST.ANT(:,p), ...
            'VariableNames', {'subject', 'state1', 'state2', 'state3'});

        rm = fitrm(in_table, 'state1-state3~1');
        ranovatbl = ranova(rm);
        stats(p,:) = table2array(ranovatbl(1,:));
    end
    rois = fdr_bh(stats(:,5)); % FDR corrected across all cortex at BL
    table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
    parcels_stats = [num2cell(table.Var1) table.Var2 num2cell(stats(:,4:5))];
    uncorr_rois = parcels_stats(stats(:,5)<0.05,:); % Uncorrected changes in activations
    p = [1:400; stats(:,5)']';
    % FDR corrected across activations at BL
    nonsigBL = ~sigBL;
    all = [1:400; stats(:,5)'];
    keep = all(:,logical(sigBL));
    holdd = all(:,nonsigBL); holdd(2,:) = 1;
    [h, crit_p, adj_ci_cvrg, keep(2,:)]=fdr_bh(keep(2,:));
    fdrp = array2table([keep holdd]');
    fdrp = sortrows(fdrp);
    fdrp = table2array(fdrp);
    vec = stats(:,4);       
end

for figureSupp2A = 1
    vec(p(:,2)>0.05)=0;
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    %a(3) = axes('position', [0.15 0.83 0.1 0.01]);
    pos = [0.39 0.8 0.1 0.1; 0.39 0.68 0.1 0.1; 0.45 0.68 0.1 0.1; 0.45 0.8 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 

    vec(fdrp(:,2)>0.05)=0; 
    toMap = zeros(1,length(unique(parc)));
            toMap(nonwall) = vec;
            OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
    clim = [0 4];
    J = customcolormap([0 0.5 1], {'#3800AF', '#8E5AFC', '#bababa'});
    pos = [0.39 0.5 0.1 0.1; 0.39 0.38 0.1 0.1; 0.45 0.38 0.1 0.1; 0.45 0.5 0.1 0.1];
    BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, J); 
end

table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
Sig.ANT = [num2cell(table.Var1) table.Var2 num2cell(vec) num2cell(p(:,2))];
Sig.ANT = Sig.ANT(p(:,2)<0.05,:);