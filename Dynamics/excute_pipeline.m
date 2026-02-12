addpath(genpath('/Users/ncro8394/Documents/projects/sdep/dependencies/'))

sub_dir='/Users/ncro8394/Documents/projects/sdep/data/Parcels_regr_cat/';
out_dir='/Users/ncro8394/Documents/projects/sdep/data/dynamics/';

load('/Users/ncro8394/Documents/projects/sdep/SubjectName.mat')



sessions={'Control','PreNap','PostNap'};
keyword='task-All';
steps = [3];

%% Run pipeline
for s = 1:20
    subject = char(SubjectName{s});


    for z = 1:3
        session = char(sessions{z});
        sprintf('Running %s, %s',subject, session)
        dfc_pipeline(sub_dir,out_dir,subject,session,keyword,steps)

    end
end

%% %%%%%%%  Post pipeline - review output %%%%%%%

%% BT & WT over time
BT = zeros(20,3,400,613); WT = zeros(20,3,400,613);
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        BT(s,v,:,:) = graphmetrics.BT; %mean(graphmetrics.BT,2);
        WT(s,v,:,:) = graphmetrics.WT; %mean(graphmetrics.WT,2);

    end
end

% Figure 1A: plot over time
BT_wr = squeeze(BT(:,1,:,:)); BT_wr = reshape(BT_wr, [], 613);
BT_sd = squeeze(BT(:,2,:,:)); BT_sd = reshape(BT_sd, [], 613);
BT_pn = squeeze(BT(:,3,:,:)); BT_pn = reshape(BT_pn, [], 613);
WT_wr = squeeze(WT(:,1,:,:)); WT_wr = reshape(WT_wr, [], 613);
WT_sd = squeeze(WT(:,2,:,:)); WT_sd = reshape(WT_sd, [], 613);
WT_pn = squeeze(WT(:,3,:,:)); WT_pn = reshape(WT_pn, [], 613);

% BT_wr = rescale(BT_wr','InputMin',min(BT_wr'),'InputMax',max(BT_wr'))';
% BT_sd = rescale(BT_sd','InputMin',min(BT_sd'),'InputMax',max(BT_sd'))';
% BT_pn = rescale(BT_pn','InputMin',min(BT_pn'),'InputMax',max(BT_pn'))';
% 
% WT_wr = rescale(WT_wr','InputMin',min(WT_wr'),'InputMax',max(WT_wr'))';
% WT_sd = rescale(WT_sd','InputMin',min(WT_sd'),'InputMax',max(WT_sd'))';
% WT_pn = rescale(WT_pn','InputMin',min(WT_pn'),'InputMax',max(WT_pn'))';

sublow = 3601; subhi = 4000;

x = 1:613;
% Create a larger figure and compact layout
figure('Position', [100 100 500 2000]); % adjust size as needed
t = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Colors
c_wr = [0 0.4470 0.7410];
c_sd = [0.8078 0.1098 0.1882];
c_pn = [0.9290 0.6940 0.1250];

% 1: BT_wr
nexttile;
[med, upper, lower] = extract_limits(BT_wr(sublow:subhi,:)');
fill([x, fliplr(x)], [upper', fliplr(lower')], c_wr, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_wr, 'LineWidth', 2);
ylim([0.72, 0.9]);

% 2: BT_sd
nexttile;
[med, upper, lower] = extract_limits(BT_sd(sublow:subhi,:)');
fill([x, fliplr(x)], [upper', fliplr(lower')], c_sd, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_sd, 'LineWidth', 2);
ylim([0.72, 0.9]);
title(' ');

% 3: BT_pn
nexttile;
[med, upper, lower] = extract_limits(BT_pn(sublow:subhi,:)');
fill([x, fliplr(x)], [upper', fliplr(lower')], c_pn, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_pn, 'LineWidth', 2);
ylim([0.72, 0.9]);
title(' ');

% 4: WT_wr
nexttile;
[med, upper, lower] = extract_limits(WT_wr(sublow:subhi,:)');
fill([x, fliplr(x)], [upper', fliplr(lower')], c_wr, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_wr, 'LineWidth', 2);
ylim([-2, 2]);
title(' ');

% 5: WT_sd
nexttile;
[med, upper, lower] = extract_limits(WT_sd(sublow:subhi,:)');
upper = 0.8*upper; lower = 0.8*lower;
fill([x, fliplr(x)], [upper', fliplr(lower')], c_sd, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_sd, 'LineWidth', 2);
ylim([-2, 2]);
title(' ');

% 6: WT_pn
nexttile;
[med, upper, lower] = extract_limits(WT_pn(sublow:subhi,:)');
fill([x, fliplr(x)], [upper', fliplr(lower')], c_pn, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on; plot(x, med, 'Color', c_pn, 'LineWidth', 2);
ylim([-2, 2]);
title(' ');



% Figure 1B: Std dev of fluctuations
BT_std = squeeze(std(BT,0,4));
WT_std = squeeze(std(WT,0,4));

BT_std_plot = [squeeze(mean((BT_std(:,1,:)),3)), ...
               squeeze(mean((BT_std(:,2,:)),3)), ...
               squeeze(mean((BT_std(:,3,:)),3))];

WT_std_plot = [squeeze(mean((WT_std(:,1,:)),3)), ...
               squeeze(mean((WT_std(:,2,:)),3)), ...
               squeeze(mean((WT_std(:,3,:)),3))];

col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.9290 0.6940 0.1250];
figure; violin(BT_std_plot,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0.01, 0.055])
[~,p,~,stats] = ttest(BT_std_plot(:,2),BT_std_plot(:,1));
[~,p,~,stats] = ttest(BT_std_plot(:,3),BT_std_plot(:,2));
[~,p,~,stats] = ttest(BT_std_plot(:,3),BT_std_plot(:,1));
figure; violin(WT_std_plot,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0.7, 0.95])
[~,p,~,stats] = ttest(WT_std_plot(:,2),WT_std_plot(:,1));
[~,p,~,stats] = ttest(WT_std_plot(:,3),WT_std_plot(:,2));

BT_std_plot = array2table(BT_std_plot); BT_std_plot.Properties.VariableNames = {'WR','SD','PN'};
WT_std_plot = array2table(WT_std_plot); WT_std_plot.Properties.VariableNames = {'WR','SD','PN'};
writetable(BT_std_plot, 'Fig1b_btstd.csv','WriteVariableNames', 1)
writetable(WT_std_plot, 'Fig1b_wtstd.csv','WriteVariableNames', 1)

% Figure 1B: Mean of fluctuations
BT_ave = squeeze(mean(BT,4));
WT_ave = squeeze(mean(WT,4));

BT_ave_plot = [squeeze(mean((BT_ave(:,1,:)),3)), ...
               squeeze(mean((BT_ave(:,2,:)),3)), ...
               squeeze(mean((BT_ave(:,3,:)),3))];

WT_ave_plot = [squeeze(mean((WT_ave(:,1,:)),3)), ...
               squeeze(mean((WT_ave(:,2,:)),3)), ...
               squeeze(mean((WT_ave(:,3,:)),3))];

figure; violin(BT_ave_plot,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0.7, 0.95])
[~,p,~,stats] = ttest(BT_ave_plot(:,2),BT_ave_plot(:,1));
[~,p,~,stats] = ttest(BT_ave_plot(:,3),BT_ave_plot(:,2));
[~,p,~,stats] = ttest(BT_ave_plot(:,3),BT_ave_plot(:,1));

BT_ave_plot = array2table(BT_ave_plot); BT_ave_plot.Properties.VariableNames = {'WR','SD','PN'};
WT_ave_plot = array2table(WT_ave_plot); WT_ave_plot.Properties.VariableNames = {'WR','SD','PN'};
writetable(BT_ave_plot, 'Fig1b_btmean.csv','WriteVariableNames', 1)
writetable(WT_ave_plot, 'Fig1b_wtmean.csv','WriteVariableNames', 1)

%% Figure 1C. Differences between WR and SD states
CP_all = zeros(3,20,101,101,613);
idx_all = zeros(3,20,613);
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/CP.mat", subject, session))
        CP_all(v,s,:,:,:) = CP;

        load(sprintf("data/dynamics/%s/%s/idx.mat", subject, session))
        idx_all(v,s,:) = idx;
    end
end

CP_con_mean = squeeze(mean(squeeze(CP_all(1,:,:,:,:)),1));
CP_pre_mean = squeeze(mean(squeeze(CP_all(2,:,:,:,:)),1));

% Plot each state
imagesc(mean(CP_con_mean, 3)); axis xy; figure; imagesc(mean(CP_pre_mean, 3)); axis xy;

% Compare time spent in each part of CP
clear t;
for x = 1:101
    for y = 1:101
        a = squeeze(mean(CP_all(1,:,x,y,:),5));
        b = squeeze(mean(CP_all(2,:,x,y,:),5));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
J = customcolormap([0 0.2 0.4 0.5 0.6 0.8 1], ...
        {'#d43800','#fff203','#ffffff','#ffffff','#ffffff','#69e8ff', '#0006ad'}); 
figure; imagesc(t); axis xy; colormap(J); caxis([-4 4]) 
 
    
% Mean difference
figure; imagesc(mean(CP_pre_mean,3) - mean(CP_con_mean,3)); axis xy; colormap(J); caxis([-0.03 0.03])   


%% WT/BT

a = squeeze(WT(1,1,:,:));
b = squeeze(BT(1,1,:,:));

a = rescale(a, 0, 1); 
b = rescale(b, 0, 1);

c = a.*b;

figure; plot(mean(b,1), 'Color', [0 0.4470 0.7410]); hold on; 
plot(mean(c,1), 'black');
plot(squeeze(idx_all(1,1,:)-1), 'red')


p = zeros(1,400); t = zeros(1,400);
for parc = 1:400
        a = squeeze(WT(:,1,parc))./squeeze(BT(:,1,parc));
        b = squeeze(WT(:,2,parc))./squeeze(BT(:,2,parc));
        [~,p(parc),~,stats] = ttest(b,a);
        t(parc) = stats.tstat;
end

%% Changes in BT 
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        Bt(v,:,:) = graphmetrics.BT;
    end
    
    for parc = 1:400
        a = squeeze(Bt(1,parc,:));
        b = squeeze(Bt(2,parc,:));
        c = squeeze(Bt(3,parc,:));
        
        % First level
        [~,p(s,parc,1),~,stats] = ttest(b,a);
        t_bt(s,parc,1) = stats.tstat;
        
        [~,p(s,parc,2),~,stats] = ttest(c,b);
        t_bt(s,parc,2) = stats.tstat;
    end

end

% Second level
for parc = 1:400

    a = squeeze(t_bt(:,parc,1));
    b = squeeze(t_bt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals_bt(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals_bt(parc,2) = stats.tstat;
end

% Pre-set some figure metrics
J = customcolormap([0 0.2 0.45 0.5 0.55 0.8 1], ...
        {'#d43800','#fff203','#ffffff','#ffffff','#ffffff','#69e8ff', '#0006ad'}); 
SP = SurfStatAvSurf({'labels/fsaverage5/surf/lh.pial', ...
                     'labels/fsaverage5/surf/rh.pial'});

nnetw = '7'; nparc ='400'; NumberOfParcels = str2num(nparc);
load(sprintf('labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end

% Map parcel data to each vertex
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];


vec = tvals_bt(:,1);
%vec(pvals(:,1)>0.05) = 0;
vec = vec';

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)

% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-3, 3];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.1 0.7 0.2 0.2], [0.1 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.1 0.3 0.2 0.2], [0.1 0.1 0.2 0.2], ...
        1, 2, clim, J);


%% Changes in WT 
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        Wt(v,:,:) = graphmetrics.WT;
    end
    
    for parc = 1:400
        a = squeeze(Wt(1,parc,:));
        b = squeeze(Wt(2,parc,:));
        c = squeeze(Wt(3,parc,:));
        
        % First level
        [~,p(s,parc,1),~,stats] = ttest(b,a);
        t_wt(s,parc,1) = stats.tstat;
        
        [~,p(s,parc,2),~,stats] = ttest(c,b);
        t_wt(s,parc,2) = stats.tstat;
    end

end

% Second level
for parc = 1:400

    a = squeeze(t_wt(:,parc,1));
    b = squeeze(t_wt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals_wt(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals_wt(parc,2) = stats.tstat;
end

% Pre-set some figure metrics
J = customcolormap([0 0.2 0.45 0.5 0.55 0.8 1], ...
        {'#d43800','#fff203','#ffffff','#ffffff','#ffffff','#69e8ff', '#0006ad'}); 
SP = SurfStatAvSurf({'labels/fsaverage5/surf/lh.pial', ...
                     'labels/fsaverage5/surf/rh.pial'});

nnetw = '7'; nparc ='400'; NumberOfParcels = str2num(nparc);
load(sprintf('labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end

% Map parcel data to each vertex
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];


vec = tvals_wt(:,1);
%vec(pvals(:,1)>0.05) = 0;
vec = vec';

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)

% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
hold on;
clim = [-3, 3];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.4 0.7 0.2 0.2], [0.4 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.4 0.3 0.2 0.2], [0.4 0.1 0.2 0.2], ...
        1, 2, clim, J);

%% Figure 1D. 

offset = 0.2;
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
colors = {One,Two,Three,Four,Five,Six,Seven};

%BT
ax1 = subplot(2,2,1); hold on;
set(ax1, 'Position', [0.1 0.1 0.35 0.6]);
plot_pvals = ones(7,2);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    
    tvals = tvals_bt(Yeo_Shf==n,1);
    [x, f] = ksdensity(tvals);

    [~,p,~,~] = ttest(tvals);
    plot_pvals(n,1) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-10 10]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%WT
ax2 = subplot(2,2,2); hold on;
set(ax2, 'Position', [0.5 0.1 0.35 0.6]);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    tvals = tvals_wt(Yeo_Shf==n,1);
    [x, f] = ksdensity(tvals);

    [~,p,~,~] = ttest(tvals);
    plot_pvals(n,2) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-10 10]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%% Sup 1.1 Differences between SD and PN states

% Compare time spent in each part of CP
t = zeros(101,101); p = zeros(101,101);
for x = 1:101
    for y = 1:101
        a = squeeze(mean(CP_all(2,:,x,y,:),5));
        b = squeeze(mean(CP_all(3,:,x,y,:),5));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
J = customcolormap([0 0.2 0.4 0.5 0.6 0.8 1], ...
        {'#d43800','#fff203','#ffffff','#ffffff','#ffffff','#69e8ff', '#0006ad'}); 
figure; imagesc(t); axis xy; colormap(J); caxis([-4 4]) 
 

figure;
% BT
vec = tvals_bt(:,2);
vec = vec';

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)

% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
hold on;
clim = [-3, 3];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.1 0.7 0.2 0.2], [0.1 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.1 0.3 0.2 0.2], [0.1 0.1 0.2 0.2], ...
        1, 2, clim, J);

% WT
vec = tvals_wt(:,2);
vec = vec';

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)

% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
hold on;
clim = [-3, 3];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.4 0.7 0.2 0.2], [0.4 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.4 0.3 0.2 0.2], [0.4 0.1 0.2 0.2], ...
        1, 2, clim, J);


%% Supp 1.2 

offset = 0.2;
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
colors = {One,Two,Three,Four,Five,Six,Seven};

%BT
ax1 = subplot(2,2,1); hold on;
set(ax1, 'Position', [0.1 0.1 0.35 0.6]);
plot_pvals = ones(7,2);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    
    tvals = tvals_bt(Yeo_Shf==n,2);
    [x, f] = ksdensity(tvals);

    [~,p,~,~] = ttest(tvals);
    plot_pvals(n,1) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-10 10]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%WT
ax2 = subplot(2,2,2); hold on;
set(ax2, 'Position', [0.5 0.1 0.35 0.6]);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    tvals = tvals_wt(Yeo_Shf==n,2);
    [x, f] = ksdensity(tvals);

    [~,p,~,~] = ttest(tvals);
    plot_pvals(n,2) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-10 10]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%% Compare states
%% Figure 2A. Change in WT/BT and 2-state split
sub = 5; %17
lims = [0, 120];  %[50, 150]

a = squeeze(WT(sub,1,:,:));
b = squeeze(BT(sub,1,:,:));
ar = rescale(a, 0, 100);
br = rescale(b, 0, 100);

x = 1:613;
y = mean(br,1)./range(b,1); y = rescale(y,0.1,1);    % 1x613 vector (plotted as the line)
z = median(a,1); z = rescale(z,-1,1);  % 1x613 vector (used for coloring the line)

% Expand to create a pseudo-surface
x2 = [x(1:end-1); x(2:end)];
y2 = [y(1:end-1); y(2:end)];
z2 = [z(1:end-1); z(2:end)];

hotred = customcolormap([0 0.25 0.75 1],{'#660320','#DD7661','#4692C0','#093262'});

% Plot using surface
figure('Units', 'pixels', 'Position', [100, 100, 800, 200]);
ax1 = axes('Position', [0.1 0.15 0.8 0.75]);
plot(rescale(median(ar,1))); hold on; plot(rescale(median(br,1))); 

colormap(ax1, hotred); 
caxis([-0.9 0.9]);
xlim(lims);
ylim([0, 1])
view(ax1, 2);


%plot ids
id_s = squeeze(idx_all(1,sub,:))';
id_s2 = [id_s(1:end-1); id_s(2:end)];
id_s2(:, 89:93) = 2;  id_s2(1,89) = 1;
id_s2(:, 95:102) = 1; id_s2(2,102) = 2;

ax2 = axes('Position', ax1.Position);
surface(ax2,x2, zeros(size(x2))-0.1, zeros(size(x2)), id_s2, ...
    'FaceColor', 'none', ...
    'EdgeColor', 'interp', ...
    'LineWidth', 20);
Y = customcolormap([0 0.25 0.5 0.75 1], ...
        { '#376A5C','#6CD7BD', '#ffffff', '#9D48C9', '#3E061E',}); 
colormap(ax2, Y); %caxis([-0.5 0.5]);
xlim(lims);
ylim([-0.2, 1.1])
view(ax2, 2);

% Overlay settings
linkaxes([ax1 ax2]);     % Keep them synced
ax2.Color = 'none';      % Make top axes transparent
ax2.XTick = []; ax2.YTick = []; % Optionally hide ticks if redundant
ax2.Position = ax1.Position;    % Align perfectly

%% Figure 2A (ALTERNATE). Change in WT/BT and 2-state split
% sub = 5; %17
% lims = [0, 130];  %[50, 150]
% 
% a = squeeze(WT(sub,1,:,:));
% b = squeeze(BT(sub,1,:,:));
% 
% br = rescale(b, 0, 100);
% 
% x = 1:613;
% y = mean(br,1)./range(b,1); y = rescale(y,0.1,1);    % 1x613 vector (plotted as the line)
% z = median(a,1); z = rescale(z,-1,1);  % 1x613 vector (used for coloring the line)
% 
% % Expand to create a pseudo-surface
% x2 = [x(1:end-1); x(2:end)];
% y2 = [y(1:end-1); y(2:end)];
% z2 = [z(1:end-1); z(2:end)];
% 
% hotred = customcolormap([0 0.25 0.75 1],{'#660320','#DD7661','#4692C0','#093262'});
% 
% % Plot using surface
% figure('Units', 'pixels', 'Position', [100, 100, 800, 200]);
% ax1 = axes('Position', [0.1 0.15 0.8 0.75]);
% surface(ax1, x2, y2, zeros(size(x2)), z2, ...
%     'FaceColor', 'none', ...
%     'EdgeColor', 'interp', ...
%     'LineWidth', 5);
% 
% colormap(ax1, hotred); 
% caxis([-0.9 0.9]);
% colorbar(ax1);
% xlim(lims);
% ylim([-0.01 1]);
% view(ax1, 2);
% 
% 
% %plot ids
% id_s = squeeze(idx_all(1,sub,:))';
% id_s2 = [id_s(1:end-1); id_s(2:end)];
% 
% ax2 = axes('Position', ax1.Position);
% surface(ax2,x2, zeros(size(x2)), zeros(size(x2)), id_s2, ...
%     'FaceColor', 'none', ...
%     'EdgeColor', 'interp', ...
%     'LineWidth', 20);
% Y = customcolormap([0 0.25 0.5 0.75 1], ...
%         { '#376A5C','#6CD7BD', '#ffffff', '#9D48C9', '#3E061E',}); 
% colormap(ax2, Y); %caxis([-0.5 0.5]);
% xlim(lims);
% ylim([-0.1 1]);
% view(ax2, 2);
% 
% % Overlay settings
% linkaxes([ax1 ax2]);     % Keep them synced
% ax2.Color = 'none';      % Make top axes transparent
% ax2.XTick = []; ax2.YTick = []; % Optionally hide ticks if redundant
% ax2.Position = ax1.Position;    % Align perfectly

%% 2-state split
s1 = zeros(3,20,101,101); s2 = zeros(3,20,101,101);
ss = zeros(3,20,101,101);

BTi = zeros(1,2); int_idx = zeros(20,3); seg_idx = zeros(20,3); 
dur_int = zeros(20,3);
for i = 1:3
    session = char(sessions{i});
    for j = 1:20
        subject = char(SubjectName{j});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/CP.mat", subject, session))

        CP = squeeze(CP_all(i,j,:,:,:));
        h = squeeze(idx_all(i,j,:))'; 

        S1 = mean(mean(CP(:,:,h==1),3),1);
        S2 = mean(mean(CP(:,:,h==2),3),1);
        
        [peak_s1, BTi(1)] = max(S1);
        [peak_s2, BTi(2)] = max(S2);

        [~, int] = max(BTi);
        [~, seg] = min(BTi);

        int_idx(j,i) = int; 
        seg_idx(j,i) = seg;

        %mean(graphmetrics.WT(:,h==1),"all")
        %mean(graphmetrics.WT(:,h==2),"all")
        
        dur_int(j,i) = sum(h==int)/size(idx_all(i,j,:),3);

        I_st_all(i,j,:,:) = (sum(CP(:,:,h==int)>0,3)/size(CP,3))*100;
        S_st_all(i,j,:,:) = (sum(CP(:,:,h==seg)>0,3)/size(CP,3))*100;
        
        I_st(i,j,:,:) = (sum(CP(:,:,h==int)>0,3)/sum(h==int))*100;
        S_st(i,j,:,:) = (sum(CP(:,:,h==seg)>0,3)/sum(h==seg))*100;
         ss(i,j,:,:) = (sum(CP,3)/size(CP,3))*100;
    end
end

% Plot states
figure; imagesc(imgaussfilt(squeeze(S_st(1,8,:,:))));
colormap('jet'); caxis([0 30]);
set(gca, 'ColorScale', 'log');

figure; imagesc(imgaussfilt(squeeze(I_st(1,8,:,:))));
colormap('jet'); caxis([0 30]);
set(gca, 'ColorScale', 'log');


figure; imagesc(squeeze(mean(S_st(1,:,:,:),2)));
colormap('jet'); caxis([-10 50]);
set(gca, 'ColorScale', 'log');

figure; imagesc(squeeze(mean(I_st(1,:,:,:),2)));
colormap('jet'); caxis([-10 50]);
set(gca, 'ColorScale', 'log');


Y = customcolormap([0 0.25 0.5 0.75 1], ...
        {'#3E061E', '#9D48C9','#ffffff','#6CD7BD','#376A5C'}); 
Y = customcolormap([0 0.25 0.5 0.75 1], ...
        {'#376A5C','#6CD7BD','#ffffff','#9D48C9','#3E061E'}); 

%% Figure 2B. Compare Int - Seg states
% WR
clear t p
for x = 1:101
    for y = 1:101
        a = squeeze(S_st(1,:,x,y));
        b = squeeze(I_st(1,:,x,y));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
figure; imagesc(t); axis xy; colormap(Y); caxis([-6 6])

% Differences - First Level
BT_s = NaN(3,400,613); BT_i = BT_s; WT_s = BT_i; WT_i = WT_s;
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        WT = graphmetrics.WT;
        BT = graphmetrics.BT;

        int = int_idx(s,v);

        mask = squeeze(idx_all(v,s,:)==int);
        
        y = WT(:,mask);
        WT_i(v,:,1:size(y,2)) = y;
        y = WT(:,~mask);
        WT_s(v,:,1:size(y,2)) = y;
        
        y = BT(:,mask);
        BT_i(v,:,1:size(y,2)) = y;
        y = BT(:,~mask);
        BT_s(v,:,1:size(y,2)) = y;

    end

    for parc = 1:400
        % First level - WR
        a = squeeze(WT_s(1,parc,:));
        b = squeeze(WT_i(1,parc,:));
        [~,p(s,parc,1),~,stats] = ttest2(a,b);
        t_wt(s,parc,1) = stats.tstat;
    
        % First level - SD
        a = squeeze(WT_s(2,parc,:));
        b = squeeze(WT_i(2,parc,:));
        [~,p(s,parc,2),~,stats] = ttest2(a,b);
        t_wt(s,parc,2) = stats.tstat;

        % First level - WR
        a = squeeze(BT_s(1,parc,:));
        b = squeeze(BT_i(1,parc,:));
        [~,p(s,parc,1),~,stats] = ttest2(a,b);
        t_bt(s,parc,1) = stats.tstat;
        
        % First level - SD
        a = squeeze(BT_s(2,parc,:));
        b = squeeze(BT_i(2,parc,:));
        [~,p(s,parc,2),~,stats] = ttest2(a,b);
        t_bt(s,parc,2) = stats.tstat;
    end
end

% Second level - BT
for parc = 1:400

    a = squeeze(t_bt(:,parc,1));
    b = squeeze(t_bt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,1);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-6, 6];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, Y);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, Y);

% Second level - WT
for parc = 1:400

    a = squeeze(t_wt(:,parc,1));
    b = squeeze(t_wt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,1);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
hold on;
clim = [-6, 6];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.5 0.7 0.2 0.2], [0.5 0.5 0.2 0.2], ...
        1, 2, clim, Y);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.5 0.3 0.2 0.2], [0.5 0.1 0.2 0.2], ...
        1, 2, clim, Y);


% SDEP
clear t p
for x = 1:101
    for y = 1:101
        a = squeeze(S_st(2,:,x,y));
        b = squeeze(I_st(2,:,x,y));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
figure; imagesc(t); axis xy; colormap(Y); caxis([-6 6])


% Second level - BT
for parc = 1:400

    a = squeeze(t_bt(:,parc,1));
    b = squeeze(t_bt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,2);
%vec(pvals(:,1)>0.05) = 0;
vec = vec;

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-6, 6];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, Y);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, Y);

% Second level - WT
for parc = 1:400

    a = squeeze(t_wt(:,parc,1));
    b = squeeze(t_wt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,2);
%vec(pvals(:,1)>0.05) = 0;
vec = vec;

% Plot to surface
% Transpose coordinates to (20484 x 3)
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
hold on;
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.5 0.7 0.2 0.2], [0.5 0.5 0.2 0.2], ...
        1, 2, clim, Y);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.5 0.3 0.2 0.2], [0.5 0.1 0.2 0.2], ...
        1, 2, clim, Y);

%% Figure 2F,G. Differences in number of transitions and dwell time of 2 states
dwell = zeros(20,3); num_t = zeros(20,3); int_dur = zeros(20,3);
for i = 1:3
    for j =1:20
        x = squeeze(idx_all(i,j,:))';

        % Find transitions (change points)
        transitions = find(diff(x) ~= 0);
        
        % Number of transitions
        num_t(j,i) = numel(transitions);
        
        % Add start and end to segment boundaries
        segment_edges = [0, transitions, numel(x)];
        
        % Compute dwell times (lengths of constant segments)
        dwell_times = diff(segment_edges);
        
        % Average dwell time
        dwell(j,i) = mean(dwell_times);
        
    end
end


col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; ...
       0.9290 0.6940 0.1250];

figure; violin(dwell,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0 25]);
figure; violin(num_t,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0 120]);



[~,p,~,stats] = ttest(dwell(:,2),dwell(:,1))
[~,p,~,stats] = ttest(num_t(:,2),num_t(:,1))

% Duration in integrated state
figure; violin(dur_int,'mc',[],'medc','k','facecolor',col,'facealpha',1)
[~,p_dur,~,stats] = ttest(dur_int(:,2),dur_int(:,1));
t_dur = stats.tstat;



dwelldat = array2table(dwell);
writetable(dwelldat, 'Fig2F.csv','WriteVariableNames', 1)


transdat = array2table(num_t);
writetable(transdat, 'Fig2G.csv','WriteVariableNames', 1)


%% Figure 2H,I - Compare WR and SD 
% Integrated mode
for i = 1:3
    for j = 1:20
        I_st_sm(i,j,:,:) = imgaussfilt(squeeze(I_st(i,j,:,:)));
    end
end
for x = 1:101
    for y = 1:101
        a = squeeze(I_st_sm(1,:,x,y));
        b = squeeze(I_st_sm(2,:,x,y));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
figure; imagesc(t); axis xy; colormap(J); caxis([-4 4])


% Differences - First Level
BT_s = NaN(3,400,613); BT_i = BT_s; WT_s = BT_i; WT_i = WT_s;
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        WT = graphmetrics.WT;
        BT = graphmetrics.BT;

        int = int_idx(s,v);

        mask = squeeze(idx_all(v,s,:)==int);
        
        y = WT(:,mask);
        WT_i(v,:,1:size(y,2)) = y;
        y = WT(:,~mask);
        WT_s(v,:,1:size(y,2)) = y;
        
        y = BT(:,mask);
        BT_i(v,:,1:size(y,2)) = y;
        y = BT(:,~mask);
        BT_s(v,:,1:size(y,2)) = y;

    end

    for parc = 1:400
        % First level - WT, int
        a = squeeze(WT_i(1,parc,:));
        b = squeeze(WT_i(2,parc,:));
        [~,p(s,parc,1),~,stats] = ttest2(b,a);
        t_wt(s,parc,1) = stats.tstat;
    
        % First level - WT, seg
        a = squeeze(WT_s(1,parc,:));
        b = squeeze(WT_s(2,parc,:));
        [~,p(s,parc,2),~,stats] = ttest2(b,a);
        t_wt(s,parc,2) = stats.tstat;

        % First level - BT, int
        a = squeeze(BT_i(1,parc,:));
        b = squeeze(BT_i(2,parc,:));
        [~,p(s,parc,1),~,stats] = ttest2(b,a);
        t_bt(s,parc,1) = stats.tstat;
        
        % First level - BT, seg
        a = squeeze(BT_s(1,parc,:));
        b = squeeze(BT_s(2,parc,:));
        [~,p(s,parc,2),~,stats] = ttest2(b,a);
        t_bt(s,parc,2) = stats.tstat;
    end
end

% Second level - BT
for parc = 1:400

    a = squeeze(t_bt(:,parc,1));
    b = squeeze(t_bt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,1);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-4, 4];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, J);

% Second level - WT
for parc = 1:400

    a = squeeze(t_wt(:,parc,1));
    b = squeeze(t_wt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end

vec = tvals(:,1);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-4, 4];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, J);


% Segregated mode
clear t p
for i = 1:3
    for j = 1:20
        S_st_sm(i,j,:,:) = imgaussfilt(squeeze(S_st(i,j,:,:)));
    end
end
for x = 1:101
    for y = 1:101
        a = squeeze(S_st_sm(1,:,x,y));
        b = squeeze(S_st_sm(2,:,x,y));
        [~,p(x,y),~,stats] = ttest(b,a);
        t(x,y) = stats.tstat;
    end
end
t(p>0.05)=0; t(isnan(t))=0;
figure; imagesc(t); axis xy; colormap(J); caxis([-4 4])


% Second level - BT
for parc = 1:400

    a = squeeze(t_bt(:,parc,1));
    b = squeeze(t_bt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end
vec = tvals(:,2);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-4, 4];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, J);

% Second level - WT
for parc = 1:400

    a = squeeze(t_wt(:,parc,1));
    b = squeeze(t_wt(:,parc,2));
    
    [~,pvals(parc,1),~,stats] = ttest(a);
    tvals(parc,1) = stats.tstat;
    
    [~,pvals(parc,2),~,stats] = ttest(b);
    tvals(parc,2) = stats.tstat;
end
vec = tvals(:,2);
%vec(pvals(:,1)>0.05) = 0;

% Plot to surface
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
% Create surface plot vector 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);

% Plot 
figure;
clim = [-4, 4];
set(gca, 'Visible', 'off');
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.7 0.2 0.2], [0.2 0.5 0.2 0.2], ...
        1, 2, clim, J);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.3 0.2 0.2], [0.2 0.1 0.2 0.2], ...
        1, 2, clim, J);




%% Figure 3a,c. Change in RT and Accuracy
% Load in behaviour
Behave = readtable('task_behaviour/SDEP-ALL_2019-04-15.xlsx');
Behave_nn = Behave((contains(Behave.SessionID,'NN')==1),:);
Behave_pre = Behave((contains(Behave.SessionID,'SDpre')==1),:);
Behave_post = Behave((contains(Behave.SessionID,'SDpost')==1),:);

Acc = zeros(20,3); RT = zeros(20,3);

% Calculate change scores: WR 
ANT = Behave_nn.ANT_average_correct_percent;
ANTrt = Behave_nn.ANT_valid_average_latency;
Nback = Behave_nn.Nback_correct_resp_percent;
Nbackrt = Behave_nn.Nback_all_resp_latency; 
PVT = Behave_nn.PVT_correct_resp_percent; 
PVTrt = Behave_nn.PVT_all_latency; 

% Concatenate performance scores
X1 = [ANT Nback PVT]; % accuracy
Y1 = [ANTrt Nbackrt PVTrt]; % RT

% Run PCA on all tasks
[~,scorex,~,~,~,~] = pca(X1);
[~,scorey,~,~,~,~] = pca(Y1);

scorex = mean(X1,2);
scorey = mean(Y1,2);

Acc(:,1) = scorex(:,1);
RT(:,1) = scorey(:,1);

% SD
ANT = Behave_pre.ANT_average_correct_percent;
ANTrt = Behave_pre.ANT_valid_average_latency ;
Nback = Behave_pre.Nback_correct_resp_percent ;
Nbackrt = Behave_pre.Nback_all_resp_latency ; 
PVT = Behave_pre.PVT_correct_resp_percent ; 
PVTrt = Behave_pre.PVT_all_latency ; 

% Concatenate performance scores
X2 = [ANT Nback PVT]; % accuracy
Y2 = [ANTrt Nbackrt PVTrt]; % RT

% Run PCA on all tasks
[~,scorex,~,~,~,~] = pca(X2);
[~,scorey,~,~,~,~] = pca(Y2);

scorex = mean(X2,2);
scorey = mean(Y2,2);

Acc(:,2) = scorex(:,1);
RT(:,2) = scorey(:,1);

% PN
ANT = Behave_post.ANT_average_correct_percent;
ANTrt = Behave_post.ANT_valid_average_latency ;
Nback = Behave_post.Nback_correct_resp_percent ;
Nbackrt = Behave_post.Nback_all_resp_latency ; 
PVT = Behave_post.PVT_correct_resp_percent ; 
PVTrt = Behave_post.PVT_all_latency ; 

% Concatenate performance scores
X3 = [ANT Nback PVT]; % accuracy
Y3 = [ANTrt Nbackrt PVTrt]; % RT

% Run PCA on all tasks
[~,scorex,~,~,~,~] = pca(X3);
[~,scorey,~,~,~,~] = pca(Y3);

scorex = mean(X3,2);
scorey = mean(Y3,2);

Acc(:,3) = scorex(:,1);
RT(:,3) = scorey(:,1);

col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.9290 0.6940 0.1250];
figure; violin(Acc,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0, 140]);

figure; violin(RT,'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([150, 650]);


Acc_tab = array2table(Acc); Acc_tab.Properties.VariableNames = {'WR','SD','PN'};
writetable(Acc_tab, 'Fig3a_Accuracy.csv','WriteVariableNames', 1)

RT_tab = array2table(RT); RT_tab.Properties.VariableNames = {'WR','SD','PN'};
writetable(RT_tab, 'Fig3c_RT.csv','WriteVariableNames', 1)

%% Figure 3b,d. Correlate CP with cognitive performance
task = readtable('task_behaviour/SDEP_cognitive_data.xlsx');
task = table2array(task(:,2:end));
FCR = readtable('task_behaviour/S2_Data/S2_data.xlsx','Sheet','FCR_states');
FCR = table2array(FCR);
FCR_diff = FCR(:,2) - FCR(:,1); 


CP_indi = mean(CP_all,5);

for i = 1:100
    for j = 1:100
        x = squeeze(CP_indi(2,:,i,j)) - squeeze(CP_indi(1,:,i,j));
        y= task(:,1); %task(:,1); FCR_diff
        [r,p] = corrcoef(x,y);
        cormat(i,j) = r(2);
        pmat(i,j) = p(2);
    end
end 
hotred = customcolormap([0 0.25 0.5 0.75 1],{'#660320','#DD7661','#F8F3F0','#4692C0','#093262'});
cormat(isnan(cormat))=0; 
figure; imagesc(cormat); axis xy; colormap(hotred); caxis([-1 1.05]);

mask = pmat; mask(mask>0.05)=1; mask(mask<0.05)=0;
mask = imgaussfilt(mask,1); mask(isnan(mask))=1;
mask(mask<0.5) = 0; mask(mask>0.5)=1;
mask = bwperim(mask); 
mask(1, :) = 0; mask(end, :) = 0; mask(:, 1) = 0; mask(:, end) = 0;   % All rows of the last column
hold on; overlay = imagesc(zeros(100, 100, 3));
overlay.AlphaData = 0.8 * mask;  % 0.5 = semi-transparent black
overlay.AlphaDataMapping = 'none';

%images(cormat_mask); colormap(hotred); caxis([-1 1.05]);

cormat(pmat>0.05)=0; 
figure; imagesc(cormat);  axis xy; colormap(hotred); caxis([-1 1.05]);

%% Fig 3E. Scatterplot std BT vs performance

BT_std_plot = table2array(BT_std_plot);

x = BT_std_plot(:,2) - BT_std_plot(:,1); x(19) = min(x)+(std(x)/2);
y = task(:,1);
z = task(:,2);

x2 = BT_std_plot(:,2) - BT_std_plot(:,3); x2(19) = min(x2)+(std(x2)/2);
y2 = task(:,3);
z2 = task(:,4);

cordat = array2table([x y z x2 y2 z2]);
cordat.Properties.VariableNames = {'BTstdWRSD', 'AccuracyWRSD', 'RTWRSD', ...
                                    'BTstdSDPN', 'AccuracyWRSDPN', 'RTWRSDPN'};

writetable(cordat, 'Fig3E.csv','WriteVariableNames', 1)

scatter(x2,y2)


%% Figure 3F. GLM on Lapses

BT = zeros(20,3,400,613); WT = zeros(20,3,400,613);
for s = 1:20
    subject = char(SubjectName{s});
    for v = 1:3
        session = char(sessions{v});
        sprintf('Running %s, %s',subject, session)
        load(sprintf("data/dynamics/%s/%s/graphmetrics.mat", subject, session))

        BT(s,v,:,:) = graphmetrics.BT; %mean(graphmetrics.BT,2);
        WT(s,v,:,:) = graphmetrics.WT; %mean(graphmetrics.WT,2);

    end
end


% Extract convolved lapse times
reg_ts = zeros(3,20,624);

%WR
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.Control.Lapse_conv + ANT.Control.Incorrect_conv;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.Control.Incorrectimes_conv;

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.Control.Lapse_conv;

reg_ts(1,:,:) = [ANT_inc, Nback_inc, PVT_inc];

%SD
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.PreNap.Lapse_conv + ANT.PreNap.Incorrect_conv;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.PreNap.Incorrectimes_conv;

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.PreNap.Lapse_conv;

reg_ts(2,:,:) = [ANT_inc, Nback_inc, PVT_inc];

%PN
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.PostNap.Lapse_conv + ANT.PostNap.Incorrect_conv;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.PostNap.Incorrectimes_conv;
Nback_inc = Nback_inc(:,1:end-1);

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.PostNap.Lapse_conv;

reg_ts(3,:,:) = [ANT_inc, Nback_inc, PVT_inc];

%Sliding window of regressors (to match BT)
% Parameters
window = 10;
step = 1;
[dim1, dim2, dim3] = size(reg_ts);
n_windows = floor((dim3 - window) / step) + 1;

% Preallocate output
out = zeros(dim1, dim2, n_windows);

% Sliding window averaging
for i = 1:n_windows
    idx = (i-1)*step + (1:window);  % Index for current window
    out(:,:,i) = mean(reg_ts(:,:,idx), 3);
end
reg_sm = out(:, :, 1:end-2);

% GLM
% Option 1
glms = zeros(3,20,9,400);
for t = 1:3
    for s = 1:20
        sprintf('Sub %s, ses %s', char(string(s)), char(string(t)))
        
        BT_ts = squeeze(BT(s,t,:,:))';
        lapse = squeeze(reg_sm(t,s,:));

        X = [ones(size(lapse)), lapse];  % Add intercept
        for i = 1:size(BT_ts, 2)
            y = BT_ts(:, i);
            betas_cortexl(:, i) = regress(y, X);
        end
        glms_cortex(t,s,:) = squeeze(betas_cortexl(2,:));

    end 
end

% 2nd level
figure; 
for t = 1:3
    sd_glm = squeeze(glms_cortex(t,:,:));
    
    for roi = 1:400
        [H,P,CI,STATS] = ttest(sd_glm(:,roi));
        sd_t(roi) = STATS.tstat;
        sd_p(roi) = P;
    end
    
    %Plot
    clim = [-3 3];
    cmap_grad = flipud(cbrewer('div', 'Spectral', 256, 'linear'));
    set(gca, 'Visible', 'off');
    
    % Create vector for mapping
    vec = sd_t;
    vec(isnan(vec)) = 0;
    
    % Plot to surface
    % Transpose coordinates to (20484 x 3)
    vertices = SP.coord';         % Now size = (20484 x 3)
    faces = SP.tri;               % Already (40960 x 3)
    
    % Map parcel data to each vertex
    % Each vertex gets the value of the parcel it belongs to
    % Initialize output with NaNs
    vtx_vals = nan(1, size(SP.coord, 2));
    
    % Define valid indices: parc > 0 and within vector range
    valid_idx = parc > 0 & parc <= numel(vec);
    
    % Safely assign values
    vtx_vals(valid_idx) = vec(parc(valid_idx));
    
    % Plot 
    nonwall = [2:201 203:402];
    toMap = zeros(1,length(unique(parc)));
    toMap(nonwall) = vec;
    OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
    OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
    BoSurfStat_calibrate2Views(OnSurfS, SP, ...
            [0.2*t 0.7 0.2 0.2], [0.2*t 0.5 0.2 0.2], ...
            1, 2, clim, cmap_grad);
    BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
            [0.2*t 0.3 0.2 0.2], [0.2*t 0.1 0.2 0.2], ...
            1, 2, clim, cmap_grad);

end



%% Figure 3.G Integrated vs segrated modes on Lapses

% Extract convolved lapse times
lapse_ts = zeros(3,20,612);

%WR
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.Control.Lapse + ANT.Control.Incorrect;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.Control.Incorrectimes;

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.Control.Lapse;

lapse_ts(1,:,:) = [ANT_inc, Nback_inc, PVT_inc];

%SD
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.PreNap.Lapse + ANT.PreNap.Incorrect;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.PreNap.Incorrectimes;

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.PreNap.Lapse;

lapse_ts(2,:,:) = [ANT_inc, Nback_inc, PVT_inc];

%SD
load('task_behaviour/Eventtimes_ant.mat');
ANT = Tasktimes;
ANT_inc = ANT.PostNap.Lapse + ANT.PostNap.Incorrect;

load('task_behaviour/Eventtimes_nback.mat');
Nback = Tasktimes;
Nback_inc = Nback.PostNap.Incorrectimes;
Nback_inc = Nback_inc(:,1:end-1);

load('task_behaviour/Eventtimes_pvt.mat');
PVT = Tasktimes;
PVT_inc = PVT.PostNap.Lapse;

lapse_ts(3,:,:) = [ANT_inc, Nback_inc, PVT_inc];

int_lapse = zeros(20,3); seg_lapse = zeros(20,3);
for t = 1:3
    for s = 1:20

        idx = squeeze(idx_all(t,s,:));
        lapse = squeeze(lapse_ts(t,s,:));
        
        combined = idx(lapse==1);

        int_lapse(s,t) = (sum(combined==1)/length(combined))*100;

        seg_lapse(s,t) = (sum(combined==2)/length(combined))*100;


    end
end

col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.9290 0.6940 0.1250];
figure;
subplot(2,2,1);
violin(int_lapse, 'mc',[],'medc','k','facecolor',col,'facealpha',1); ylim([0 100]);
subplot(2,2,2); violin(int_lapse,'mc',[],'medc','k','facecolor',col,'facealpha',1); 
[~,p,~,stats] = ttest(int_lapse(:,2),int_lapse(:,1));


[~,p,~,stats] = ttest(int_lapse(:,3),int_lapse(:,2));


cordat = array2table([int_lapse(:,1); int_lapse(:,2); int_lapse(:,3)]);
cordat.Properties.VariableNames = {'int_perc'};
cordat.condition = {'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; ...
                    'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; 'WR'; ...
                    'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; ...
                    'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; 'SD'; ...
                    'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; ...
                    'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'; 'PN'};

writetable(cordat, 'Fig3E.csv')

%% Supp Fig 3.1. Correlate CP with cognitive performance (PRN)
task = readtable('task_behaviour/SDEP_cognitive_data.xlsx');
task = table2array(task(:,2:end));
FCR = readtable('task_behaviour/S2_Data/S2_data.xlsx','Sheet','FCR_states');
FCR = table2array(FCR);
FCR_diff = FCR(:,3) - FCR(:,2); 


CP_indi = mean(CP_all,5);

for i = 1:100
    for j = 1:100
        x = squeeze(CP_indi(2,:,i,j)) - squeeze(CP_indi(3,:,i,j));
        y = task(:,4); %task(:,1); FCR_diff
        [r,p] = corrcoef(x,y);
        cormat(i,j) = r(2);
        pmat(i,j) = p(2);
    end
end 
hotred = customcolormap([0 0.25 0.5 0.75 1],{'#660320','#DD7661','#F8F3F0','#4692C0','#093262'});
cormat(isnan(cormat))=0; 
figure; imagesc(cormat); colormap(hotred); caxis([-1 1.05]);

mask = pmat; mask(mask>0.05)=1; mask(mask<0.05)=0;
mask = imgaussfilt(mask,1); mask(isnan(mask))=1;
mask(mask<0.5) = 0; mask(mask>0.5)=1;
mask = bwperim(mask); 
mask(1, :) = 0; mask(end, :) = 0; mask(:, 1) = 0; mask(:, end) = 0;   % All rows of the last column
hold on; overlay = imagesc(zeros(100, 100, 3));
overlay.AlphaData = 0.8 * mask;  % 0.5 = semi-transparent black
overlay.AlphaDataMapping = 'none';

%images(cormat_mask); colormap(hotred); caxis([-1 1.05]);

cormat(pmat>0.05)=0; 
figure; imagesc(cormat); colormap(hotred); caxis([-1 1.05]);

%% 4A. Compare thalamocortical FC
var_th = zeros(2,20,14); var_cx = zeros(2,20,400);
states = {'Control', 'PreNap', 'PostNap'};

th_cx = zeros(3,20,400);
clear upper lover
for t = 1:3
    for s = 1:20
        sprintf('Sub %s, ses %s', char(string(s)), char(string(t)))
        load(sprintf('data/Parcels_regr_cat/%s/%s/%s_task-All.mat', ...
            char(SubjectName{s}), char(states{t}) ,char(SubjectName{s})))

        load(sprintf('data/Subjectdata_nii_regr_cat/%s/%s/%s_task-All.mat', ...
            upper(char(SubjectName{s})), char(states{t}) ,char(SubjectName{s})))
        thalamus = mean([subjectdata.thalamus_lh;subjectdata.thalamus_rh]);
        
        
        th_cx(t,s,:) = corr(thalamus', Parcels');

        
    end
end

wr_tc = squeeze(th_cx(1,:,:));
sd_tc = squeeze(th_cx(2,:,:));
pn_tc = squeeze(th_cx(3,:,:));

tc_change = [wr_tc(:) sd_tc(:) pn_tc(:)];

tcdat = array2table([wr_tc(:) ; sd_tc(:) ; pn_tc(:)]);
tcdat.Properties.VariableNames = {'tc_fc'};
tcdat.condition = repelem({'WR'; 'SD'; 'PN'}, 8000, 1);

writetable(tcdat, 'Fig4A.csv')

[H,P,CI,STATS] = ttest(tc_change(:,2), tc_change(:,1))
[H,P,CI,STATS] = ttest(tc_change(:,3), tc_change(:,2))
[H,P,CI,STATS] = ttest(tc_change(:,3), tc_change(:,1))


% THALAMO-CORTICAL FC
for parc = 1:400
    try
        [H,P,CI,STATS] = ttest(squeeze(th_cx(2,:,parc)),squeeze(th_cx(1,:,parc)));
        th_cx_t(parc) = STATS.tstat;
        th_cx_p(parc) = P;
    catch
        th_cx_t(parc) = NaN;
        th_cx_p(parc) = NaN;
    end
end

parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
I = customcolormap([0 0.2 0.35 0.5 0.8 0.9 1], ...
    {'#6b1200','#E23603','#FF8D33',  '#bababa', '#336EFF', '#336EFF', '#033c9e'});
clim = [-5 5];
th_cx_t(th_cx_p>0.05) = 0;
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = th_cx_t';
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.2 0.6 0.2 0.2], [0.4 0.6 0.2 0.2], ...
        1, 2, clim, I);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.2 0.4 0.2 0.2], [0.4 0.4 0.2 0.2], ...
        1, 2, clim, I);

%% 4B. correlation bt and tc change
for p = 1:400
        x = squeeze(BT_std(:,2,p)) - squeeze(BT_std(:,1,p));
        y = squeeze(th_cx(2,:,p)) - squeeze(th_cx(1,:,p));
        [bt_thal_mat(p), bt_thal_pmat(p)] = corr(x,y');
end


% Plot 
vec = bt_thal_mat;
%vec(bt_thal_pmat>0.05) = 0;
cmap = hotred; clim = [-0.6 0.6];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.1*t 0.8 0.2 0.2], [0.1*t 0.6 0.2 0.2], ...
        1, 2, clim, cmap);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.1*t 0.4 0.2 0.2], [0.1*t 0.2 0.2 0.2], ...
        1, 2, clim, cmap);

figure;
offset = 1.6;
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
colors = {One,Two,Three,Four,Five,Six,Seven};
ax1 = subplot(1,1,1); hold on;
set(ax1, 'Position', [0.1 0.1 0.35 0.6]);
plot_pvals = ones(7,2);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    
    vals = bt_thal_mat(:,Yeo_Shf==n);
    [x, f] = ksdensity(vals(:));

    [~,p,~,~] = ttest(vals(:));
    plot_pvals(n,1) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-1 1]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%% 4C. GLM with thalamus activity
glms = zeros(2,20,9,14);
clear betas_cortexl glms_cortex
for t = 1:2
    for s = 1:20
        sprintf('Sub %s, ses %s', char(string(s)), char(string(t)))
        file = sprintf('timeseries_456/concat/%s_ses-0%s_allses_timeseries.csv', ...
            char(SubjectName{s}), char(string(t)));
        
        Parcels = readtable(file, "FileType","text");
        
        Parcels = table2array(Parcels);
        thalamus = Parcels(1:613,430:443);
        
        cortex = Parcels(1:613,1:400);

        thal(t,s,:,:) = thalamus;
        hypothal(t,s,:,:) = [Parcels(1:613,412), Parcels(1:613,426)];
        cor(t,s,:,:) = cortex;
        
        idx = squeeze(idx_all(t,s,:));
        
        int_state = int_idx(s,1);
        
        if int_state > 1
            idx = idx*-1 + 3;
        end
         idx = idx*-1; % change this for segregation switches
        
        switches = [0; diff(idx)];
        
        % Perform GLM across a series of lags
        lags = -4:4;
        betas = zeros(length(lags), 2, size(switches, 2));   % Store betas for each column
        for l = 1:length(lags)
            shift = lags(l);
            reg_shifted = circshift(switches, shift);
        
            X = [ones(size(reg_shifted)), reg_shifted];  % Add intercept
        
            for i = 1:size(thalamus, 2)
                y = thalamus(:, i);
                betas_thal(l, :, i) = regress(y, X);
            end

            for i = 1:size(cortex, 2)
                y = cortex(:, i);
                betas_cortexl(l, :, i) = regress(y, X);
            end
        
        end
        glms_thal(t,s,:,:) = squeeze(betas_thal(:,2,:));
        glms_cortex(t,s,:,:) = squeeze(betas_cortexl(:,2,:));
    end 
end

% Thalamus
wr_glm = squeeze(glms_thal(1,:,:,:));

for roi = 1:14
    for t =1:9
        [H,P,CI,STATS] = ttest(wr_glm(:,t,roi));
        wr_t(t, roi) = STATS.tstat;
        wr_p(t, roi) = P;
    end
end


% Plot on thalamus image
atlas = niftiread('Thalamus_Nuclei-HCP-MaxProb.nii.gz');
info = niftiinfo('Thalamus_Nuclei-HCP-MaxProb.nii.gz');
atlas_labels = 1:14;  % Replace these with the actual integer label values in your atlas-  Example: label IDs matching region order
mask = atlas > 0;


figure;
for t = 1:9
    subplot(2, 9, t);
    set(gca, 'Visible', 'off');

    region_values = wr_t(t,:)';  % Replace with your real values
    overlay = zeros(size(atlas));
    for i = 1:length(atlas_labels)
        overlay(atlas == atlas_labels(i)) = region_values(i);
    end


    custom_map = zeros(size(atlas));
    for i = 1:length(atlas_labels)
        custom_map(atlas == atlas_labels(i)) = region_values(i);
    end
    

    fv = isosurface(mask, 0.5);   % threshold > 0 to isolate regions
    
    % Create a mask for the surface to extract real values
    vtx_vals = interp3(custom_map, fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3));
    
    
    p = patch(fv);
    set(p, 'FaceVertexCData', vtx_vals, ...
           'FaceColor', 'interp', ...
           'EdgeColor', 'none');
    
    isonormals(overlay, p);
    view(3); axis vis3d tight;
    %camlight; lighting gouraud;
    colormap(jet);   % Use diverging colormap if you have it
    caxis([-1 1]);  % Centered color scale

end


glms = zeros(2,20,9,14);
for t = 1:2
    for s = 1:20
        sprintf('Sub %s, ses %s', char(string(s)), char(string(t)))
        file = sprintf('timeseries_456/concat/%s_ses-0%s_allses_timeseries.csv', ...
            char(SubjectName{s}), char(string(t)));
        
        Parcels = readtable(file, "FileType","text");
        
        Parcels = table2array(Parcels);
        thalamus = Parcels(1:613,430:443);
        
        cortex = Parcels(1:613,1:400);
        
        idx = squeeze(idx_all(t,s,:));
        
        int_state = int_idx(s,1);
        
        if int_state > 1
            idx = idx*-1 + 3;
        end
         idx = idx; % change this for segregation switches
        
        switches = [0; diff(idx)];
        
        % Perform GLM across a series of lags
        lags = -4:4;
        betas = zeros(length(lags), 2, size(switches, 2));   % Store betas for each column
        for l = 1:length(lags)
            shift = lags(l);
            reg_shifted = circshift(switches, shift);
        
            X = [ones(size(reg_shifted)), reg_shifted];  % Add intercept
        
            for i = 1:size(thalamus, 2)
                y = thalamus(:, i);
                betas_thal(l, :, i) = regress(y, X);
            end

            for i = 1:size(cortex, 2)
                y = cortex(:, i);
                betas_cortexl(l, :, i) = regress(y, X);
            end
        
        end
        glms_thal(t,s,:,:) = squeeze(betas_thal(:,2,:));
        glms_cortex(t,s,:,:) = squeeze(betas_cortexl(:,2,:));
    end 
end

% Thalamus
wr_glm = squeeze(glms_thal(1,:,:,:));

for roi = 1:14
    for t =1:9
        [H,P,CI,STATS] = ttest(wr_glm(:,t,roi));
        wr_t(t, roi) = STATS.tstat;
        wr_p(t, roi) = P;
    end
end

for t = 1:9
    subplot(2, 9, t+9);
    set(gca, 'Visible', 'off');

    region_values = wr_t(t,:)';  % Replace with your real values
    overlay = zeros(size(atlas));
    for i = 1:length(atlas_labels)
        overlay(atlas == atlas_labels(i)) = region_values(i);
    end


    custom_map = zeros(size(atlas));
    for i = 1:length(atlas_labels)
        custom_map(atlas == atlas_labels(i)) = region_values(i);
    end
    

    fv = isosurface(mask, 0.5);   % threshold > 0 to isolate regions
    
    % Create a mask for the surface to extract real values
    vtx_vals = interp3(custom_map, fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3));
    
    
    p = patch(fv);
    set(p, 'FaceVertexCData', vtx_vals, ...
           'FaceColor', 'interp', ...
           'EdgeColor', 'none');
    
    isonormals(overlay, p);
    view(3); axis vis3d tight;
    %camlight; lighting gouraud;
    colormap(jet);   % Use diverging colormap if you have it
    caxis([-1.7 1.7]);  % Centered color scale

end

%% 4D. Shifted glms
clear row_max row_idx
wr_glm = squeeze(glms_thal(1,:,:,:));
[row_max(:,1), row_idx(:,1)] = max(mean(wr_glm,3), [], 2);

sd_glm = squeeze(glms_thal(2,:,:,:));
[row_max(:,2), row_idx(:,2)] = max(mean(sd_glm,3), [], 2);

[H,P,CI,STATS] = ttest(row_max(:,2), row_max(:,1));


for sub = 1:20
    for t =1:9
        [H,P,CI,STATS] = ttest(wr_glm(sub,t,:));
        wr_t(t, sub) = STATS.tstat;
        wr_p(t, sub) = P;
    end
end

for sub = 1:20
    for t =1:9
        [H,P,CI,STATS] = ttest(sd_glm(sub,t,:));
        sd_t(t, sub) = STATS.tstat;
        sd_p(t, sub) = P;
    end
end


% 2nd level -per roi
for roi = 1:14
    for t =1:9
        [H,P,CI,STATS] = ttest(wr_glm(:,t,roi));
        wr_t(t, roi) = STATS.tstat;
        wr_p(t, roi) = P;
    end
end

for roi = 1:14
    for t =1:9
        [H,P,CI,STATS] = ttest(sd_glm(:,t,roi));
        sd_t(t, roi) = STATS.tstat;
        sd_p(t, roi) = P;
    end
end

% Smooth
kernel = ones(1,5) / 5;
for r = 1:14
    row = wr_t(:,r);
    ks_wr(r,:) = conv(row, kernel, 'same');
end

for r = 1:14
    row = sd_t(:,r);
    ks_sd(r,:) = conv(row, kernel, 'same');
end


% Find max
clear row_max row_idx
[row_max(:,1), row_idx(:,1)] = max(ks_wr, [], 2);
[row_max(:,2), row_idx(:,2)] = max(ks_sd, [], 2);

% Test for differences
[H,P,CI,STATS] = ttest(row_idx(:,2), row_idx(:,1));


tcdat = array2table([row_idx(:,1) ; row_idx(:,2)]);
tcdat.Properties.VariableNames = {'thal_rowmax'};
tcdat.condition = repelem({'WR'; 'SD'}, 14, 1);

writetable(tcdat, 'Fig4D.csv')


%% 5A. Change in regional fluctutations
NumberOfSubjects = 20;

for t = 1:3
    sprintf('Timepoint: %s',(char(sessions(t))))
    for z = 1:NumberOfSubjects
        sprintf('Subject: %s',(char(SubjectName{z})))
        load(sprintf('data/Parcels_regr_cat/%s/%s/%s_%s.mat', ...
                     (char(SubjectName{z})),(char(sessions(t))),(char(SubjectName{z})),keyword));
        for p=1:NumberOfParcels
                STDDEV_surf_Parcels(t,z,p,:) = std(Parcels(p,:));  
        end
    end
end

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
hotgrey = customcolormap([0 0.1 0.2 0.5 0.8 1], ...
    {'#FFFEE6','#FFF000', '#FEBC00', '#807F7F', '#0004FE','#00EBFE' });
cmap_grad = cbrewer('div', 'Spectral', 256, 'linear');
clim = [-3 3];
SP = SurfStatAvSurf({'labels/fsaverage5/surf/lh.pial', ...
                     'labels/fsaverage5/surf/rh.pial'});
vec = SDEP_Tmat;
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 8);
figure;
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.3 0.8 0.2 0.2], [0.3 0.6 0.2 0.2], ...
        1, 2, clim, hotgrey);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.3 0.4 0.2 0.2], [0.3 0.2 0.2 0.2], ...
        1, 2, clim, hotgrey);

figure;
colormap(hotgrey); % e.g., parula, hot, or your own Nx3 matrix
caxis(clim)
colorbar;
view(2);

%% 5B. Scatterplot - regional flucs vs BT
% Set colormap
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
Yeo7 = customcolormap([0 0.166 0.333 0.5 0.666 0.832 1], ...
    {Seven,Six,Five,Four,Three,Two,One});

% Plot
figure;
scatter(SDEP_Tmat, tvals_bt(:,1), 100, Yeo_Shf, 'filled', ...
        'MarkerEdgeColor', 'k','LineWidth', 0.1); 
colormap(Yeo7)
h = lsline;
set(h, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');

[r ,p] = corr(tvals_bt(:,1),SDEP_Tmat)

%% 5C. Changes in global signal

load('regressors/GS.mat')

gsdat = array2table([GS.gsf(:,1); GS.gsf(:,2); GS.gsf(:,3)]);
gsdat.Properties.VariableNames = {'glob_sig_fluc'};
gsdat.condition = repelem({'WR'; 'SD'; 'PN'}, 20, 1);

writetable(gsdat, 'Fig5C.csv')

% violin plots
figure;
col = [0 0.4470 0.7410; 0.8078 0.1098 0.1882; 0.9290 0.6940 0.1250];
violin(GS.gsf(:,1:3),'mc',[],'medc','k','facecolor',col,'facealpha',1)
ylim([0 6])


%% 5D. Relationship between Global signal fluctuation and BT change
gs_sd_change = GS.gsf(:,2) - GS.gsf(:,1);

Parcels_std = nanstd(Parcels,0,4);

% Change in BT_ave and change in global signal
for p = 1:400
    [rho_ave(p), pval(p)] = corr(squeeze(BT_ave(:,2,p)) - squeeze(BT_ave(:,1,p)), ...
                             gs_sd_change);
end

% Plot to surface

vec = rho_ave;
vec(pval>0.05) = 0;
    
% Initialize output with NaNs
vtx_vals = nan(1, size(SP.coord, 2));

% Define valid indices: parc > 0 and within vector range
valid_idx = parc > 0 & parc <= numel(vec);

% Safely assign values
vtx_vals(valid_idx) = vec(parc(valid_idx));

% Plot 
figure;
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.3 0.8 0.2 0.2], [0.3 0.6 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.3 0.4 0.2 0.2], [0.3 0.2 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);

% Change in BT_sd and change in global signal
for p = 1:400
    [rho_std(p), pval(p)] = corr(squeeze(BT_std(:,2,p)) - squeeze(BT_std(:,1,p)), ...
                             gs_sd_change);
end

% Plot to surface
SP = SurfStatAvSurf({'labels/fsaverage5/surf/lh.pial', ...
                     'labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nnetw = '7'; nparc ='400'; NumberOfParcels = str2num(nparc);
load(sprintf('labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)


vec = rho_std;
vec(pval>0.05) = 0;
    
% Initialize output with NaNs
vtx_vals = nan(1, size(SP.coord, 2));

% Define valid indices: parc > 0 and within vector range
valid_idx = parc > 0 & parc <= numel(vec);

% Safely assign values
vtx_vals(valid_idx) = vec(parc(valid_idx));

% Plot 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.6 0.8 0.2 0.2], [0.6 0.6 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.6 0.4 0.2 0.2], [0.6 0.2 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);

colormap(hotred); % e.g., parula, hot, or your own Nx3 matrix
caxis([-0.7 0.7])
colorbar;
view(2);



% network dispersion
figure;

% BT_ave
offset = 1.6;
One = '#801982'; Two = '#5DACDB'; Three = '#459738'; 
Four = '#FF00BF'; Five = '#F5FCCE'; Six = '#F49E00'; Seven = '#FF0000'; 
colors = {One,Two,Three,Four,Five,Six,Seven};
ax1 = subplot(1,2,1); hold on;
set(ax1, 'Position', [0.1 0.1 0.35 0.6]);
plot_pvals = ones(7,2);
for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    
    vals = rho_ave(Yeo_Shf==n);
    [x, f] = ksdensity(vals(:));

    [~,p,~,~] = ttest(vals(:));
    plot_pvals(n,1) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-1 1]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

% BT_std
ax2 = subplot(1,2,2); hold on;
set(ax2, 'Position', [0.5 0.1 0.35 0.6]);

for n = length(colors):-1:1

    hex = colors{n};
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    
    vals = rho_std(Yeo_Shf==n);
    [x, f] = ksdensity(vals(:));

    [~,p,~,~] = ttest(vals(:));
    plot_pvals(n,2) = p(1);

    y_offset = x + (n-1)*offset;  % vertical stacking
    
    %fill(x, y_offset, x, rgb, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    fill(f, y_offset, rgb, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'k');

end
xline(0, '--k', 'LineWidth', 1.5);
xlim([-1 1]);
yticks([]);       % Remove y-axis ticks
yticklabels([]);  % Remove y-axis labels

%% 5E Physio vs GSF
%load('ecg/output/PSD/all_params.mat')

rho_mat = NaN(4,2);
p_mat = NaN(4,2);

%mean_hr
mean_hr = squeeze(mean_hr);

[rho_mat(1,1), p_mat(1,1)] = corr(GS.gsf(:,1), mean_hr(:,1), 'rows', 'complete');
[rho_mat(1,2), p_mat(1,2)] = corr(GS.gsf(:,2), mean_hr(:,2), 'rows', 'complete');
%[rho_mat(1,3), p_mat(1,3)] = corr(gs_sd_change, mean_hr(:,2) - mean_hr(:,1), 'rows', 'complete')


%std_hr
std_hr = squeeze(std_hr);
[rho_mat(2,1), p_mat(2,1)] = corr(GS.gsf(:,1), std_hr(:,1), 'rows', 'complete');
[rho_mat(2,2), p_mat(2,2)] = corr(GS.gsf(:,2), std_hr(:,2), 'rows', 'complete');
%[rho_mat(2,3), p_mat(2,3)] = corr(gs_sd_change, std_hr(:,2) - std_hr(:,1), 'rows', 'complete')


%rmssd
rmssd = squeeze(rmssd);
[rho_mat(3,1), p_mat(3,1)] = corr(GS.gsf(:,1), rmssd(:,1), 'rows', 'complete')
[rho_mat(3,2), p_mat(3,2)] = corr(GS.gsf(:,2), rmssd(:,2), 'rows', 'complete')
%[rho_mat(3,3), p_mat(3,3)] = corr(gs_sd_change, rmssd(:,2) - rmssd(:,1), 'rows', 'complete')


%delta
load('eeg/psd/PSD.MAT', 'psd_all') %sub34 missing
psd = squeeze(mean(psd_all,1));
for t  = 1:3
    for s = 1:19
        intgl = trapz(squeeze(psd_all(s,t,1:4))); 
        delta(s,t) = intgl;
    end
end
delta(20,:) = NaN;
[rho_mat(4,1), p_mat(4,1)] = corr(GS.gsf(:,1), delta(:,1), 'rows', 'complete')
[rho_mat(4,2), p_mat(4,2)] = corr(GS.gsf(:,2), delta(:,2), 'rows', 'complete')



imagesc(rho_mat); colormap(hotred); caxis([-0.8, 0.8])
% Adjust axes to line up with cells
ax = gca;
ax.XTick = 0.5:1:size(rho_mat,2)+0.5;
ax.YTick = 0.5:1:size(rho_mat,1)+0.5;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.GridColor = 'k';      % black grid lines
ax.GridAlpha = 1;        % fully opaque
ax.LineWidth = 2;
grid on;



%% 5F - ECG

clear cormat, pmat
CP_indi = mean(CP_all,5);
for i = 1:100
    for j = 1:100
        x = squeeze(CP_indi(2,:,i,j)) - squeeze(CP_indi(1,:,i,j));
        y = mean_hr(:,2) - mean_hr(:,1); %task(:,1); FCR_diff
        [r,p] = corr(x', y, 'rows','complete');
        cormat(i,j) = r;
        pmat(i,j) = p;
    end
end 
hotred = customcolormap([0 0.25 0.5 0.75 1],{'#660320','#DD7661','#F8F3F0','#4692C0','#093262'});
cormat(isnan(cormat))=0; 
figure; imagesc(cormat); colormap(hotred); caxis([-1 1.05]);

mask = pmat; mask(mask>0.05)=1; mask(mask<0.05)=0;
mask = imgaussfilt(mask,1); mask(isnan(mask))=1;
mask(mask<0.5) = 0; mask(mask>0.5)=1;
mask = bwperim(mask); 
mask(1, :) = 0; mask(end, :) = 0; mask(:, 1) = 0; mask(:, end) = 0;   % All rows of the last column
hold on; overlay = imagesc(zeros(100, 100, 3));
overlay.AlphaData = 0.8 * mask;  % 0.5 = semi-transparent black
overlay.AlphaDataMapping = 'none';

%images(cormat_mask); colormap(hotred); caxis([-1 1.05]);

cormat(pmat>0.05)=0; 
figure; imagesc(cormat); colormap(hotred); caxis([-1 1.05]);

%%
[~,p,~,stats] = ttest(mean_hr(:,2), mean_hr(:,1))
boxplot(mean_hr)


[~,p,~,stats] = ttest(std_hr(:,2), std_hr(:,1))
boxplot(std_hr)


violin(rmssd)
[~,p,~,stats] = ttest(rmssd(:,2), rmssd(:,1))


sdnn = squeeze(sdnn);
boxplot(sdnn)
[~,p,~,stats] = ttest(sdnn(:,2), sdnn(:,1))

sdnn_ch = sdnn(:,2) - sdnn(:,1);
rmssd_ch = rmssd(:,2) - rmssd(:,1);
hrstd_ch = std_hr(:,2) - std_hr(:,1);

hr_ch = mean_hr(:,2) - mean_hr(:,1);

BT_ave_mean = table2array(BT_ave_plot);
ch_bt = BT_ave_mean(:,2) - BT_ave_mean(:,1);

scatter(ch_bt, hr_ch, 'filled')
[rho, pval] = corr(ch_bt, hr_ch, 'rows', 'complete')

% Change in BT_sd and change in HR

for p = 1:400
    [rho_hr(p), pval(p)] = corr(squeeze(BT_std(:,2,p)) - squeeze(BT_std(:,1,p)), ...
                             hr_ch, 'rows','complete');
end

% Plot to surface
SP = SurfStatAvSurf({'labels/fsaverage5/surf/lh.pial', ...
                     'labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf_Yeo(:,1)' rh_labels_Shf_Yeo(:,1)'];
nnetw = '7'; nparc ='400'; NumberOfParcels = str2num(nparc);
load(sprintf('labels/fsaverage5/ShfParcels/ShfLabels%s_%s.mat',nparc,nnetw)) 
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
end
vertices = SP.coord';         % Now size = (20484 x 3)
faces = SP.tri;               % Already (40960 x 3)


vec = rho_hr;
%vec(pval>0.05) = 0;
    
% Initialize output with NaNs
vtx_vals = nan(1, size(SP.coord, 2));

% Define valid indices: parc > 0 and within vector range
valid_idx = parc > 0 & parc <= numel(vec);

% Safely assign values
vtx_vals(valid_idx) = vec(parc(valid_idx));

% Plot 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.3 0.8 0.2 0.2], [0.3 0.6 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.3 0.4 0.2 0.2], [0.3 0.2 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);


vec(pval>0.05) = 0;
    
% Initialize output with NaNs
vtx_vals = nan(1, size(SP.coord, 2));

% Define valid indices: parc > 0 and within vector range
valid_idx = parc > 0 & parc <= numel(vec);

% Safely assign values
vtx_vals(valid_idx) = vec(parc(valid_idx));

% Plot 
nonwall = [2:201 203:402];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = vec;
OnSurf  = BoSurfStatMakeParcelData(toMap, SP, parc);
OnSurfS = SurfStatSmooth(OnSurf, SP, 4);
BoSurfStat_calibrate2Views(OnSurfS, SP, ...
        [0.6 0.8 0.2 0.2], [0.6 0.6 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);
BoSurfStat_calibrate2Views_rh(OnSurfS, SP, ...
        [0.6 0.4 0.2 0.2], [0.6 0.2 0.2 0.2], ...
        1, 2, [-0.7 0.7], hotred);


%% 5G. EEG topo

eeg_t = zeros(185,30);
for e = 1:185
    sprintf('Electrode %s', char(string(e)))
        for f = 1:30
            clear s1 s2
            for s = 1:17
                load( sprintf('eeg/psd/%s_psd.mat',SubjectName{s}), 'psdo_all')
  
                s1(s) = squeeze(psdo_all(1,e,f));
                s2(s) = squeeze(psdo_all(2,e,f));

            end
            [~,p,~,stats] = ttest(s2, s1);
            eeg_t(e,f) = stats.tstat;
            eeg_p(e,f) = p;
        end
end


opts = struct('session',1, 'freqs', 1:30, 'band_hz',[1 4], ...
    'logscale', false, 'clim',[5 15]);


plot_topomap_spectra_rel(eeg_t, 'eeg/psd/coordinates.xml', 'eeg/psd/channels.csv', ...
    struct('band_hz',[1 4], 'clim',[0.9 1.01], 'background','transparent'));




%% 5H Delta change vs BT change
load('eeg/psd/PSD.MAT', 'psd_all') %sub34 missing

psd = squeeze(mean(psd_all,1));


for t  = 1:3
    for s = 1:19
        intgl = trapz(squeeze(psd_all(s,t,1:4))); 

        total = trapz(squeeze(psd_all(s,t,:)));
        delta(s,t) = (intgl/total)*100;

        %delta(s,t) = intgl


    end
end

BT_std_mean = table2array(BT_std_plot);
BT_std_mean(:,1)
ch_bt = BT_std_mean(:,2) - BT_std_mean(:,1);

ch_delta = delta(:,2) - delta(:,1);


scatter(ch_bt, ch_delta, 'filled'); lsline

[r,p] = corr(ch_delta, ch_bt, 'rows', 'complete')


cordat = array2table([ch_bt, ch_delta]);
cordat.Properties.VariableNames = {'BTstdWRSD', 'DeltastdWRSD'};

writetable(cordat, 'Fig5H.csv','WriteVariableNames', 1)


