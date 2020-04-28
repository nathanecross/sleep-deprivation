%
% Figure 3
%

%% Load dependencies and setup
addpath('../../dependencies/surfstat');
addpath('../../dependencies/shadedErrorBar.m');
load('../../labels/fsaverage5/ShfParcels/ShfLabels400_17.mat') 
load('../embedding.mat');
tasks = {'ANT', 'Nback', 'PVT'};
NumberOfParcels = 400;

% Import Freesurfer surface for plotting
FS = SurfStatAvSurf({'../../labels/fsaverage5/surf/lh.pial', '../../labels/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf(:,1)' rh_labels_Shf(:,1)'];
nonwall = [2:201 203:402];

% Squeeze first 3 gradients and flip for viewing
u(:,1) = 1:400;
u(:,2:4) = embedding(:,1:3,1)*-1;

%% 3A. Obtain active areas in response to task onset
% all tasks
clear BL PRE POST
for t = 1:length(tasks)
    load(sprintf('data/activations_%s.mat',tasks{t}));
    BL.(char(sprintf('%s',tasks{t}))) = Activations.BL; 
    PRE.(char(sprintf('%s',tasks{t}))) = Activations.PRE; 
    POST.(char(sprintf('%s',tasks{t}))) = Activations.POST; 
    [~,~,~,stats] = ttest(BL.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(1,:) = stats.tstat;
    [~,~,~,stats] =ttest(PRE.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(2,:) = stats.tstat;
    [~,~,~,stats] = ttest(POST.(char(sprintf('%s',tasks{t}))));
    actis.(char(sprintf('group%svec',tasks{t})))(3,:) = stats.tstat;
    for p = 1:NumberOfParcels
        [~,pv,~,stats] = ttest(BL.(char(sprintf('%s',tasks{t})))(:,p), PRE.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,1,:) = [0 0 stats.sd stats.tstat pv];
        [~,pv,~,stats] = ttest(PRE.(char(sprintf('%s',tasks{t})))(:,p), POST.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,2,:) = [0 0 stats.sd stats.tstat pv];
        [~,pv,~,stats] = ttest(BL.(char(sprintf('%s',tasks{t})))(:,p), POST.(char(sprintf('%s',tasks{t})))(:,p));
        actis.(char(sprintf('%s_stats',tasks{t})))(p,3,:) = [0 0 stats.sd stats.tstat pv];
    end
end

%% Figure 3A
% PVT Task activations at WR state
vec = actis.(sprintf('group%svec',string(tasks(3))))(1,:)';
[~,p] = ttest(BL.(char(sprintf('%s',string(tasks(3))))));
p = [(1:400)' p'];
toMap = zeros(1,length(unique(parc)));
        toMap(nonwall) = vec;
        OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
f = figure;
OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
clim = [-4 4];
I = customcolormap([0 0.2 0.5 0.8 1], {'#FF8D33', '#FFE633', '#bababa', '#33C7FF', '#336EFF'});
pos = [0.01 0.6 0.4 0.4; 0.01 0.25 0.4 0.4; 0.35 0.25 0.4 0.4; 0.35 0.6 0.4 0.4];
BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, I); 

figure;
x = u(:,2);
y = u(:,3);
z = u(:,4);
s = 70;
c = vec;
h = scatter3(x,y,z,s,c,'filled', 'MarkerEdgeColor', 'k');
set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
xlabel('G1','FontSize',12,'FontWeight','bold')
ylabel('G2','FontSize',12,'FontWeight','bold')
zlabel('G3','FontSize',12,'FontWeight','bold')
colormap(gca, I)
view(-6,5)
scatter3(embedding(:,1),embedding(:,2),embedding(:,3),s+10, ...
    [0.5 0.5 0.5], 'filled')
set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
alpha 0.2
view(-6,5)
hold on
 scatter3(embedding(p(:,2)<0.05,1),embedding(p(:,2)<0.05,2),embedding(p(:,2)<0.05,3),s+10, ...
     vec(p(:,2)<0.05,1),'filled', 'MarkerEdgeColor', 'k')
 set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
colormap(gca, I)
caxis([-4 4])

%% Figure 3B
% PVT Task activations long gradient bins
f = figure('units', 'centimeters', 'position', [0 0 20 20]);
cmap_grad = cbrewer('div', 'Spectral', 256, 'linear');
cmap_shift = cbrewer('seq', 'YlOrRd', 256, 'linear');
clim_shift = [6 12];
scatter_axes = [-0.12 0.16 -0.09 0.1 -0.1 0.11];

for g = 1:3
    toMap = zeros(1,length(unique(parc)));
    toMap(nonwall) = embedding(:,g,1);
    clim = [min(toMap) max(toMap)];
    OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 4);
    BoSurfStat_calibrate2Views(OnSurfS, FS, ...
        [0.12+((g-1)*0.278) 0.5 0.15 0.15], [0.21+((g-1)*0.278) 0.5 0.15 0.15], ...
        1, 2, clim, cmap_grad); 
end      
hold on;
tasktrace = actis.(sprintf('group%svec',string(tasks(3))));
stats_trace = actis.(sprintf('%s_stats',string(tasks(3))));
I = customcolormap([0 0.2 0.5 0.8 1], {'#FF8D33', '#FFE633', '#7C7C7C', '#33C7FF', '#336EFF'});
J = customcolormap([0 0.2 0.3 0.7 0.8 1], {'#AF0000','#ea3030','#7C7C7C','#7C7C7C','#8E5AFC','#3800AF'});
for g=1:3
    %WR activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(1,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,1,3))'/2)),u(:,1+g)', 50, 10);
    p1 = color_line(1:50,fity,fity,'FaceColor',[15/225, 146/225, 225/255]);
    p1.LineWidth = 2;
    p1.EdgeColor = [15/225, 146/225, 225/255];
    ylim([-5 5])
    xlim([0 50])
    pa = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[0, 90/225, 204/255]});
    colormap(gca,J)
    caxis([-2 2])
    title(sprintf('Gradient %s',string(g)))
    hold on;
    %SD activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(2,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,2,3))'/2)),u(:,1+g)', 50, 10);
    p2 = color_line(1:50,fity,fity,'FaceColor',[211/255, 42/255, 67/255]);
    p2.LineWidth = 2;
    ylim([-5 5])
    xlim([0 50])
    pb = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[211/255, 42/255, 67/255]});
    colormap(gca,J)
    caxis([-2 2])
    p2.EdgeColor = [211/255, 42/255, 67/255];
    hold on;
    %PRN activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(3,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,2,3))'/2)),u(:,1+g)', 50, 10);
    p3 = color_line(1:50,fity,fity,'FaceColor',[1, 192/255, 0]);
    p3.LineWidth = 2;
    p3.EdgeColor = [1, 192/255, 0];
    ylim([-5 5])
    yticks([-5 -2.5 0 2.5 5])
    xlim([0 50])
    pc = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[1, 192/255, 0]});
    colormap(gca,J)
    caxis([-2 2])
    plot(zeros(1,50),'--','Color','k');
    if g == 3
        legend([pa.edge(1) pb.edge(1) pc.edge(1)],{'WR','SD','PRN'})
        legend('boxoff')
    end
    hold off;
    
end

%% Figure 3C
% 2-back Task activations at WR state
vec = actis.(sprintf('group%svec',string(tasks(2))))(1,:)';
toMap = zeros(1,length(unique(parc)));
        toMap(nonwall) = vec;
        OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
f = figure;
OnSurfS = SurfStatSmooth(OnSurf, FS, 8);
clim = [-4 4];
I = customcolormap([0 0.2 0.5 0.8 1], {'#FF8D33', '#FFE633', '#bababa', '#33C7FF', '#336EFF'});
pos = [0.01 0.6 0.4 0.4; 0.01 0.25 0.4 0.4; 0.35 0.25 0.4 0.4; 0.35 0.6 0.4 0.4];
BoSurfStat_calibrate4Views(OnSurfS, FS, pos, [1;2;3;4], clim, I); 

figure;
x = u(:,2);
y = u(:,3);
z = u(:,4);
s = 70;
c = vec;
h = scatter3(x,y,z,s,c,'filled', 'MarkerEdgeColor', 'k');
set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
xlabel('G1','FontSize',12,'FontWeight','bold')
ylabel('G2','FontSize',12,'FontWeight','bold')
zlabel('G3','FontSize',12,'FontWeight','bold')
cmap_shift = cbrewer('seq', 'Purples', 256, 'linear');
colormap(gca, I)
view(-6,5)
scatter3(embedding(:,1),embedding(:,2),embedding(:,3),s+10, ...
    [0.5 0.5 0.5], 'filled')
set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
alpha 0.2
view(-6,5)
hold on
 scatter3(embedding(p(:,2)<0.05,1),embedding(p(:,2)<0.05,2),embedding(p(:,2)<0.05,3),s+10, ...
     vec(p(:,2)<0.05,1),'filled', 'MarkerEdgeColor', 'k')
 set(gca,'XLim',[-0.13 0.16],'YLim',[-0.09 0.1],'ZLim',[-0.12 0.1])
colormap(gca, I)
caxis([-4 4])

%% Figure 3D
% 2-back Task activations long gradient bins
f = figure('units', 'centimeters', 'position', [0 0 20 20]);
cmap_grad = cbrewer('div', 'Spectral', 256, 'linear');
cmap_shift = cbrewer('seq', 'YlOrRd', 256, 'linear');
clim_shift = [6 12];
scatter_axes = [-0.12 0.16 -0.09 0.1 -0.1 0.11];

for g = 1:3
    toMap = zeros(1,length(unique(parc)));
    toMap(nonwall) = embedding(:,g,1);
    clim = [min(toMap) max(toMap)];
    OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 4);
    BoSurfStat_calibrate2Views(OnSurfS, FS, ...
        [0.12+((g-1)*0.278) 0.5 0.15 0.15], [0.21+((g-1)*0.278) 0.5 0.15 0.15], ...
        1, 2, clim, cmap_grad); 
end      
hold on;
tasktrace = actis.(sprintf('group%svec',string(tasks(2))));
stats_trace = actis.(sprintf('%s_stats',string(tasks(2))));
I = customcolormap([0 0.2 0.5 0.8 1], {'#FF8D33', '#FFE633', '#7C7C7C', '#33C7FF', '#336EFF'});
J = customcolormap([0 0.2 0.3 0.7 0.8 1], {'#AF0000','#ea3030','#7C7C7C','#7C7C7C','#8E5AFC','#3800AF'});
for g=1:3
    %WR activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(1,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,1,3))'/2)),u(:,1+g)', 50, 10);
    p1 = color_line(1:50,fity,fity,'FaceColor',[15/225, 146/225, 225/255]);
    p1.LineWidth = 2;
    p1.EdgeColor = [15/225, 146/225, 225/255];
    ylim([-5 5])
    xlim([0 50])
    pa = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[0, 90/225, 204/255]});
    colormap(gca,J)
    caxis([-2 2])
    title(sprintf('Gradient %s',string(g)))
    hold on;
    %SD activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(2,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,2,3))'/2)),u(:,1+g)', 50, 10);
    p2 = color_line(1:50,fity,fity,'FaceColor',[211/255, 42/255, 67/255]);
    p2.LineWidth = 2;
    ylim([-5 5])
    xlim([0 50])
    pb = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[211/255, 42/255, 67/255]});
    colormap(gca,J)
    caxis([-2 2])
    p2.EdgeColor = [211/255, 42/255, 67/255];
    hold on;
    %PRN activations on each gradient
    subplot(2,3,g+3);
    fity = slidingWindowEqualDistrib(tasktrace(3,:), u(:,1+g)', 50, 10);
    fitsd = slidingWindowEqualDistrib(((squeeze(stats_trace(:,2,3))'/2)),u(:,1+g)', 50, 10);
    p3 = color_line(1:50,fity,fity,'FaceColor',[0.5, 192/255, 0]);
    p3.LineWidth = 2;
    p3.EdgeColor = [1, 192/255, 0];
    ylim([-5 5])
    yticks([-5 -2.5 0 2.5 5])
    xlim([0 50])
    pc = shadedErrorBar(1:50,fity,fitsd,'lineProps',{'color',[1, 192/255, 0]});
    colormap(gca,J)
    caxis([-2 2])
    plot(zeros(1,50),'--','Color','k');
    if g == 3
        legend([pa.edge(1) pb.edge(1) pc.edge(1)],{'WR','SD','PRN'})
        legend('boxoff')
    end
    hold off;
    
end
