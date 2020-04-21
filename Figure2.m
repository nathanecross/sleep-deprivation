% 
% Code for Figure 2
%

%% Load dependencies
addpath('dependencies/surfstat');
load('gradients/embedding.mat');
load('labels/fsaverage5/ShfParcels/ShfLabels400_17.mat') 
load('labels/Yeo17_Shf400.mat');

% Import Freesurfer surface for plotting
FS = SurfStatAvSurf({'../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/lh.pial', '../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf(:,1)' rh_labels_Shf(:,1)'];
nonwall = [2:201 203:402];

% Load Yeo 17 Network template on fsaverage5 surface
Yeo_17Clusters_ref = load('labels/fsaverage5/YeoNetworks/1000subjects_clusters017_ref.mat'); 


% Squeeze first 3 gradients
u(:,1) = 1:400;
u(:,2:4) = embedding(:,1:3,1);


%% Figure 2A
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
        [0.01 0.95-(g*0.3) 0.4 0.4], [0.45 0.95-(g*0.3) 0.4 0.4], ...
        1, 2, clim, cmap_grad); 
end      
a(10) = axes('position', [0.18 0.01 0.5 0.03]);
colormap(a(10), cmap_grad)
imagesc(1:100); axis off

% fname = char('../Gradients/Figure2/BL_grads.png');
% exportfig(f, fname, 'Format', 'png', 'Resolution', 600, 'FontSize', 0.8)

%% Figure 2B
% Figure in the paper was created in Graphpad Prism

for g = 1:3
    figure;
    grp = Yeo17_Shf400;
    G = embedding(:,g,1);
    if g ==1
        G1 = embedding(:,g,1);
    elseif g ==2
        G2 = embedding(:,g,1);
    elseif g ==3
        G3 = embedding(:,g,1);
    end
    boxplot(G, grp);    
end

%% Figure 2C
% Diagram of Figure 2C was created in Adobe Illustrator

zcmap1 = 1-(((max(u(:,2))-(u(:,2)))/(max(u(:,2))-min(u(:,2)))));
zcmap2 = 1-(((max(u(:,3))-(u(:,3)))/(max(u(:,3))-min(u(:,3)))));
zcmap3 = 1-(((max(u(:,4))-(u(:,4)))/(max(u(:,4))-min(u(:,4)))));
zcmap = [zcmap1 zcmap2 zcmap3];

x = u(:,2);
y = u(:,3);
z = u(:,4);
s = 50;
c = zcmap;
h = scatter3(x,y,z,s,c,'filled', 'MarkerEdgeColor', 'k');
set(gca,'XLim',[-0.12 0.16],'YLim',[-0.09 0.15],'ZLim',[-0.1 0.1])
xlabel('G1','FontSize',12,'FontWeight','bold')
ylabel('G2','FontSize',12,'FontWeight','bold')
zlabel('G3','FontSize',12,'FontWeight','bold')
cmap = zcmap;
colormap(gca, cmap)
view(-6,5)

view(85,5)
view(-131,5)

cmap = [0.9 0.9 0.9; zcmap];
toMap = zeros(1,length(unique(parc)));
toMap(nonwall) = u(:,1);
clim = [min(toMap) max(toMap)];
OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
BoSurfStat_calibrate2Views(OnSurf, FS, ...
        [0.01 0.5 0.4 0.4], [0.45 0.5 0.4 0.4], ...
        1, 2, clim, cmap); 

