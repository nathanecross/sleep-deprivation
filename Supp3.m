%
% Supplementary Figure 3
%
% Changes in gradients following GSR

%% Load dependencies and setup

% Add dependencies
addpath('dependencies/customcolormap');
addpath('dependencies/fireice');
addpath('dependencies/cbrewer');
addpath('dependencies/DNorm2');
addpath('dependencies/surfstat');
addpath('gradients');

% Load data
load('../code_output_gsr/gradients/embedding.mat');
load('../code_output_gsr/gradients/indi_embed.mat');
load('labels/Yeo17_Shf400.mat');
load('SubjectName');
load('labels/fsaverage5/ShfParcels/ShfLabels400_17.mat') 

times = {'Control','PreNap','PostNap'};
% Import Freesurfer surface for plotting
FS = SurfStatAvSurf({'../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/lh.pial', '../../Parcellations/FreeSurfer5.3/subjects/fsaverage5/surf/rh.pial'});
parc = [lh_labels_Shf(:,1)' rh_labels_Shf(:,1)'];
nonwall = [2:201 203:402];

%% Figure 4A 
% Range of gradients across states
for g = 1:3
    for z = 1:length(SubjectName)
        for t = 1:length(times)
            maxx(g,z,t) = max(indi_embed(:,g,z,t));
            minn(g,z,t) = min(indi_embed(:,g,z,t));
            siz(g,z,t) = max(indi_embed(:,g,z,t)) - min(indi_embed(:,g,z,t));
            bin_val(t,g,:,z) = slidingWindowEqualDistrib((indi_embed(:,g,z,t))', (indi_embed(:,g,z,t))', 50, 10);
            siz_bin(g,z,t) = bin_val(t,g,50,z) - bin_val(t,g,1,z);
        end
    end
end

for g=1:3
    figure;
    b = barh([mean(maxx(g,:,1)),mean(maxx(g,:,2)),mean(maxx(g,:,3))]);
    xlim([-0.2 0.18]);
    hold on; 
    er = errorbar([mean(maxx(g,:,1)),mean(maxx(g,:,2)),mean(maxx(g,:,3))],b.XData,[std(maxx(g,:,1)), ...
        std(maxx(g,:,2)),std(maxx(g,:,3))],[std(maxx(g,:,1)), ...
        std(maxx(g,:,2)),std(maxx(g,:,3))],'horizontal');    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    b.FaceColor = 'flat';
    hold on; barh([mean(minn(g,:,1)),mean(minn(g,:,2)),mean(minn(g,:,3))]);
    er = errorbar([mean(minn(g,:,1)),mean(minn(g,:,2)),mean(minn(g,:,3))],b.XData,[std(minn(g,:,1)), ...
        std(minn(g,:,2)),std(minn(g,:,3))],[std(minn(g,:,1)), ...
        std(minn(g,:,2)),std(minn(g,:,3))],'horizontal');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
end

%% Figure 4B
% Distribution plots of Gradients
for g=1:3
    figure;
    ksdensity(embedding(:,g,1)); hold on; ksdensity(embedding(:,g,2)); hold on; ksdensity(embedding(:,g,3));
    ylim([0 18])
    xlim([-0.2 0.2])
end

%% Figure 4C,D,E
% Colorsphere in Figure 4E was created in Adobe Illustrator
subs = SubjectName;
%embedding = embedding*-1;

for prep_figure4C = 1
    % Calculate distances to all other nodes in gradient space    
    for s = 1:3
        for ii = 1:length(subs)
            for p = 1:400
                grad_dist(p,:,ii,s) = DNorm2(indi_embed(:,1:3,ii,s) - indi_embed(p,1:3,ii,s), 2);
            end
        end    
    end
end

for prep_figure4D = 1
    clear stats
    for p = 1:400 % for each node
        in_table = table([1:size(indi_embed,3)]', squeeze(mean(grad_dist(p,:,:,1))), ...
            squeeze(mean(grad_dist(p,:,:,2))), squeeze(mean(grad_dist(p,:,:,3))), ...
            'VariableNames', {'subject', 'state1', 'state2', 'state3'});
        rm = fitrm(in_table, 'state1-state3~1');
        ranovatbl = ranova(rm);
        stats(p,:) = table2array(ranovatbl(1,:));
    end
    rois = fdr_bh(stats(:,5)); 


    subs = SubjectName;
    % calculate nodal shifts between states
    % Euclidean distance in 3D space
    for ii = 1:length(subs)
        shift(:,ii,1) = DNorm2(indi_embed(:,1:3,ii,2) - indi_embed(:,1:3,ii,1), 2);
        shift(:,ii,2) = DNorm2(indi_embed(:,1:3,ii,3) - indi_embed(:,1:3,ii,2), 2);
    end

    % test whether the shift is non-zero
%     for jj = 1:2 % for each state shift
%         [~,shift_p(:,jj),~,stats] = ttest(shift(:,:,jj)');
%         shift_t(:,jj) = stats.tstat;
%     end
%     [~,~,~,shift_p] = fdr_bh(shift_p); % fdr correction

    % cross-correlate spatial t-stat maps
%     corr(shift_t); 
end
   
f = figure('units', 'centimeters', 'position', [0 0 20 20]);
cmap_cent = hot;
cmap_shift = cbrewer('seq', 'Purples', 256, 'linear');
scatter_axes = [-0.12 0.16 -0.09 0.1 -0.1 0.11];

for figure_4C = 1
    a(1) = axes('position', [0.05 0.5 0.3 0.3]);
    scatter3(embedding(:,1),embedding(:,2),embedding(:,3),50, ...
        squeeze(mean(mean(grad_dist(:,:,:,1),3))),'filled', 'MarkerEdgeColor', 'k')
    view(-6,5)
    colormap(a(1), cmap_cent)
    axis(scatter_axes)

    toMap = ones(1,length(unique(parc)));
    toMap(nonwall) = squeeze(mean(mean(grad_dist(:,:,:,1),3)));
    OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 4);
    BoSurfStat_calibrate2Views(OnSurfS, FS, ...
        [0.05 0.82 0.15 0.15], [0.2 0.82 0.15 0.15], ...
        1, 2, [prctile(toMap,2) prctile(toMap,98)], cmap_cent);

    a(3) = axes('position', [0.15 0.83 0.1 0.01]);
    imagesc(1:100); axis off
    colormap(a(3), flip(cmap_cent))
end
for figure_4D = 1
    a(2) = axes('position', [0.5 0.5 0.3 0.3]);
    scatter3(embedding(:,1),embedding(:,2),embedding(:,3),50, ...
        [0.5 0.5 0.5], 'filled')
    alpha 0.2
    view(-6,5)
    hold on
    scatter3(embedding(rois,1),embedding(rois,2),embedding(rois,3),50, ...
        stats(rois,4) .* (fdr_bh(stats(rois,5))),'filled', 'MarkerEdgeColor', 'k')
    colormap(a(2), cmap_shift)
    axis(scatter_axes)

    toMap = ones(1,length(unique(parc)));
    toMap(nonwall) = stats(:,4) .* (fdr_bh(stats(:,5)));
    OnSurf  = BoSurfStatMakeParcelData(toMap, FS, parc);
    OnSurfS = SurfStatSmooth(OnSurf, FS, 4);
    BoSurfStat_calibrate2Views(OnSurfS, FS, ...
        [0.5 0.82 0.15 0.15], [0.65 0.82 0.15 0.15], ...
        1, 2, [prctile(toMap,2) prctile(toMap,98)], cmap_shift);

    a(4) = axes('position', [0.6 0.83 0.1 0.01]);
    imagesc(1:100); axis off
    colormap(a(4), cmap_shift)

    %exportfig(f, [mainDir '\Figure2A.png'], 'Format'  , 'png', ...
    %        'color', 'cmyk', 'Resolution', 600, 'FontSize', 0.01)
end

for prep_figure4E = 1
    
    % find difference along each gradient
    gdiff(:,:,:,1) = indi_embed(:,:,:,2) - indi_embed(:,:,:,1);
    gdiff(:,:,:,2) = indi_embed(:,:,:,3) - indi_embed(:,:,:,2);
    gdiff(:,:,:,3) = indi_embed(:,:,:,3) - indi_embed(:,:,:,1);

    % create colour axes
    cmap_g1 = customcolormap([0 0.5 1], [1 1 0; 1 1 1; 0 0 1], 50);
    cmap_g2 = customcolormap([0 0.5 1], [1 0 0; 1 1 1; 0 1 1], 50);
    cmap_g3 = customcolormap([0 0.5 1], [1 0 1; 1 1 1; 0 1 0], 50);

       % test whether non-zero and get t-stat
    for jj = 1:3
        for g = 1:3
            [~,gdiff_p(:,jj,g),~,stats] = ttest(squeeze(gdiff(:,g,:,jj))');
            gdiff_t(:,jj,g) = stats.tstat;
            gdiff_sd(:,jj,g) = stats.sd;
        end
    end

    % convert shifts to colours
    these_rois = find(rois);
    ruler = rescale(1:length(cmap_g1),-4,4); 
    ruler_idx = ones(size(gdiff_t)) * length(ruler)/2; % set at white 
    for jj = 1:3
        for g = 1:3
            for n = 1:length(these_rois)
                % match each t-stat to a cell in the colour axes
                [~, ruler_idx(these_rois(n),jj,g)] = min(abs(diff([repmat(double(gdiff_t(these_rois(n),jj,g)), 1, length(ruler)); ruler])));
            end
        end
        gdiff_colour(:,:,jj) = ((cmap_g1(ruler_idx(:,jj,1),:) + cmap_g2(ruler_idx(:,jj,2),:) + cmap_g3(ruler_idx(:,jj,3),:))/3);
        
    end

    % impute colours to the cortical surface
    to_colour = ones(length(unique(parc)), 3);
    to_colour(nonwall,:) = gdiff_colour(:,:,1);
    for col = 1:3
        gdiff1_onsurf(:,col) = BoSurfStatMakeParcelData(to_colour(:,col), FS, parc);
        gdiff1_onsurfS(:,col) = SurfStatSmooth(gdiff1_onsurf(:,col)', FS, 3);
    end
    to_colour(nonwall,:) = gdiff_colour(:,:,2);
    for col = 1:3
        gdiff2_onsurf(:,col) = BoSurfStatMakeParcelData(to_colour(:,col), FS, parc);
        gdiff2_onsurfS(:,col) = SurfStatSmooth(gdiff2_onsurf(:,col)', FS, 3);
    end
    to_colour(nonwall,:) = gdiff_colour(:,:,3);
    for col = 1:3
        gdiff3_onsurf(:,col) = BoSurfStatMakeParcelData(to_colour(:,col), FS, parc);
        gdiff3_onsurfS(:,col) = SurfStatSmooth(gdiff3_onsurf(:,col)', FS, 3);
    end
end   

hold on;
for figure_4E = 1
    
    colourSurface(gdiff1_onsurf, FS, [0.295 0.31 0.08 0.08; 0.38 0.31 0.08 0.08; ...
        0.295 0.25 0.08 0.08; 0.38 0.25 0.08 0.08]);
    
    colourSurface(gdiff2_onsurf, FS, [0.48 0.31 0.08 0.08; 0.565 0.31 0.08 0.08; ...
        0.48 0.25 0.08 0.08; 0.565 0.25 0.08 0.08]);
    
    colourSurface(gdiff3_onsurf, FS, [0.665 0.31 0.08 0.08; 0.75 0.31 0.08 0.08; ...
        0.665 0.25 0.08 0.08; 0.75 0.25 0.08 0.08]);
    
    indi_embed = indi_embed.*-1;
    for g = 1:3
        a(g) = axes('position', [0.3 0.25-(g*0.06) 0.15 0.05])
        for r = 1:length(these_rois)
            errorbar(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,1:2))), ...
                squeeze(std(indi_embed(these_rois(r),g,:,1:2))), ...
                'Color', gdiff_colour(these_rois(r),:,1))
            hold on
            scatter(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,1:2))), 10, ...
                gdiff_colour(these_rois(r),:,1), 'filled');
            set(gca,'xtick',[]);
        end
        xlim([0.5 2.5])
        if g == 1
            ylim([-0.1 0.2])
        end
        a(g) = axes('position', [0.485 0.25-(g*0.06) 0.15 0.05])
        for r = 1:length(these_rois)
            errorbar(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,2:3))), ...
                squeeze(std(indi_embed(these_rois(r),g,:,2:3))), ...
                'Color', gdiff_colour(these_rois(r),:,2))
            hold on
            scatter(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,2:3))), 10, ...
                gdiff_colour(these_rois(r),:,2), 'filled')
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        end
        xlim([0.5 2.5])
        
        
        a(g) = axes('position', [0.67 0.25-(g*0.06) 0.15 0.05])
        for r = 1:length(these_rois)
            errorbar(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,[1 3]))), ...
                squeeze(std(indi_embed(these_rois(r),g,:,[1 3]))), ...
                'Color', gdiff_colour(these_rois(r),:,3))
            hold on
            scatter(1:2, squeeze(mean(indi_embed(these_rois(r),g,:,[1 3]))), 10, ...
                gdiff_colour(these_rois(r),:,3), 'filled')
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        end
        xlim([0.5 2.5])
    end
      
%     exportfig(f, [mainDir '\Figure2B_watercolour.png'], 'Format'  , 'png', ...
%         'color', 'cmyk', 'Resolution', 600, 'FontSize', 0.001)
    
end