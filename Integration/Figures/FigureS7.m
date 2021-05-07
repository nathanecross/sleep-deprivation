%% Figure S7.A Connectivity Matrix
keyword = 'task-All'; NumberOfSubjects=20;
NumberOfParcels = 400; nparc = '400'; nnetw= '17';
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}]; 
contrasts = {'Control','PreNap'}; %'#9e1d03'
cmap = customcolormap([0 0.2 0.4 0.5 0.8 0.9 1], {'#6b1200','#E23603','#FF8D33',  '#bababa', '#336EFF', '#336EFF', '#033c9e'});

load('/Volumes/NATHAN/SDEP/code_and_output/data/GLM400_task-All.mat')

matrix = GLM.Recovery_Tmat;
[matrix_ordered, matrix_bordered, L] = reorder_matrices(matrix, nparc, nnetw); 

figure;   
imagesc(matrix_bordered); 
set(gca, 'XTick', 1:length(L)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(L)); % center y-axis ticks on bins
set(gca, 'XTickLabel', L, 'FontSize', 5); % set x-axis labels
set(gca, 'YTickLabel', L, 'FontSize', 5);
caxis([-7, 7]);
colormap(cmap); % <--- go to town and play around on this 

matrix(GLM.SDEP_Pmat>0.05)=0;
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



%% Figure S7.C. Hierarchical Model
clear a
load('data/Integration400_17_task-All.mat');
load('data/Integration400_7_task-All.mat');
data = NaN(17,3); data2 = NaN(17,3);
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
data = table2array(readtable('/Volumes/NATHAN/SDEP/Integration_activations/HI_prn.xlsx','ReadVariableNames',false));
data(1:7,2) = [data(6,2);data(7,2);data(3,2);data(5,2);data(4,2);data(1,2);data(2,2)];
data(1:17,3) = [data(11,3);data(12,3);data(13,3);data(17,3);data(14,3);data(15,3);data(16,3);data(5,3);data(6,3); ...
             data(9,3);data(10,3);data(7,3);data(8,3);data(3,3);data(4,3);data(1,3);data(2,3)];

ass_statz = HI_17.Assemblies.stats_alphabet_order_prn';
ass_statz(ass_statz==0)=NaN;

data = [data [ass_statz; nan(length(data)-size(ass_statz,1),17)]];
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
clim = [-30 30];


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

% networks 1:17
for i = 9:25
    if i < 13
        a(i) = axes('position', [0.47 1.36-(0.05*i) 0.04 0.037]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    elseif i < 16
        a(i) = axes('position', [0.47 1.32-(0.05*i) 0.04 0.04]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    elseif i < 18
        a(i) = axes('position', [0.47 1.215-(0.045*i) 0.04 0.03]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    elseif i < 20
        a(i) = axes('position', [0.47 1.205-(0.045*i) 0.04 0.035]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    elseif i < 22
        a(i) = axes('position', [0.47 1.195-(0.045*i) 0.04 0.03]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)   
    elseif i < 24
        a(i) = axes('position', [0.47 1.185-(0.045*i) 0.04 0.03]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
    else
        a(i) = axes('position', [0.47 1.175-(0.045*i) 0.04 0.03]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
        
    end
end

% assemblies

for i = 26:82
    if i < 38
        if i == 26 || i == 27 || i == 28  
            a(i) = axes('position', [0.6 1.198-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        elseif  i == 29 || i == 30 || i == 31 
           a(i) = axes('position', [0.6 1.176-(0.01*i) 0.06 0.0045]);
           imagesc(dat(i)); axis off;
           colormap(a(i), cmap)
           caxis(a(i), clim)
        elseif i == 32 || i == 33 || i == 34 
            a(i) = axes('position', [0.6 1.156-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 1.136-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end

    elseif i < 49
        if i == 38 || i == 39 || i == 40 || i == 41 || i == 42
            a(i) = axes('position', [0.6 1.01-(0.008*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        elseif  i == 43 || i == 44 || i == 45
            a(i) = axes('position', [0.6 1.08-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 1.06-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    elseif i < 55
        if i == 49 || i == 50 ||i == 51 
            a(i) = axes('position', [0.6 1.01-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 0.995-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    elseif i < 63
        if i == 55 || i == 56 || i == 57
            a(i) = axes('position', [0.6 0.97-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 0.785-(0.007*i) 0.06 0.004]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    elseif i < 70
        if i == 63 || i == 64 || i == 65 || i == 66
            a(i) = axes('position', [0.6 0.825-(0.008*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 0.943-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    elseif i < 76
        if i == 70 || i == 71 || i == 72
            a(i) = axes('position', [0.6 0.92-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 0.903-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    else
        if i == 76 || i == 77 || i == 78 || i == 79
            a(i) = axes('position', [0.6 0.729-(0.008*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        else
            a(i) = axes('position', [0.6 0.873-(0.01*i) 0.06 0.0045]);
            imagesc(dat(i)); axis off;
            colormap(a(i), cmap)
            caxis(a(i), clim)
        end
    end
end



j = 70;
a(j) = axes('position', [0.9 0.05 0.05 0.9]);
imagesc((flipud([clim(1):0.01:clim(end)]'))); axis off
colormap(a(j), cmap)
caxis(a(j), clim);