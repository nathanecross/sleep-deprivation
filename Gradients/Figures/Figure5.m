%
% Figure 5
%

%% Load dependencies and setup
addpath('../../dependencies');
NOGSR = load('../../data/ConMat400_task-All.mat');
GSR = load('../../data/ConMat400_task-All.mat');
load('../../data/GLM.mat'); NOGSR.GLM = GLM;
load('../../data/gsr/GLM.mat'); GSR.GLM = GLM;
load('../../gradients/embedding.mat'); NOGSR.embedding = embedding;
load('../../gradients/gsr/embedding.mat'); GSR.embedding = embedding;
load('../../labels/Yeo17_Shf400.mat')
times = fieldnames(GSR.ConMat.times);

%% Figure 5A + B
f = figure('units', 'centimeters', 'position', [0 0 20 40]);
cmap = fireice;
cmap = [1 1 1; cmap];
for figure5a = 1
    for t = 1:2
        a(1) = axes('position', [-0.2+t*0.23 0.8 0.2 0.11]);
        matrix = NOGSR.ConMat.Session_mean_z.(char(times(3)));
        [~, matrix_bordered] = reorder_matrices(matrix, '400', '17'); 
        matrix_plot = tril(matrix_bordered); matrix_plot(matrix_plot==0)=NaN;
        L = matrix_bordered(2,2:end);   
        imagesc(matrix_plot); 
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        caxis([-1, 1]);
        colormap(cmap); 

        a(2) = axes('position', [0.3+t*0.23 0.8 0.2 0.11]);
        matrix = GSR.ConMat.Session_mean_z.(char(times(3)));
        [~, matrix_bordered] = reorder_matrices(matrix, '400', '17'); 
        matrix_plot = tril(matrix_bordered); matrix_plot(matrix_plot==0)=NaN;
        L = matrix_bordered(2,2:end);   
        imagesc(matrix_plot); 
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        caxis([-1, 1]);
        colormap(cmap); 

    end

    matrix = NOGSR.GLM.SDEP_Tmat;
    [~, matrix_bordered] = reorder_matrices(matrix, '400', '17'); 
    L = matrix_bordered(2,2:end);
    matrix_plot = tril(matrix_bordered); matrix_plot(matrix_plot==0)=NaN;
    a(1) = axes('position', [0.09 0.61 0.4 0.17]);
    imagesc(matrix_plot); 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([-5 5]);
    colormap(fireice); 

    matrix = GSR.GLM.SDEP_Tmat;
    [~, matrix_bordered] = reorder_matrices(matrix, '400', '17'); 
    matrix_plot = tril(matrix_bordered); matrix_plot(matrix_plot==0)=NaN;
    a(4) = axes('position', [0.58 0.61 0.3 0.17]);
    imagesc(matrix_plot); 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([-5 5]);
    colormap(cmap); 
end
for Figure5B = 1
    load('Yeo18cmap.mat');
    a(5) = axes('position', [0.1 0.38 0.33 0.18]);
    x = squeeze(NOGSR.embedding(:,1,1));
    y = squeeze(NOGSR.embedding(:,2,1));
    z = squeeze(NOGSR.embedding(:,3,1));
    s = 60;
    c = Yeo17_Shf400;
    h = scatter3(x,y,z,s,c,'filled','LineWidth', 0.00000000000000000001, 'MarkerEdgeColor', 'k');
    set(gca,'XLim',[-0.12 0.16],'YLim',[-0.09 0.15],'ZLim',[-0.1 0.1])
    xlabel('G1','FontSize',12,'FontWeight','bold')
    ylabel('G2','FontSize',12,'FontWeight','bold')
    zlabel('G3','FontSize',12,'FontWeight','bold')
    cmap = Yeo17;
    colormap(gca, cmap)
    view(-6,5)

    a(6) = axes('position', [0.58 0.38 0.33 0.18]);
    figure;
    x = squeeze(GSR.embedding(:,1,1))*-1;
    y = squeeze(GSR.embedding(:,2,1))*-1;
    z = squeeze(GSR.embedding(:,3,1))*-1;
    s = 60;
    c = Yeo17_Shf400;
    h = scatter3(x,y,z,s,c,'filled', 'LineWidth', 0.00000000000000000001, 'MarkerEdgeColor', 'k');
    set(gca,'XLim',[-0.12 0.16],'YLim',[-0.09 0.15],'ZLim',[-0.1 0.1])
    xlabel('G1','FontSize',12,'FontWeight','bold')
    ylabel('G2','FontSize',12,'FontWeight','bold')
    zlabel('G3','FontSize',12,'FontWeight','bold')
    cmap = Yeo17;
    colormap(gca, cmap)
    view(-6,5)

    a(4) = axes('position', [0.46 0.38 0.03 0.18]);
    imagesc((1:100)'); axis off
    colormap(a(4), Yeo17)
end
%% Figure 5C
figure;
% Without GSR
% Pairwise distance in gradient space
for t = 1:3
        for p = 1:400
            grad_dist(p,:,t) = DNorm2(NOGSR.embedding(p,1:3,t) - NOGSR.embedding(:,1:3,t), 2);
        end
end

gd_change = grad_dist(:,:,2) - grad_dist(:,:,1);
I = customcolormap([0 0.2 0.5 0.8 1], {'#D53517','#FBE837', '#DDFEB2' , '#7EE8FA','#0C32F8'});
matrix = gd_change;
[~, matrix_bordered] = reorder_matrices(matrix, '400', '17');
matrix_plot = tril(matrix_bordered); 
a(5) = axes('position', [0.1 0.55 0.3 0.395]);
imagesc(matrix_plot); 
colormap(I);
caxis([-0.04 0.04]);

% With GSR
% Pairwise distance in gradient space
for t = 1:3
        for p = 1:400
            grad_dist(p,:,t) = DNorm2(GSR.embedding(p,1:3,t) - GSR.embedding(:,1:3,t), 2);
        end
end

gd_change = grad_dist(:,:,2) - grad_dist(:,:,1);
matrix = gd_change;
[~, matrix_bordered] = reorder_matrices(matrix, '400', '17');
matrix_plot = tril(matrix_bordered); 
a(5) = axes('position', [0.58 0.55 0.3 0.395]);
imagesc(matrix_plot); 
colormap(I);
caxis([-0.04 0.04]);
