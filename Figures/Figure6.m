%
% Figure 6
%

%% Load dependencies and setup
addpath('gradients');
addpath('dependencies/customcolormap');
load('../code_output_gsr/gradients/indi_embed.mat');
load('data/Behave.mat');
load('SubjectName.mat');
times = {'Control','PreNap','PostNap'};
%% EXTENT OF EACH GRADIENT AT BASELINE AND RELATIONSHIP TO BEHAVIOUR

for g = 1:3
    for z = 1:length(SubjectName)
        for t = 1:length(times)
            siz(g,z,t) = max(indi_embed(:,g,z,t)) - min(indi_embed(:,g,z,t));
            bin_val(t,g,:,z) = slidingWindowEqualDistrib(indi_embed(:,g,z,t)', indi_embed(:,g,z,t)', 50, 10);
            siz_bin(g,z,t) = bin_val(t,g,50,z) - bin_val(t,g,1,z);
        end
    end
    
end


%% Relationship to behaviour

% BEHAVIOUR
Behave_nn = Behave((contains(Behave.SessionID,'NN')==1),:);
Behave_pre = Behave((contains(Behave.SessionID,'SDpre')==1),:);
Behave_post = Behave((contains(Behave.SessionID,'SDpost')==1),:);
Nback = Behave_pre.x2back_correct_resp_percent - Behave_nn.x2back_correct_resp_percent;
Nback_rec = Behave_post.x2back_correct_resp_percent - Behave_pre.x2back_correct_resp_percent;

for Figure6 = 1
siz_combine = squeeze(siz(1,:,1)) + squeeze(siz(1,:,3));
I = customcolormap([0 0.25 0.5 0.75 1], {'#00AA03', '#50BF70', '#bababa', '#FF9B9B', '#A80505'});
f = figure;
    a(1) = axes('position', [0.05 0.5 0.4 0.4]);
    scatter(Nback, siz_combine, 20, Nback, 'filled', 'MarkerEdgeColor', 'k');
    colormap(a(1), I)
    xlim([min(Nback) max(Nback)])
    ylim([0.38 0.55])
    l = lsline;
    l.Color = 'k';
    l.LineWidth = 2;
    l.LineStyle = '--';
    a(2) = axes('position', [0.55 0.5 0.4 0.4]);
    scatter(Nback_rec, siz_combine, 20, Nback_rec, 'filled', 'MarkerEdgeColor', 'k');
    colormap(a(2), I)
    xlim([min(Nback_rec) max(Nback_rec)])
    ylim([0.38 0.55])
    l = lsline;
    l.Color = 'k';
    l.LineWidth = 2;
    l.LineStyle = '--';
    fname = char('../code_output_gsr/gradients/F7_corr.png');
    exportfig(f, fname, 'Format', 'png', 'Resolution', 600, 'FontSize', 0.8)
end