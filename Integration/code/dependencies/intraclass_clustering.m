function ShfLabel = intraclass_clustering(matrix,nparc,nnetw,plot)

% Clusters parcels within each of a set number of networks based on 
% intraclass correlation from a connectivity matrix
% Inputs: 
% matrix - Correlation matrix. Each matrix should be thresholded at p<0.05.
% nparc - No. of parcellations
% nnetw - No. of networks


load(sprintf('labels/Yeo%snames.mat',nnetw))
Yeomask = load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc)); 
Yeomask = Yeomask.(char(string(fieldnames(Yeomask))));
ShfLabel = (1:str2double(nparc))';
 % group average correlation matrix

figure;
for n = 1:str2double(nnetw)
    % mask network
    mask = Yeomask==(n);
    tempmat = matrix(mask,mask);
    tempmat = tempmat.^2; %take R-squared
    labels = [(1:length(tempmat))' ShfLabel(mask)];
    % calculate intraclass similarity
    class = sqrt(1-tempmat);
    D = linkage(class,'ward','euclidean');
    nn = 1:length(tempmat);
    np = length(nn);
    % threshold by greatest increase of intracluster distance
    for i = 1:length(D)
        d = get_cluster_dists(D, i);
        dist(i,1) = mean(d);
    end
    cutoff = max(dist);
    % take only clusters below threshold
    TT = cluster(D,'Cutoff',cutoff,'Criterion','distance');
    labels(:,3) = str2num([repmat([char(string(Networks(n,1))) char("0")],length(labels),1) char(string(TT))]);
    for p = 1:length(labels)
        ShfLabel(labels(p,2),2) = labels(p,3);
    end
    
    if contains(plot,'on')
        % plot clusters
        if n < 5
            pos1 = 4;
        elseif n < 9
            pos1 = 3;
        elseif n <13
            pos1 = 2;
        else
            pos1 = 1;
        end
        if n == 1 || n == 5 || n == 9 || n == 13
            pos2 = 1;
        elseif n == 2 || n == 6 || n == 10 || n == 14
            pos2 = 2;
        elseif n == 3 || n == 7 || n == 11 || n == 15
            pos2 = 3;
        elseif n == 4 || n == 8 || n == 12 || n == 16
            pos2 = 4;
        else
            pos2 = 5;
        end 
        a(n) = axes('position', [0.17*pos2-0.05 0.22*pos1-0.15 0.15 0.15]);
        h.(sprintf('tree%s',string(n))) = dendrogram(D,'ColorThreshold', cutoff);
        thistree = h.(sprintf('tree%s',string(n)));
            for l = 1:length(thistree)
                if sum(thistree(l).Color == [0 1 1]) == 3
                    thistree(l).Color = [0 0.06 1];
                elseif sum(thistree(l).Color == [0.5 1 0]) == 3 || sum(thistree(l).Color == [0 1 0]) == 3 || sum(round(thistree(l).Color,1) == [0.8 1.0 0.0]) == 3 
                    thistree(l).Color = [0 0.49 0.25];
                elseif sum(thistree(l).Color == [0 1 1]) == 3 || sum(round(thistree(l).Color,1) == [0.0 1.0 0.4]) == 3
                    thistree(l).Color =  [0.89 0.59 0];
                end
            end

        title(sprintf('%s',char(string(Networks(n,2)))))
        clear dist 
    else
    end
end

