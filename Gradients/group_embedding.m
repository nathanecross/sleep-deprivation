    
function embedding = group_gradient(ConMat, times)
    % organise connectivity matrices
    state_mats = cat(3, ConMat.Session_mean_z.(char(times(1,3))),  ConMat.Session_mean_z.(char(times(2,3))), ConMat.Session_mean_z.(char(times(3,3))));
    state_mats(isnan(state_mats)) = 1; % change diagonal from nan to 1
    
    % perform embedding
    for ii = 1:3
        normangle = connectivity2normangle(state_mats(:,:,ii), 90);
        [embedding(:,:,ii), results] = mica_diffusionEmbedding(normangle);
        if ii == 1
            embedding(:,2,ii) = embedding(:,2,ii).* -1; % flip G2 so DMN is at the top
            lambdas = results.lambdas/sum(results.lambdas);
        else
            % align other states to control
            [U,~,V] = svd(embedding(:,:,1).' * embedding(:,:,ii),0);
            U = U(:, 1:min(size(U,2), size(V,2)));
            V = V(:, 1:min(size(U,2), size(V,2)));
            xfms = V * U';
            embedding(:,:,ii) = embedding(:,:,ii) * xfms;  
        end
    end
end