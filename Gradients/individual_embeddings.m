
function indi_embed = individual_embeddings(subjdata, embedding)
    
    states = fieldnames(subjdata.times);
    for s = 1:3
        this_state = subjdata.times.(states{s});
        subs = fieldnames(this_state);
        for ii = 1:length(subs)
            mat         = this_state.(subs{ii}).ShfCorr_z;
            mat(isnan(mat)) = 1;
            normangle   = connectivity2normangle(mat, 90);
            [indi_embed(:,:,ii,s), results] = mica_diffusionEmbedding(normangle);

            % align gradients to group-level embedding
            [U,~,V] = svd(embedding(:,:,1).' * indi_embed(:,:,ii,s),0);
            U = U(:, 1:min(size(U,2), size(V,2)));
            V = V(:, 1:min(size(U,2), size(V,2)));
            xfms = V * U';
            indi_embed(:,:,ii,s) = indi_embed(:,:,ii,s) * xfms;
        end
    end
end