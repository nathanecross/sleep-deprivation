% tasks
nparc = '400';
keyword = 'task-All';
nnetw = '17';
load(sprintf('data/Subject_vols_%s.mat', keyword));
load(sprintf('data/Parcels%s_%s_%s.mat',nparc,nnetw,keyword));
load(sprintf('data/ConMat%s_%s_%s.mat',nparc,nnetw,keyword));
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
load('SubjectName.mat');
NumberOfSubjects = 20; NumberOfNets = str2num(nnetw); NumberOfParcels = str2num(nparc);
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap'];
     
for t = 1:length(times)
    for z = 1:NumberOfSubjects
        sprintf('Running Subject %s for %s', (char(SubjectName{z})), (char(times(t,3))))
        vols = fieldnames(Subject.vols.(char(times(t,3))).(char(SubjectName{z})));
        ConMat.times.(char(times(t,3))).(char(SubjectName{z})).corr_sub = ...
            ConMat.times.(char(times(t,3))).(char(SubjectName{z})).corr;
        for v = 1:length(vols)
            for p = 1:NumberOfParcels
                [rr, pp] = corrcoef(Subject.vols.(char(times(t,3))).(char(SubjectName{z})).(char(vols(v))), ...
                   Parcels.(char(times(t,3))).(char(SubjectName{z})).ShfMeanclean.(char("p"+p)));
                submat(t,z,v,p) = rr(2);
                ConMat.times.(char(times(t,3))).(char(SubjectName{z})).corr_sub(NumberOfParcels+v,p) = ...
                    rr(2);
            end
        end
        % z-score
        ConMat.times.(char(times(t,3))).(char(SubjectName{z})).Shfcorr_sub_z = ...
            atanh(ConMat.times.(char(times(t,3))).(char(SubjectName{z})).corr_sub);
        ConMat.times.(char(times(t,3))).(char(SubjectName{z})).Shfcorr_sub_z(ConMat.times.(char(times(t,3))).(char(SubjectName{z})).Shfcorr_sub_z>8)=NaN;
        submat_z(t,z,:,:) = ConMat.times.(char(times(t,3))).(char(SubjectName{z})).Shfcorr_sub_z(401:end,:);
    end
    % session mean
    for i=1:NumberOfParcels+length(vols)
        for j=1:NumberOfParcels
            for z=1:NumberOfSubjects
                if isfield(ConMat.times.(char(times(t,3))),(char(SubjectName{z}))) == 1
                    holdr(1,z) = ConMat.times.(char(times(t,3))).(char(SubjectName{z})).corr_sub(i,j);
                    ConMat.Session_mean_sub.(char(times(t,3)))(i,j) = mean(holdr);
                    holdz(1,z) = ConMat.times.(char(times(t,3))).(char(SubjectName{z})).Shfcorr_sub_z(i,j);
                    ConMat.Session_mean_sub_z.(char(times(t,3)))(i,j) = mean(holdz);
                else
                    continue
                end
            end
        end
    end
end

save(sprintf('data/ConMat%s_%s_%s_vols.mat',nparc,nnetw,keyword), 'ConMat');

%% nap
load(sprintf('data/Subject_nap-Allchunks.mat'));
load(sprintf('data/Parcels%s_%s_nap-Allchunks.mat',nparc,nnetw));
load(sprintf('data/ConMat%s_%s_nap-Allchunks.mat',nparc,nnetw));
load(sprintf('labels/Yeo%s_Shf%s.mat',nnetw,nparc));
load('SubjectName.mat');
NumberOfSubjects = 20; NumberOfNets = str2num(nnetw); NumberOfParcels = str2num(nparc);

for z = 1:NumberOfSubjects
    if z == 9 || z == 20
        submat(4,z,:,:) = NaN;
        submat_z(4,z,:,:) = NaN;
    else 
        sprintf('Running Subject %s', (char(SubjectName{z})))
        vols = fieldnames(Subject.vols.(char(SubjectName{z})));
        ConMat.times.nap.(char(SubjectName{z})).corr_sub = ...
            ConMat.times.nap.(char(SubjectName{z})).corr;
        for v = 1:length(vols)
            for p = 1:NumberOfParcels
                [rr, pp] = corrcoef(Subject.vols.(char(SubjectName{z})).(char(vols(v))), ...
                   Parcels.nap.(char(SubjectName{z})).ShfNet(p,:));
                submat(4,z,v,p) = rr(2);
                ConMat.times.nap.(char(SubjectName{z})).corr_sub(NumberOfParcels+v,p) = ...
                    rr(2);
            end
        end
            % z-score
        ConMat.times.nap.(char(SubjectName{z})).Shfcorr_sub_z = ...
            atanh(ConMat.times.nap.(char(SubjectName{z})).corr_sub);
        ConMat.times.nap.(char(SubjectName{z})).Shfcorr_sub_z(ConMat.times.nap.(char(SubjectName{z})).Shfcorr_sub_z>8)=NaN;
        submat_z(4,z,:,:) = ConMat.times.nap.(char(SubjectName{z})).Shfcorr_sub_z(401:end,:);
    end
end
for i=1:NumberOfParcels+length(vols)
    for j=1:NumberOfParcels
        for z=1:NumberOfSubjects
            if isfield(ConMat.times.nap,(char(SubjectName{z}))) == 1
                holdr(1,z) = ConMat.times.nap.(char(SubjectName{z})).corr_sub(i,j);
                ConMat.Session_mean_sub.nap(i,j) = mean(holdr);
                holdz(1,z) = ConMat.times.nap.(char(SubjectName{z})).Shfcorr_sub_z(i,j);
                ConMat.Session_mean_sub_z.nap(i,j) = mean(holdz);
            else
                continue
            end
        end
    end
end
save(sprintf('data/ConMat%s_%s_nap-Allchunks.mat',nparc,nnetw), 'ConMat');


thal = submat(:,:,3:4,:);
thal = mean(thal(:,:,:,:),4);
thal = mean(thal(:,:,:,:),3);
thal = squeeze(thal)';

thal2 = thal;
thal = [thal2(:,1:2) thal2(:,4) thal2(:,3)] ;
submat = cat(1,submat(1,:,:,:), submat(2,:,:,:), submat(4,:,:,:), submat(3,:,:,:));

save('data/Thalamus.mat', 'thal', 'submat');


%% Create group mean and plot
submat_mean = squeeze(mean(submat,2));
vol_names = vols(1:2:end);
for i = 1:length(vol_names)
    vol_names(i) = cellstr(vol_names{i}(1:end-3));
end
for t = 1:length(times)
    matrix = squeeze(submat_mean(t,:,:));
    for l = 1:length(matrix(:,1))
         new = reorder_vectors(matrix(l,:)',nparc,nnetw);
         matrix(l,:) = new';
    end
    figure;   
    imagesc(matrix);
    set(gca, 'YTickLabel', vol_names, 'FontSize', 5);
    caxis([-1, 1]);
    colormap fireice; % <--- go to town and play around on this 

end



%% Conduct t-ttests
load('data/Thalamus.mat','submat')
clear holdmat
for j=1:NumberOfParcels
    for i = 1:length(vols)
            for z = 1:NumberOfSubjects
                        holdmat(z,1) = submat(1,z,i,j);
                        holdmat(z,2) = submat(2,z,i,j);
                        holdmat(z,3) = submat(3,z,i,j);
                        holdmat(z,4) = submat(4,z,i,j);
            end
            a = array2table(holdmat);
            a.Properties.VariableNames = {'WR','SD','NREM','PRN'}; % add in the co-variates here
            s1 = holdmat(:,1);
            s2 = holdmat(:,2);
            s3 = holdmat(:,3);
            s4 = holdmat(:,4);
            [~,p,~,stats] = ttest(s2, s1);
            SDEP_vols_Tmat(i,j) = stats.tstat;
            SDEP_vols_Pmat(i,j) = p;
            [~,p,~,stats] = ttest(s4, s2);
            Recovery_vols_Tmat(i,j) = stats.tstat;
            Recovery_vols_Pmat(i,j) = p;
            [~,p,~,stats] = ttest(s4, s1);
            ConRec_vols_Tmat(i,j) = stats.tstat;
            ConRec_vols_Pmat(i,j) = p;
            [~,p,~,stats] = ttest(s3, s1);
            NREM_vols_Tmat(i,j) = stats.tstat;
            NREM_vols_Pmat(i,j) = p;

    end 
end


matrix = SDEP_vols_Tmat;
for l = 1:length(matrix(:,1))
         new = reorder_vectors(matrix(l,:)',nparc,nnetw);
         matrix(l,:) = new';
end
figure;   
imagesc(matrix');
set(gca, 'XTickLabel', vol_names, 'FontSize', 5);
caxis([-8, 8]);
colormap fireice; % <--- go to town and play around on this 



%%
load('/Volumes/NATHAN/SDEP/code_and_output/data/Integration400_7.mat')
deltaI = HI.Itot(:,2) - HI.Itot(:,1);
deltaFCR = HI.FCR(:,2) - HI.FCR(:,1);

thal = submat(:,:,3:4,:);

thal_change = (thal(2,:,:,:)-thal(1,:,:,:))./thal(1,:,:,:);
thal_change(thal_change<-50)=NaN;
thal_change_mean = squeeze(nanmean(thal_change,3));
thal_change_mean = mean(thal_change_mean,2);


[r,p] = corrcoef(thal_change_mean,deltaFCR)
figure;
scatter((thal_change_mean),(deltaI))

thal_indi = squeeze(nanmean(thal,3));
thal_indi = mean(thal_indi,3);
plot([thal_indi(1,:)',thal_indi(2,:)',thal_indi(3,:)']')
