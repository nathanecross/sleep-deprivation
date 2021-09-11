
% /Volumes/NATHAN/SDEP/EEG_tasks

load('ica/psd_cat.mat')
chans = table2cell(readtable('channels.csv'));

mask = ~startsWith(chans(:,2),'E');

coi = psd(:,:,mask,:);

coi_ave = squeeze(mean(coi,3));


figure;
for p = 1:19
    plot(squeeze(coi_ave(p,1,:))); hold on
end

coi_ave(coi_ave(:,1,2)>1000,:,:) = NaN;
coi_ave(coi_ave(:,1,4)>1000,:,:) = NaN;


for f = 1:10
    
    NN = squeeze(coi_ave(:,1,f));
    SD = squeeze(coi_ave(:,2,f));
    
    delta_f(:,f) = SD-NN;
    
    [~,pv(f),~,stats] = ttest(SD,NN);
    
    tstat(f) = stats.tstat;
    
end


%% Relationships with Integration
load('Integration400_7_task-All_stdinf.mat')

delta_f(20,:) = NaN;
FCR_delta = HI.FCR(:,2)- HI.FCR(:,1);
for f = 1:10
    [r,p] = corrcoef(FCR_delta,delta_f(:,f),'rows','complete');
    cor(f) = r(2);
    pf(f) = p(2);
end

I_delta = HI.Itot(:,2)- HI.Itot(:,1);
for f = 1:10
    [r,p] = corrcoef(I_delta,delta_f(:,f),'rows','complete');
    cor(f) = r(2);
    pf(f) = p(2);
end


