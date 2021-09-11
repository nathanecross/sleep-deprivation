%% Plotting changes for Louvain communities

load('data/Integration400_7_task-All_louv.mat')
%within each network
for n = 1:6
    if n<5
        Iwsn(n,:,1) = HI.Control{21, 1}.int_intra.samples(n,n,:);
        Iwsn(n,:,2) = HI.PreNap{21, 1}.int_intra.samples(n,n,:);
        sig(n+1,1) = sum(squeeze(Iwsn(n,:,1))>squeeze(Iwsn(n,:,2)))/1000;
    elseif n ==6
        Iwsn(n-1,:,1) = HI.Control{21, 1}.int_intra.samples(n,n,:);
        Iwsn(n-1,:,2) = HI.PreNap{21, 1}.int_intra.samples(n,n,:);
        sig(n,1) = sum(squeeze(Iwsn(n-1,:,1))>squeeze(Iwsn(n-1,:,2)))/1000;
    end 
end

Iws_WR = squeeze(mean(Iwsn(:,:,1),2));
Iws_SD = squeeze(mean(Iwsn(:,:,2),2));

Iws_ch = ((HI.PreNap{21,1}.int_total.mean - HI.Control{21,1}.int_total.mean)/HI.PreNap{21,1}.int_total.mean)*100;
Iwsn_ch = ((Iws_SD-Iws_WR)./Iws_WR)*100;
Iws_ch = [Iws_ch; NaN(length(Iwsn_ch)-length(Iws_ch),1)]


data = [Iws_ch,Iwsn_ch]
cmap = customcolormap([0 0.2 0.4 0.5 0.8 0.9 1], {'#6b1200','#E23603','#FF8D33',  '#bababa', '#336EFF', '#336EFF', '#033c9e'});

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
clim = [-50 50];

f = figure;

% cortex
a(1) = axes('position', [0.1 0.05 0.1 0.9]);
imagesc(dat(1)); axis off;
colormap(a(1), cmap)
caxis(a(1), clim)

% communities 
for i = 2:6
        a(i) = axes('position', [0.3 1.1-(0.17*i) 0.1 0.15]);
        imagesc(dat(i)); axis off;
        colormap(a(i), cmap)
        caxis(a(i), clim)
end

% colorbar
j = 7;
a(j) = axes('position', [0.9 0.05 0.05 0.9]);
imagesc((flipud([clim(1):0.01:clim(end)]'))); axis off
colormap(a(j), cmap)
caxis(a(j), clim);