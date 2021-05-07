
% Make Schaefer Parcels from annotations file (downloaded from:
% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3)

NumberofParcels='300';
NumberofNetworks='17';

% Left Hemisphere
filename = sprintf('annot%s/lh.Schaefer2018_%sParcels_%sNetworks_order.annot', NumberofParcels,NumberofParcels,NumberofNetworks);
[v, lab, ct] = read_annotation(filename,1);
n = length(lab);

for i=1:length(lab)
   annot = lab(i); 
   index = find(ct.table(:,5) == annot);
   lh_labels_Shf(i,1) = index - 1; 
end

% Right Hemisphere
filename = sprintf('annot%s/rh.Schaefer2018_%sParcels_%sNetworks_order.annot', NumberofParcels,NumberofParcels, NumberofNetworks);
[v, lab, ct] = read_annotation(filename);
n = length(lab);

for i=1:length(lab)
   annot = lab(i); 
   index = find(ct.table(:,5) == annot);
   rh_labels_Shf(i,1) = index + (str2double(NumberofParcels)/2); 
end

rh_labels_Shf(rh_labels_Shf==50)=0;
save(sprintf('ShfLabels%s_%s.mat',NumberofParcels,NumberofNetworks),'lh_labels_Shf','rh_labels_Shf')