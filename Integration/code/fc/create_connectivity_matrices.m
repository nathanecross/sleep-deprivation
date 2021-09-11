


% ### ### ### IMPORTANT !!! ### ### ###
%
% If the data has not been imported previously, first run the script 
% import_data.m
%

NumberOfParcels = 400;
sub_dir = '../../data/Parcels/';
out_dir = '../../data/fc';
load('../../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = length(SubjectName);
nparc = char(string(NumberOfParcels));
keyword = 'task-PVT';
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}];


% Correlation between parcels for each subject
for t = 1:length(times(:,1))
    for z= 1:NumberOfSubjects 
        fprintf(['Correlating timeseries between parcels for ' char(times(t,3)) ' - ' char(SubjectName{z}) '\n']);
        load([ sub_dir char(SubjectName{z}) '/' char(times(t,3)) '/' char(SubjectName{z}) '_' keyword '.mat'])
        for i=1:NumberOfParcels
            for j = 1:NumberOfParcels                 
                [a, pv] = corrcoef(Parcels(i,:),Parcels(j,:));
                ConMat.corr(i,j) = a(1,2);
                ConMat.pval(i,j) = pv(1,2);
            end 
        end
        fprintf( ['Performing Fisher r to z transformations for ' char(times(t,3)) ' - ' char(SubjectName{z}) '\n' ]);
        ConMat.corr_z = atanh(ConMat.corr);
        ConMat.corr_z(ConMat.corr_z>=8)=NaN;
        if ~exist(sprintf('%s/',out_dir), 'dir')
            mkdir(sprintf('%s/',out_dir))
        end
        if ~exist(sprintf('%s/%s',out_dir,char(SubjectName{z})), 'dir')
            mkdir(sprintf('%s/%s',out_dir,char(SubjectName{z})))
        end
        if ~exist(sprintf('%s/%s/%s',out_dir,char(SubjectName{z}),string(times{t,3})), 'dir')
            mkdir(sprintf('%s/%s/%s',out_dir,char(SubjectName{z}),string(times{t,3})))
        end
        % And save the data
        save(sprintf('%s/%s/%s/%s_%s.mat',out_dir,char(SubjectName{z}),string(times{t,3}),char(SubjectName{z}),keyword),'ConMat');
    end
end

% Create group-average correlation matrices per session
holdmat1=NaN(NumberOfSubjects, NumberOfParcels, NumberOfParcels);
holdmat2=NaN(NumberOfSubjects, NumberOfParcels, NumberOfParcels);
for t = 1:length(times(:,1))
    for z = 1:NumberOfSubjects
        load( [out_dir '/' char(SubjectName{z}) '/' char(times(t,3)) '/' char(SubjectName{z}) '_' keyword '.mat']) 
        holdmat1(z,:,:) = ConMat.corr; 
        holdmat2(z,:,:) = ConMat.corr_z; 
    end
end
clear ConMat
for i=1:NumberOfParcels
    for j=1:NumberOfParcels
        ConMat.Session_mean.(char(times(t,3)))(i,j) = nanmean(holdmat1(:,i,j));
        ConMat.Session_meanz.(char(times(t,3)))(i,j) = nanmean(holdmat1(:,i,j));
    end
end


% Save the data
if ~exist([out_dir '_grp/'], 'dir')
          mkdir([out_dir '_grp/']);
end
save([out_dir sprintf('_grp/ConMat%s_%s_%s.mat',nparc,nnetw,keyword)], 'ConMat'); 




