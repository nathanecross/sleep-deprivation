
function ConMat = create_connectivity_matrices(filename,nparc,lh_labels_Shf,rh_labels_Shf)

% ### ### ### IMPORTANT !!! ### ### ###
%
% If the data has not been imported previously, first run the script 
% import_data.m
%

load(sprintf('data/%s.mat',filename)); 
%load('times/TaskTimes.mat');
times = fieldnames(Schaeffer.times);
NumberOfSubjects = length(fieldnames(Schaeffer.times.(char(string(times(1))))));
SubjectName = fieldnames(Schaeffer.times.(char(string(times(1)))));

task = extractAfter(filename,'-');
NumberOfParcels = str2num(nparc);


% Correlation between parcels for each subject
for t = 1:length(times)
    for z= 1:NumberOfSubjects 
        fprintf( 'Correlating timeseries between parcels for ' + string(times(t)) + ' - ' + SubjectName{z} + '\n');
        if isfield(Schaeffer.times.(char(times(t))),(SubjectName{z})) == 1
            for x=1:NumberOfParcels
                for y = 1:NumberOfParcels                 
                    a = corrcoef(Schaeffer.times.(char(times(t))).(char(SubjectName{z})).ShfNetMean.(char("p"+num2str(x))),Schaeffer.times.(char(times(t))).(char(SubjectName{z})).ShfNetMean.(char("p"+num2str(y))));
                    ConMat.times.(char(times(t))).(char(SubjectName{z})).corr(x,y) = a(1,2);
                end 
            end
        else
            continue
        end
    end
end

% Apply Fisher A to Z transformation on connectivity matrix
for t =1:length(times)
    for z= 1:NumberOfSubjects 
        fprintf( 'Performing Fisher r to z transformations for ' + string(times(t)) + ' - ' + SubjectName{z} + '\n');
        if isfield(Schaeffer.times.(char(times(t))),(char(SubjectName{z}))) == 1 
            ConMat.times.(char(times(t))).(char(SubjectName{z})).ShfCorr_z = atanh(ConMat.times.(char(times(t))).(char(SubjectName{z})).corr);
            ConMat.times.(char(times(t))).(char(SubjectName{z})).ShfCorr_z(ConMat.times.(char(times(t))).(char(SubjectName{z})).ShfCorr_z>=8)=NaN;
        else
           continue
        end
    end
end


% Create mean of matrices per session
for t =1:length(times)
    holdr=zeros(1,NumberOfSubjects);
    holdz=zeros(1,NumberOfSubjects);
    ConMat.Session_mean.(char(times(t))) = zeros(NumberOfParcels);
    ConMat.Session_mean_z.(char(times(t))) = zeros(NumberOfParcels);
    fprintf( 'Creating mean for ' + string(times(t)) + '\n')
    for i=1:NumberOfParcels
        for j=1:NumberOfParcels
            for z=1:NumberOfSubjects
                if isfield(ConMat.times.(char(times(t))),(char(SubjectName{z}))) == 1
                    holdr(1,z) = ConMat.times.(char(times(t))).(char(SubjectName{z})).corr(i,j);
                    ConMat.Session_mean.(char(times(t)))(i,j) = mean(holdr);
                    holdz(1,z) = ConMat.times.(char(times(t))).(char(SubjectName{z})).ShfCorr_z(i,j);
                    ConMat.Session_mean_z.(char(times(t)))(i,j) = mean(holdz);
                else
                    continue
                end
            end
        end
    end
end


% Save data
% keyword1 = char(extractAfter(filename,"_"));



