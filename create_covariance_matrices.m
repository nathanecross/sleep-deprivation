function CovMat = create_coraviance_matrices(filename,nparc,lh_labels_Shf,rh_labels_Shf)

% ### ### ### IMPORTANT !!! ### ### ###
%
% If the data has not been imported previously, first run the script 
% import_data.m
%

load(sprintf('data/%s.mat',filename)); 
%load('times/TaskTimes.mat');
NumberOfSubjects = length(fieldnames(Schaeffer.times.Control));
SubjectName = fieldnames(Schaeffer.times.Control);
times = fieldnames(Schaeffer.times);
task = extractAfter(filename,'-');
NumberOfParcels = str2num(nparc);


% Correlation between parcels for each subject
for t =1:length(times)
    for z= 1:NumberOfSubjects 
        fprintf( 'Calculating covariance between parcels for ' + string(times(t)) + ' - ' + SubjectName{z} + '\n');
        if isfield(Schaeffer.times.(char(times(t))),(SubjectName{z})) == 1
            for x=1:NumberOfParcels
                for y = 1:NumberOfParcels                 
                    a = cov(Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfMeanclean.("p"+num2str(x)),Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfMeanclean.("p"+num2str(y)));
                    CovMat.times.(char(times(t))).(SubjectName{z}).cov(x,y) = a(1,2);
                end 
            end
        else
            continue
        end
    end
end


% Create mean of z matrices per session
for t =1:length(times)
    holdz=zeros(1,NumberOfSubjects);
    CovMat.Session_mean.(char(times(t))) = zeros(NumberOfParcels);
    fprintf( 'Creating mean for ' + string(times(t)) + '\n')
    for i=1:NumberOfParcels
        for j=1:NumberOfParcels
            for z=1:NumberOfSubjects
                if isfield(CovMat.times.(char(times(t))),(SubjectName{z})) == 1
                    holdz(1,z) = CovMat.times.(char(times(t))).(SubjectName{z}).cov(i,j);
                    CovMat.Session_mean.(char(times(t)))(i,j) = mean(holdz);
                else
                    continue
                end
            end
        end
    end
end