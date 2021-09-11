function Schaeffer = task_onset_regression(filename,nparc)

load(sprintf('data/%s.mat',filename));
NumberOfSubjects = length(fieldnames(Schaeffer.times.Control));
SubjectName = fieldnames(Schaeffer.times.Control);
times = fieldnames(Schaeffer.times);
task = extractAfter(filename,'-');
NumberOfParcels = str2double(nparc);

% Load task times
load('times/Tasktimes.mat');

% Regress out confounds (e.g. task effects)
for t =1:length(times)
    for z= 1:NumberOfSubjects 
        if isfield(Schaeffer.times.(char(times(t))),(SubjectName{z})) == 1
            fprintf( 'Regressing out confounds for ' + string(times(t)) + ' - ' + SubjectName{z} + '\n');
            for x=1:NumberOfParcels
                    %controls = [Subject.(char(times(t))).(SubjectName{z}).csf Subject.(char(times(t))).(SubjectName{z}).white_matter Subject.(char(times(t))).(SubjectName{z}).trans_x Subject.(char(times(t))).(SubjectName{z}).trans_y Subject.(char(times(t))).(SubjectName{z}).trans_z Subject.(char(times(t))).(SubjectName{z}).rot_x Subject.(char(times(t))).(SubjectName{z}).rot_y Subject.(char(times(t))).(SubjectName{z}).rot_z];
                    controls = Tasktimes.(sprintf('%s_%s_conv',task,string(times(t))))(z,:);
                    [b,dev,stats] = glmfit(controls(1:length(Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfNetMean.(char(("p"+num2str(x)))))), Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfNetMean.(char(("p"+num2str(x)))));
                    Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfMeanclean.(char(("p"+num2str(x)))) = stats.resid';
                    
            end
        else
            fprintf( 'Sequence for ' + string(times(t)) + ' - ' + SubjectName{z} + "doesn't exist \n");
        end
    end
end
