function Subject = concatenate_timeseries_vols(file,times)

table = readtable("labels/roi.txt",'Delimiter','\t');
roi = table2array(table(:,1));    

    for f = 1:length(file)
        filename = char(file(f));
        load(sprintf('data/Subject_vols_%s.mat',filename)); 
        NumberOfSubjects = numel(fieldnames(Subject.vols.Control));
        SubjectName = fieldnames(Subject.vols.Control);
        task = extractAfter(filename,'-');
        for t=1:length(times)
            for z= 1:NumberOfSubjects
                if isfield(Subject.vols.(char(times(t,3))),(SubjectName{z})) == 1
                    fprintf( 'Concatenating lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
                    for r=1:numel(roi)
                        if exist('Task') == 0

                                name = sprintf('Task.vols.%s.%s.%s',(char(times(t,3))),(SubjectName{z}),char(roi(r)));
                                Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) = Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r)));
                        else
                            if isfield(Task.vols,(char(times(t,3)))) == 1
                                if isfield(Task.vols.(char(times(t,3))),(SubjectName{z})) == 1
                                    if isfield(Task.vols.(char(times(t,3))).(SubjectName{z}),(char(roi(r)))) == 1
                                        Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) = [Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r)))];
                                    else
                                        Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) = Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r)));
                                    end
                                else
                                    Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) = Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r)));
                                end
                            else
                                Task.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r))) = Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(roi(r)));
                            end
                        end
                    end
                else
                    fprintf( 'No lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
                end
            end
        end
        clear Subject
    end
      
    Subject = Task;
    clear Task
    Subject.NumberOfSubjects = NumberOfSubjects;
    Subject.SubjectName = SubjectName;
    Subject.times = times;
    
    
          
end