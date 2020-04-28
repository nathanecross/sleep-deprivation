function Subject = concatenate_timeseries(file)

% 
%   
%
    for f = 1:length(file)
        filename = char(file(f));
        load(sprintf('data/Subject_%s.mat',filename)); 
        NumberOfSubjects = Subject.NumberOfSubjects;
        SubjectName = Subject.SubjectName;
        times = Subject.times;
        task = extractAfter(filename,'-');
        for t=1:length(times)
            for z= 1:NumberOfSubjects
                if isfield(Subject.(char(times(t,3))),(SubjectName{z})) == 1
                    fprintf( 'Concatenating lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
                    name = sprintf('Task.%s.%s.lh_data',(char(times(t,3))),(SubjectName{z}));
                    if exist('Task') == 0
                        Task.(char(times(t,3))).(SubjectName{z}).lh_data = Subject.(char(times(t,3))).(SubjectName{z}).lh_data;
                    else
                        if isfield(Task,(char(times(t,3)))) == 1
                            if isfield(Task.(char(times(t,3))),(SubjectName{z})) == 1
                                if isfield(Task.(char(times(t,3))).(SubjectName{z}),'lh_data') == 1
                                    Task.(char(times(t,3))).(SubjectName{z}).lh_data = [Task.(char(times(t,3))).(SubjectName{z}).lh_data Subject.(char(times(t,3))).(SubjectName{z}).lh_data];
                                else
                                    Task.(char(times(t,3))).(SubjectName{z}).lh_data = Subject.(char(times(t,3))).(SubjectName{z}).lh_data;
                                end
                            else
                                Task.(char(times(t,3))).(SubjectName{z}).lh_data = Subject.(char(times(t,3))).(SubjectName{z}).lh_data;
                            end
                        else
                            Task.(char(times(t,3))).(SubjectName{z}).lh_data = Subject.(char(times(t,3))).(SubjectName{z}).lh_data;
                        end
                    end
                else
                    fprintf( 'No lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
                end
            end
        end
        clear Subject
    end
    
    for f = 1:length(file)
        filename = char(file(f));
        load(sprintf('data/Subject_%s.mat',filename)); 
        NumberOfSubjects = Subject.NumberOfSubjects;
        SubjectName = Subject.SubjectName;
        times = Subject.times;
        task = extractAfter(filename,'-');
        for t=1:length(times)
            for z= 1:NumberOfSubjects
                if isfield(Subject.(char(times(t,3))),(SubjectName{z})) == 1
                    fprintf( 'Concatenating rh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
                    name = sprintf('Task.%s.%s.rh_data',(char(times(t,3))),(SubjectName{z}));
                    if exist('Task') == 0
                        Task.(char(times(t,3))).(SubjectName{z}).rh_data = Subject.(char(times(t,3))).(SubjectName{z}).rh_data;
                    else
                        if isfield(Task,(char(times(t,3)))) == 1
                            if isfield(Task.(char(times(t,3))),(SubjectName{z})) == 1
                                if isfield(Task.(char(times(t,3))).(SubjectName{z}),'rh_data') == 1
                                    Task.(char(times(t,3))).(SubjectName{z}).rh_data = [Task.(char(times(t,3))).(SubjectName{z}).rh_data Subject.(char(times(t,3))).(SubjectName{z}).rh_data];
                                else
                                    Task.(char(times(t,3))).(SubjectName{z}).rh_data = Subject.(char(times(t,3))).(SubjectName{z}).rh_data;
                                end
                            else
                                Task.(char(times(t,3))).(SubjectName{z}).rh_data = Subject.(char(times(t,3))).(SubjectName{z}).rh_data;
                            end
                        else
                            Task.(char(times(t,3))).(SubjectName{z}).rh_data = Subject.(char(times(t,3))).(SubjectName{z}).rh_data;
                        end
                    end
                else
                    fprintf( 'No rh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
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