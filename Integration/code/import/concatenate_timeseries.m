function concatenate_timeseries(subject, inpath, outpath, times, conditions, format, outname)

% 
%   
%

SubjectList = subject.SubjectList; 
NumberOfSubjects = subject.NumberOfSubjects;

% Concatenate cortical surface data
if contains(format,'gii') == 1 
    for z= 1:NumberOfSubjects
        SubjectName = char(SubjectList{z});
        for t=1:length(times(:,1))
            fprintf( ['Concatenating surface timeseries mean for ' SubjectName ' - ' char(string(times(t,3))) '\n']);
            subject.surf.lh_data = []; subject.surf.rh_data = [];         
            for f = 1:length(conditions)
                filename = [SubjectName '_' char(conditions(f))];
                try
                    load(sprintf('%s/%s/%s/%s.mat',inpath,SubjectName,char(times(t,3)),filename)); 
                    subject.surf.lh_data = [subject.surf.lh_data subjectdata.lh_data];
                    subject.surf.rh_data = [subject.surf.rh_data subjectdata.rh_data]; 
                catch
                    continue
                end
            end
            % Lets make some output directories
            if ~exist(sprintf('%s/',outpath), 'dir')
               mkdir(sprintf('%s/',outpath))
            end
            if ~exist(sprintf('%s/%s',outpath,SubjectName), 'dir')
               mkdir(sprintf('%s/%s',outpath,SubjectName))
            end
            if ~exist(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})), 'dir')
               mkdir(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})))
            end
            % And save the data
            subjectdata.lh_data = subject.surf.lh_data;
            subjectdata.rh_data = subject.surf.rh_data;
            subjectdata.SubjectName = SubjectName;
            save(sprintf('%s/%s/%s/%s_%s.mat',outpath,SubjectName,string(times{t,3}),SubjectName,outname),'subjectdata');
        end
    end
end

% Concatenate parcels
if contains(format,'Parcels') == 1
    for z= 1:NumberOfSubjects
        SubjectName = char(SubjectList{z});
        for t=1:length(times(:,1))
            fprintf( ['Concatenating surface timeseries mean for ' SubjectName ' - ' char(string(times(t,3))) '\n']);
            Parcels_cat = [];         
            for f = 1:length(conditions)
                filename = [SubjectName '_' char(conditions(f))];
                try
                    load(sprintf('%s/%s/%s/%s.mat',inpath,SubjectName,char(times(t,3)),filename)); 
                    Parcels_cat = [Parcels_cat Parcels];
                catch
                    continue
                end
            end
            Parcels = Parcels_cat;
            
            % Lets make some output directories
            if ~exist(sprintf('%s/',outpath), 'dir')
               mkdir(sprintf('%s/',outpath))
            end
            if ~exist(sprintf('%s/%s',outpath,SubjectName), 'dir')
               mkdir(sprintf('%s/%s',outpath,SubjectName))
            end
            if ~exist(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})), 'dir')
               mkdir(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})))
            end
            
            % And save the data
            save(sprintf('%s/%s/%s/%s_%s.mat',outpath,SubjectName,string(times{t,3}),...
                 SubjectName,outname),'Parcels');
        end
    end
end

% Concatenate subcortical volume data
if contains(format,'nii') == 1 
    for z= 1:NumberOfSubjects
        SubjectName = char(SubjectList{z});
        for t=1:length(times(:,1))
            fprintf( ['Concatenating volume timeseries mean for ' SubjectName ' - ' char(string(times(t,3))) '\n']);       
            for c = 1:length(conditions)
                filename = [SubjectName '_' char(conditions(c))];
                try
                    load(sprintf('%s/%s/%s/%s.mat',inpath,SubjectName,char(times(t,3)),filename)); 
                    fields = fieldnames(subjectdata);
                    for f = 1:length(fields)
                        if contains(fields(f),'SubjectName')==0
                            if c==1
                                subject.vol.(char(fields(f))) = subjectdata.(char(fields(f)));
                            else
                                subject.vol.(char(fields(f))) = [subject.vol.(char(fields(f))) subjectdata.(char(fields(f)))];
                            end
                        end
                    end
                catch
                    continue
                end
            end
            % Lets make some output directories
            if ~exist(sprintf('%s/',outpath), 'dir')
               mkdir(sprintf('%s/',outpath))
            end
            if ~exist(sprintf('%s/%s',outpath,SubjectName), 'dir')
               mkdir(sprintf('%s/%s',outpath,SubjectName))
            end
            if ~exist(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})), 'dir')
               mkdir(sprintf('%s/%s/%s',outpath,SubjectName,string(times{t,3})))
            end
            % And save the data
            subjectdata = subject.vol;
            save(sprintf('%s/%s/%s/%s_%s.mat',outpath,SubjectName,string(times{t,3}),SubjectName,outname),'subjectdata');
        end
    end
end  
          
end