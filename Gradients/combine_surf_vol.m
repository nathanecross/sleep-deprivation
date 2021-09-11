function Global_timeseries = combine_surf_vol(keyword)


    load((sprintf('Schaeffer400_%s.mat',keyword)));
    load((sprintf('Subject_vols_%s.mat',keyword)));
    NumberOfSubjects = Subject.NumberOfSubjects;
    SubjectName = Subject.SubjectName;
    times = Subject.times;
    
    for t=1:length(times)
            for z= setdiff(1:NumberOfSubjects, 4) % FIX THIS WHEN SUB-06 Nback fixed
                for p = 1:numel(fields(Schaeffer.times.(char(times(t,3))).(SubjectName{z}).ShfNetMean))
                    if p==1
                        Global_timeseries.(char(times(t,3))).(SubjectName{z}).data = Schaeffer.times.(char(times(t,3))).(SubjectName{z}).ShfNetMean.(char("p"+p));
                    else
                        Global_timeseries.(char(times(t,3))).(SubjectName{z}).data = [Global_timeseries.(char(times(t,3))).(SubjectName{z}).data; Schaeffer.times.(char(times(t,3))).(SubjectName{z}).ShfNetMean.(char("p"+p))];               
                    end
                end
                region = (fields(Subject.vols.(char(times(t,3))).(SubjectName{z})));
                for r = 1:numel(region)
                    Global_timeseries.(char(times(t,3))).(SubjectName{z}).data = [Global_timeseries.(char(times(t,3))).(SubjectName{z}).data; Subject.vols.(char(times(t,3))).(SubjectName{z}).(char(region(r)))];
                end
            end
    end
end
