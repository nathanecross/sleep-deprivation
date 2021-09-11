function Schaeffer = map_Schaeffer_parcellations(filename,nparc,lh_labels_Shf,rh_labels_Shf)
       
    % ### ### ### IMPORTANT !!! ### ### ###
    %
    % If the data has not been imported previously, first run the script 
    % import_data.m
    %


    load(sprintf('data/%s.mat',filename)); 
    %load('times/TaskTimes.mat');
    NumberOfSubjects = Subject.NumberOfSubjects;
    SubjectName = Subject.SubjectName;
    times = Subject.times;
    task = extractAfter(filename,'-');


    % Combining data into single dataset per subject + applying Schaeffer et al. 2019 parcellations 
    NumberOfParcels = str2num(nparc);
    Shf_lh = logical(zeros(length(lh_labels_Shf),(NumberOfParcels/2)));
    % Left hemisphere timeseries mean
    for t = 1:length(times(:,1))
        for z = 1:NumberOfSubjects
            fprintf( 'Creating lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
            for i = 1:(NumberOfParcels/2)
                for y = 1:length(lh_labels_Shf)
                    Shf_lh(y,i) = lh_labels_Shf(y,1)==i;
                end
                if isfield(Subject.(char(times(t,3))),(char(SubjectName{z}))) == 1
                Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNet_lh.(char("p"+num2str(i))) =  Subject.(char(times(t,3))).(char(SubjectName{z})).lh_data(Shf_lh(:,i),:);
                Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean_lh.(char("p"+num2str(i))) = mean(Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNet_lh.(char("p"+num2str(i))));
                else
                    fprintf( 'lh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + "doesn't exist \n");
                    continue
                end
            end
        end 
    end

    % Right hemisphere timeseries mean
    Shf_rh = logical(zeros(length(rh_labels_Shf),(NumberOfParcels-((NumberOfParcels/2)))));
    for t =1:length(times(:,1))
        for z= 1:NumberOfSubjects
            fprintf( 'Creating rh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n');
            for i =(NumberOfParcels/2)+1:NumberOfParcels
                for y = 1:length(rh_labels_Shf)
                    Shf_rh(y,i-50) = rh_labels_Shf(y,1) == i;
                end
                if isfield(Subject.(char(times(t,3))),(char(SubjectName{z}))) == 1
                Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNet_rh.(char("p"+num2str(i))) =  Subject.(char(times(t,3))).(char(SubjectName{z})).rh_data(Shf_rh(:,i-50),:);
                Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean_rh.(char("p"+num2str(i))) = mean(Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNet_rh.(char("p"+num2str(i))));
                else
                    fprintf( 'rh timeseries mean for ' + string(times(t,3)) + ' - ' + SubjectName{z} + "doesn't exist \n");
                    continue
                end
            end
        end
    end

    % Concatenate Left and Right hemispheres into ShfNetMean
    for t =1:length(times(:,1))
        for z= 1:NumberOfSubjects
            fprintf( 'Concatenate Left and Right hemispheres for ' + string(times(t,3)) + ' - ' + SubjectName{z} + '\n')
            if isfield(Subject.(char(times(t,3))),(char(SubjectName{z}))) == 1
                for i = 1:(NumberOfParcels/2)
                     Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean.(char("p"+num2str(i))) = Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean_lh.(char("p"+num2str(i)));
                end
                for i = (NumberOfParcels/2)+1:NumberOfParcels
                     Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean.(char("p"+num2str(i))) = Schaeffer.times.(char(times(t,3))).(char(SubjectName{z})).ShfNetMean_rh.(char("p"+num2str(i)));
                end
            else
                continue
            end
        end
    end

end
