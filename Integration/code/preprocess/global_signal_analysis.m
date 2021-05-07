%% Import global signal 
sub_dir='../../data/Subjectdata_gii_cat/';
keyword = 'task-All'; %<-------------------------------- Set to specific task - this is case sensitive (for all tasks: keyword1 = 'task-All')
format = '.gii';  %<--- Set to format in which you would like to import (e.g. volumetric = '.nii'; surface = '.gii'; both = ['.nii','.gii'])
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02'}, 'PostNap']; 
load('../../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = length(SubjectName);

%% Tasks (cortex only)
for t = 1:length(times)
    for z = 1:NumberOfSubjects
        load([sub_dir char(SubjectName{z}) '/' (char(times(t,3))) '/' ...
              char(SubjectName{z}) '_' keyword '.mat']);
        sub = erase((char(SubjectName{z})),'-');
        gs.(char(times(t,3))).(sub) = mean([subjectdata.lh_data; subjectdata.rh_data]);
        gsf(z,t) = std(gs.(char(times(t,3))).(sub));
    end
end

%% Nap
nap_dir = '../../data/Subjectdata_gii_cat/';
keyword = 'nap';
for z = 1:NumberOfSubjects
    if z == 9  
    elseif z == 20
    else
    sub = erase((char(SubjectName{z})),'-');
    load([nap_dir char(SubjectName{z}) '/nap/' ...
          char(SubjectName{z}) '_' keyword '.mat']);
    gs.nap.(sub) = mean([subjectdata.lh_data; subjectdata.rh_data]);  
    gsf(z,4) = nanstd(gs.nap.(sub));
    end
end

GS.gs = gs; GS.gsf = gsf;
save('../../regressors/GS.mat', 'GS');

%% Tasks (cortex + subcortex)
clear gs gsf GS
vol_dir = '../../data/Subjectdata_nii_cat/';
keyword = 'task-All';
for t = 1:length(times)
    for z = 1:NumberOfSubjects
        sub = erase((char(SubjectName{z})),'-');
        load([sub_dir char(SubjectName{z}) '/' (char(times(t,3))) '/' ...
              char(SubjectName{z}) '_' keyword '.mat']);
        all = [subjectdata.lh_data; subjectdata.rh_data];  
        load([vol_dir char(SubjectName{z}) '/' (char(times(t,3))) '/' ...
              char(SubjectName{z}) '_' keyword '.mat']); 
        volumes = fieldnames(subjectdata);
        for f = 1:length(volumes)
            if contains(volumes(f), 'SubjectName') == 0
                all = [all; subjectdata.(char(volumes(f)))];
            end
        end
        gs.(char(times(t,3))).sub = mean(all);
        gsf(z,t) = std(gs.(char(times(t,3))).sub);
    end
end

%% Nap (cortex + subcortex)
nap_dir = '../../data/Subjectdata_gii_cat/';
vol_dir = '../../data/Subjectdata_nii_cat/';
keyword = 'nap';
for z = 1:NumberOfSubjects
    if z == 9  
    elseif z == 20
    else
    sub = erase((char(SubjectName{z})),'-');    
    load([nap_dir char(SubjectName{z}) '/nap/' ...
          char(SubjectName{z}) '_' keyword '.mat']);
    all = [subjectdata.lh_data; subjectdata.rh_data];  
    load([vol_dir char(SubjectName{z}) '/nap/' ...
          char(SubjectName{z}) '_' keyword '.mat']); 
    volumes = fieldnames(subjectdata);
    for f = 1:length(volumes)
        if contains(volumes(f), 'SubjectName') == 0
            all = [all; subjectdata.(char(volumes(f)))];
        end
    end
    gs.nap.sub = mean(all);
    gsf(z,4) = nanstd(gs.nap.sub);
    end
end


GSa.gs = gs; GSa.gsf = gsf;
save('../../regressors/GSa.mat', 'GSa');

