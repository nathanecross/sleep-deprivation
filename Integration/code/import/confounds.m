sub_dir = '../../../nifti_bids/derivatives/fmriprep/';
format = '.gii';
conditions = {'task-ANT', 'task-Nback', 'task-PVT'};
times = [{'ses-01','','Control'} 
         {'ses-02','run-01','PreNap'}
         {'ses-02','run-02', 'PostNap'}];

% Load confounds for tasks
for c = 1:3     
    for t = 1:3
        [Subjectdata,NumberOfSubjects,SubjectName,noFuncFolder] = ...
            load_confounds(sub_dir,times(t,1),times(t,2),char(times(t,3)),char(conditions(c)),format);
        Subject.(char(times(t,3))) = Subjectdata.(char(times(t,3)));
    end
end

% Extract framewise displacement
for t = 1:length(times)
    for z = 1:NumberOfSubjects
        FD.(char(times(t,3))).(char(SubjectName{z})) = ...
            [Subject.ANT.surf.(char(times(t,3))).(char(SubjectName{z})).framewise_displacement' ...
            Subject.Nback.surf.(char(times(t,3))).(char(SubjectName{z})).framewise_displacement' ...
            Subject.PVT.surf.(char(times(t,3))).(char(SubjectName{z})).framewise_displacement'];
        FD.tab(z,t) = mean(FD.(char(times(t,3))).(char(SubjectName{z})));
    end
end

FDtab(1,1) = nanmean(FD.tab(:,1))
FDtab(1,2) = nanstd(FD.tab(:,1))
FDtab(1,3) = nanmean(FD.tab(:,2))
FDtab(1,4) = nanstd(FD.tab(:,2))

% Differences between sessions
[p,tbl,FD.stats] = anova1(FD.tab)

save('../regressors/FD.mat','FD');

%% Resting state

% Load confounds for tasks
for t = 1:3
    [Subjectdata,NumberOfSubjects,SubjectName,noFuncFolder] = load_confounds(sub_dir,times(t,1),times(t,2),char(times(t,3)),'task-CROSS',format);
    Subjectrs.(char(times(t,3))) = Subjectdata.(char(times(t,3)));
end

% Extract framewise displacement
for t = 1:length(times)
    
    for z = 1:NumberOfSubjects
        if isfield(Subjectrs.(char(times(t,3))),(char(SubjectName{z})))
            FDrs.(char(times(t,3))).(char(SubjectName{z})) = Subjectrs.(char(times(t,3))).(char(SubjectName{z})).framewise_displacement';
            FDrs.tab(z,t) = mean(FDrs.(char(times(t,3))).(char(SubjectName{z})));
        else
            FDrs.tab(z,t) = NaN
        end
    end
end

FDtab(2,1) = nanmean(FDrs.tab(:,1))
FDtab(2,2) = nanstd(FDrs.tab(:,1))
FDtab(2,3) = nanmean(FDrs.tab(:,2))
FDtab(2,4) = nanstd(FDrs.tab(:,2))

% Differences between sessions
[p,tbl,FDrs.stats] = anova1(FDrs.tab)

save('../regressors/FDrs.mat','FDrs');

