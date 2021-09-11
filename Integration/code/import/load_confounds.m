function  [Subjectdata,NumberOfSubjects,SubjectName,noFuncFolder] = load_confounds(sub_dir,session,run,condition,keyword,ext)
homeFolder = pwd;
[SubjectList, NumberOfFolders] = cycle_directory_contents(sub_dir);
noFuncFolder = []; FuncFolder = []; remove = [];
NumberOfSubjects=NumberOfFolders;
cd(homeFolder)
for n = 1:NumberOfFolders
    if exist(strcat(sub_dir, char((SubjectList(n))), "/", session, "/func/")) == 0
        noFuncFolder = [noFuncFolder SubjectList(n,:)]; %creates list of subjects with missing data
        NumberOfSubjects=NumberOfSubjects-1;
    else
        FuncFolder = [FuncFolder SubjectList(n,:)];
    end
end
FuncFolder = FuncFolder';
SubjectName = cell(length(FuncFolder),1);
for o = 1:length(FuncFolder)%%% change 
    SubjectNumber = split(FuncFolder(o,1),"-"); 
    SubjectName{o} = "sub"+SubjectNumber(2,1);
    path = strcat(sub_dir, (FuncFolder(o)),  "/", session, "/func/");
    filename = find_preproc_file(path, keyword, ext);
    filename_lh_all = filename(find(contains(filename, 'lh')));
    filename_rh_all = filename(find(contains(filename, 'rh')));
    if length(filename_lh_all)>1
        index_lh = contains(filename_lh_all, run);
        filename_lh = char(filename_lh_all(index_lh));
    else 
        filename_lh = char(filename_lh_all);
    end
    if length(filename_rh_all)>1
        index_rh = contains(filename_rh_all, run);
        filename_rh = char(filename_rh_all(index_rh));
    else
        filename_rh = char(filename_rh_all);
    end
    
    if isempty(filename) == 0
        filename2 = find_preproc_file(path, strcat(keyword, '_', run), 'confounds_regressors.tsv');
        filename2 = filename2(startsWith(filename2,'.')==0);
        if length(filename2) > 0
            sub = char(SubjectName{o});
            %fprintf(['Loading confounds file for ' char(SubjectName{o}) '\n'])
            filepath = char(strcat(sub_dir, (FuncFolder(o)),  "/", session, "/func/", filename2));
            confounds = tdfread(filepath);
            Subjectdata.(condition).(sub).global_signal = confounds.global_signal';
            Subjectdata.(condition).(sub).white_matter = confounds.white_matter';
            Subjectdata.(condition).(sub).csf = confounds.csf';
            framewise_displacement = confounds.framewise_displacement;
            framewise_displacement(1,:) = "0";
            framewise_displacement = str2num(framewise_displacement);   
            Subjectdata.(condition).(sub).framewise_displacement = framewise_displacement';
            Subjectdata.(condition).(sub).trans_x = confounds.trans_x';
            Subjectdata.(condition).(sub).trans_y = confounds.trans_y';
            Subjectdata.(condition).(sub).trans_z = confounds.trans_z';
            Subjectdata.(condition).(sub).rot_x = confounds.rot_x';
            Subjectdata.(condition).(sub).rot_y = confounds.rot_y';
            Subjectdata.(condition).(sub).rot_z = confounds.rot_z';
        else
            fprintf(['No confounds file for ' char(SubjectName{o}) '\n'])
            continue
        end
    else
        clear filename_lh_all filename_rh_all
    end
end
 cd(homeFolder);
end
