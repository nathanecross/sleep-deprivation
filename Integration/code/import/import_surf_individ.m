function subject = import_surf_individ(inpath, outpath, keyword, times, part)
%
% Imports surface data and saves it as a .mat file for further analyses.
%
% INPUT
% 
% inpath         path where preprocessed data is stored
%                WARNING: Data must be BIDS compatible!
%
% outpath        path to save output
%
% keyword        should be set to look for specific sequences 
%                (eg. keyword = 'task-CROSS') 
%
% times          should be a nx3 cell of study specific session, run and name
%                combinations 
%                (e.g. times = [{'ses-01','run-01','Control'}])
%
%
% OUTPUT         
%                Each surface is saved as a .mat file in << outpath >>
%
% subject        a struct containing:
%                                    - SubjectList:         List of subjects in <<< inpath >>
%                                    - FuncFoldergii:       List of subjects with surface data
%                                    - noFuncFoldergii:     List of subjects without surface data
%                                    - NumberOfSubjects:    # of subjects with surface data
%                                    - times:               Same as specified in input
%
% Author: Nathan Cross (2019)
% ##########################################################################

% Check keyword type and update if necessary
try key = split(keyword,"-",1); catch; key=keyword; end
if length(key)>1
    key = char(key(2));
else
    key = char(key);
end

FuncFoldergii = {}; noFuncFoldergii = {};
for t=1:length({times{:,1}})
    session=times{t,1};
    run=times{t,2};
    ext = '.gii';
   
    % Obtain total number of subjects, subjects with data, and subjects with missing data
    homeFolder = pwd;
    if strcmp(part,'all')
        [SubjectList, NumberOfSubjects] = cycle_directory_contents(inpath);
    else
        SubjectList = part;
        NumberOfSubjects = length(part);
    end
    
    cd(homeFolder)
    for n = 1:NumberOfSubjects
        skip = false;
        if exist(strcat(inpath, char((SubjectList(n))), "/", session, "/func/")) == 0
            noFuncFoldergii(t,n) = SubjectList(n); %creates list of subjects with missing data
            NumberOfSubjects=NumberOfSubjects-1;
        else
            FuncFoldergii(t,n) = SubjectList(n); %creates list of subjects with data
%             subjectNumber = split(SubjectList(n,1),"-"); 
%             SubjectName = "sub"+subjectNumber(2,1);
            SubjectName = char(SubjectList(n));
            path = strcat(inpath, SubjectName,  "/", session, "/func/");
            filename = find_preproc_file(path, key, ext);
            pat = {'-L','lh'};
            filename_lh_all = filename(find(contains(filename, pat)));
            pat = {'-R','rh'};
            filename_rh_all = filename(find(contains(filename, pat)));
            
            if length(filename_lh_all)>1
                index_lh = contains(filename_lh_all, run);
                filename_lh = char(filename_lh_all(index_lh));
            elseif isempty(filename_lh_all) == 1
                sprintf('WARNING: Identifier for left and right hemispheres is incorrect for %s',char(SubjectName))
                skip = true;
            else
                filename_lh = char(filename_lh_all);
            end
            if length(filename_rh_all)>1
                index_rh = contains(filename_rh_all, run);
                filename_rh = char(filename_rh_all(index_rh));
            else
                filename_rh = char(filename_rh_all);
            end
            
            if skip ~= 1
                subjectdata = loadGeneralData_gii(filename_lh,filename_rh,path,session,run,keyword);
                subjectdata.SubjectName = SubjectName;

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
                save(sprintf('%s/%s/%s/%s_%s.mat',outpath,SubjectName,string(times{t,3}),SubjectName,keyword),'subjectdata');
            end
        end
    end

end
    subject.FuncFoldergii = FuncFoldergii; 
    subject.noFuncFoldergii = noFuncFoldergii;
    subject.SubjectList = SubjectList; 
    subject.NumberOfSubjects = NumberOfSubjects;
    subject.times = times;
end

