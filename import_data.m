function Subject = import_data(subject_dir, keyword, format, times)
%
% IMPORT DATA - Imports the data and also saves it as a .mat file
% Set subject_dir to be the path where preprocessed data is stored
%                           WARNING: Data must be BIDS compatible!
% Arguments in should be set to look for specific sequences 
% (e.g. keyword = 'task-CROSS') as well as whether the data is in 
% the volumetric (e.g. format = '.nii';) or surface space 
% (e.g. format = '.gii').
% Times should be a nx3 cell of study specific session, run and name
% combinations (e.g. times = [{'ses-01','run-01','Control'}]). 
% 
z = 1;

if contains(format,'gii') == 1 %Load GIFT file (.nii) for cortical surface data
    for i=1:length({times{:,1}})
        session=times{i,1};
        run=times{i,2};
        ext = '.gii';
        sub_dir = subject_dir{z};
        [Subjectdata,NumberOfSubjects,SubjectName,noFuncFolder] = loadGeneralData_gii(sub_dir,session,run,keyword,ext);
        Time = times{i,3};
        Subject.(Time) = Subjectdata;
        clear i Time run session Subjectdata
        Subject.NumberOfSubjects = NumberOfSubjects;
        Subject.SubjectName = SubjectName;
        Subject.times = times;
        Subject.Missing_data = noFuncFolder;   
    end
    z = 2;
end

if contains(format,'nii') == 1   %Load NIFTI file (.nii) for other regions of interest
    ext = '.nii';
    sub_dir = subject_dir{z};
    table = readtable("labels/roi.txt",'Delimiter','\t');
    roi = table2array(table(:,1));  
    for i=1:length({times{:,1}})
        session=times{i,1};
        run=times{i,2};
        [Subjectdata] = loadGeneralData_nii(sub_dir,session,run,keyword,ext,roi);
        Time = times{i,3};
        Subject.vols.(Time) = Subjectdata;
        
        clear i Time run session Subjectdata  
    end
    

    
    
    

    % table = readtable(regionsOfInterestFile,'Delimiter','\t'); %% change if needed
    % brainStruct = table(:,1); brainStruct = table2array(brainStruct);
    % fields = {'srcbext','analyzehdr','bhdr','niftihdr', ...
    %     'fspec','pwd','flip_angle','tr','te','ti',...
    %     'vox2ras0','volsize','height','width','depth',...
    %     'nframes','vox2ras','nvoxels','xsize','ysize',...
    %     'zsize','x_r','x_a','x_s','y_r','y_a','y_s',...
    %     'z_r','z_a','z_s','c_r','c_a','c_s','vox2ras1',...
    %     'Mdc','volres','tkrvox2ras','vol'};
end

if contains(format,'gii') == 0 && contains(format,'nii') == 0
    fprintf('Incorrect format specification. Check if format is specified to either ".nii" or ".gii" \n');   
  
end

