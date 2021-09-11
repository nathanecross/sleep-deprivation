function  Subjectdata = loadGeneralData_gii(filename_lh,filename_rh,sub_dir,session,run,keyword)

    fprintf(['loading ' filename_lh ' + ' filename_rh '\n'])
    filepath_lh = char(strcat(sub_dir, filename_lh));
    filepath_rh = char(strcat(sub_dir, filename_rh));
    lhg = gifti(filepath_lh);
    rhg = gifti(filepath_rh);
    Subjectdata.lh_data = lhg.cdata;
    Subjectdata.rh_data = rhg.cdata;
    
    % Import confounds file (only if importing directly from fmriprep)
    filename2 = find_preproc_file(sub_dir, strcat(keyword, '_', run), 'confounds_regressors.tsv');
    if length(filename2) > 0 
        if startsWith(string(filename2),'.') == 0
            filepath = char(strcat(sub_dir, filename2))
            confounds = tdfread(filepath);
            Subjectdata.global_signal = confounds.global_signal;
            Subjectdata.white_matter = confounds.white_matter;
            Subjectdata.csf = confounds.csf;
            framewise_displacement = confounds.framewise_displacement;
            framewise_displacement(1,:) = "0";
            framewise_displacement = str2num(framewise_displacement);   
            Subjectdata.framewise_displacement = framewise_displacement;
            Subjectdata.trans_x = confounds.trans_x;
            Subjectdata.trans_y = confounds.trans_y;
            Subjectdata.trans_z = confounds.trans_z;
            Subjectdata.rot_x = confounds.rot_x;
            Subjectdata.rot_y = confounds.rot_y;
            Subjectdata.rot_z = confounds.rot_z;
        end
    else
        fprintf(['No confounds file found' '\n'])
    end


end
