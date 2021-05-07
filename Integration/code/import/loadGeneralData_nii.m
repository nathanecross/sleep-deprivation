function  Subjectdata = loadGeneralData_nii(path,filename,session,run,key,roi,roipath)
    
    file = [char(path) char(filename)];
    nifti = load_nii(file);  

    for r = 1:numel(roi)
        regionOfInterest = char(roi(r));
        sprintf('Masking ROI: %s',regionOfInterest)
        mask = load_mask(roipath, regionOfInterest);
        %b = repmat(reshape(mask, size(mask,1), size(mask,2), 1, size(mask,3)), 1, 1, size(mask,2), 1);
        for tt = 1:numel(nifti.img(1,1,1,:))
            scan = nifti.img(:,:,:,tt);
            timeseries(tt) = mean(scan(mask)); 
        end
        Subjectdata.(sprintf('%s',char(regionOfInterest))) = timeseries;
    end

    % % Import confounds file (only if importing directly from fmriprep)
    filename2 = find_preproc_file(path, strcat(key, '_', run), 'confounds_regressors.tsv');
    if length(filename2) > 0
        filepath = char(strcat(sub_dir, (FuncFolder(o)),  "/", session, "/func/", filename2));
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
    else
        fprintf('No confounds file found \n')
    end
    
end
