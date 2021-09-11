function import_pipeline(options)

% Runs Preliminary Analysis Pipeline for fMRI data:
% 1. Imports surface data
% 2. Applies parcellation template to surface data
% 3. Imports volumetric data
% 4. (Optional) Regress out any confound of interest (e.g. motion, task events)
% 5. Creates connectivty matrices
%
% All functions save separate outputs in BIDS format
%     to location specified in options.<step>.outpath

addpath(genpath(pwd));

% 1. Import surface data
if contains(options.import.surf.run,'ON')
    inpath = options.import.surf.inpath; 
    keyword = options.keyword; 
    times = options.times;  
    outpath = options.import.surf.outpath;
    part = options.subjectlist;
    for c = 1:length(keyword)
            subject = import_surf_individ(inpath, outpath, char(keyword(c)), times, part);
    end
    
    if contains(options.import.surf.concat,'ON') == 1
        outname = options.import.surf.concatname;
        concatenate_timeseries(subject, outpath, [outpath '_cat'], ...
                               times, options.keyword, 'gii',outname);
    end
    save('../../data/subject.mat','subject');
end

% 2. Parcellation of surface data (so far just Schaeffer parcels)
if contains(options.parcel.run,'ON') 
    ptemplate = char(string(options.parcel.template(1)));
    NumberOfParcels = cell2mat(options.parcel.template(2));
    nparc = char(string(options.parcel.template(2)));
    ntemplate = char(string(options.network.template(1)));
    NumberOfNetworks = cell2mat(options.network.template(2));
    nnetw = char(string(options.network.template(2)));
    fs = options.import.surf.template;
    
    % Load Parcellation Template
    load(sprintf('../../labels/%s/%s/Labels%s_%s.mat',fs,ptemplate,nparc,nnetw))
    if contains(ptemplate,'Schaeffer')==1
        lh_labels = lh_labels_Shf;
        rh_labels = rh_labels_Shf;
    end
    
    % Map Schaeffer parcellations to the Yeo Networks
    if contains(ptemplate,'Schaeffer')==1 && contains(ntemplate,'Yeo')==1 
        Yeo_Clusters_ref = load(sprintf('../../labels/%s/%sNetworks/1000subjects_clusters%s_ref.mat',...
                                    fs,ntemplate,nnetw));
        for i = 1:length(lh_labels_Shf)
            lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
            lh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.lh_labels(i);
        end
        for i = 1:length(rh_labels_Shf)
            rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
            rh_labels_Shf_Yeo(i,2) = Yeo_Clusters_ref.rh_labels(i);
        end     
    end
    
    inpath = options.import.surf.outpath; % Check inpath/outpath based on concatentation or not
    outpath = options.parcel.outpath;
    ops.nparc=nparc; ops.nnetw=nnetw; ops.fs=fs; ops.labpath='labels';
    ops.template = ptemplate; ops.times = options.times;
     
    
    for c = 1:length(options.keyword)
        keyword = options.keyword(c);
        try key = split(keyword,"-",1); catch; key=keyword; end
        if length(key)>1
            key = char(key(2));
        else
            key = char(key);
        end
        ops.key = {key};
        map_parcellations(inpath, outpath, 'all', ops, lh_labels, rh_labels)
    end
    
    if contains(options.import.surf.concat,'ON') == 1
        inpath = [inpath '_cat'];
        outpath = [options.parcel.outpath '_cat'];
        try
            ops = rmfield( ops , 'key' );
        catch
        end
        map_parcellations(inpath, outpath, 'all', ops, lh_labels, rh_labels)
    end
end

% 3. Import volumetric (i.e. subcortical) data
if contains(options.import.vol.run,'ON')
    inpath = options.import.vol.inpath; 
    keyword = options.keyword; 
    times = options.times;  
    outpath = options.import.vol.outpath;
    roi = options.import.vol.roi;
    roipath = options.import.vol.roipath;
    for c = 1:length(keyword)
            subject = import_vol_individ(inpath, outpath, char(keyword(c)), times, roi, roipath);
    end
    if contains(options.import.vol.concat,'ON') == 1
       outname = options.import.vol.concatname;
       concatenate_timeseries(subject, outpath, [outpath '_cat'], ...
                               times, options.keyword, 'nii', outname);
    end
    save('../../data/subject.mat','subject');
end

% Regress out confounds of interest
if contains(options.conf.surf.run,'ON')
    load('../../data/subject.mat')
    confounds = options.conf.confounds;
    conditions = options.keyword;
    inpath = options.conf.surf.inpath;
    outpath = options.conf.surf.outpath;
    times = options.times;
    confound_regression(subject, inpath, outpath, confounds, times, conditions, 'gii');
    
    if contains(options.import.surf.concat,'ON') == 1
       outname = options.import.surf.concatname;
       concatenate_timeseries(subject, outpath, [outpath '_cat'], ...
                               times, options.keyword, 'Parcels', outname);
    end
end
if contains(options.conf.vol.run,'ON')
    load('../../data/subject.mat')
    confounds = options.conf.confounds;
    conditions = options.keyword;
    times = options.times; 
    inpath = options.conf.vol.inpath;
    outpath = options.conf.vol.outpath;
    confound_regression(subject, inpath, outpath, confounds, times, conditions, 'nii');
    
    if contains(options.import.vol.concat,'ON') == 1
       outname = options.import.vol.concatname;
       concatenate_timeseries(subject, outpath, [outpath '_cat'], ...
                               times, options.keyword, 'nii', outname);
    end
end

end