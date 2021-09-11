function map_parcellations(inpath, outpath, subjlist, ops, lh_labels, rh_labels)
    % 
    % This script maps surface data onto the a parcellation template.
    %
    %
    % INPUT
    % 
    % inpath        path to input data containing .mat files. The input directory should be in BIDS format. 
    %               Each .mat file should contain 2 sets of arrays of size vertices x timepoints which 
    %               correspond to lh_data and rh_data surface data (i.e. imported from function import_data).
    %               
    % outpath       path to save output
    %
    % subjist       list of subjects to run. 
    %               This can either be:
    %                                  a. 'all'
    %                                  b. a cell array, where each cell contains subject ID
    %                                                           
    % ops           structure with the mapping options:
    %               nparc           - number of parcellations to use (e.g. 100)
    %               nnetw           - number of networks (i.e. 7 or 17)
    %               fs              - resolution of freesurfer (e.g. fsaverage5, default)
    %               run (optional)  - cell array containing specific run to analyse (e.g. 'run-01')
    %               key (optional)  - cell array containing specific sequence to analyse (e.g. 'CROSS')
    % 
    % lh_labels     a number of (vertices x 1) array of labels for each vertex on the left hemisphere surface
    %
    % rh_labels     a number of (vertices x 1) array of labels for each vertex on the right hemisphere surface
    %
    %
    %
    % OUTPUT
    %
    % Data is saved in << outpath >> under a BIDS format: ie. each subject's data in their own folder. 
    % 
    %
    % Author: Nathan Cross (2019)
    % ##########################################################################
    
    
    % Check subjlist and prepare list of subjects to run
    if ischar(subjlist) == 1 && strcmp(subjlist, 'all') == 1
        subjlist = dir(inpath);
        subjlist = {subjlist.name}'; 
        subjlist = subjlist(startsWith(string(subjlist),'.')==0); %remove system files
    end
    NumberOfSubjects = length(subjlist);
    
    
    % Set up the correct template
    NumberOfParcels = str2num(ops.nparc);
    parc_lh = logical(zeros(length(lh_labels),(NumberOfParcels/2)));
    parc_rh = logical(zeros(length(rh_labels),(NumberOfParcels/2)));
    
    if ischar(ops.times) == 1 && strcmp(ops.times, 'all') == 1 
            times = dir([inpath '/' char(subjlist{z})]);
            times = {times.name}'; times = times(startsWith(string(times),'.')==0);
    else
        times = ops.times(:,3);
    end
    
    for z = 1:NumberOfSubjects
        % Read session names inside subject directory
        
        for t = 1:length(times) %cycle through sessions
            sesslist = dir([inpath '/' char(subjlist{z}) '/' char(times(t))]);
            sesslist = {sesslist.name}';
            % Set files (runs) to load
            if isfield(ops, 'run') == 1
                runs = ops.run;
            else
                runs = sesslist(startsWith(string(sesslist),'.')==0);
            end
            
            if isfield(ops,'key') == 1
                runs = runs(contains(runs,ops.key));
            end
            
            for r = 1:length(runs) % cycle through files
                run = sesslist((contains(sesslist,char(runs(r)))));
                
                load([inpath '/' char(subjlist{z}) '/' char(times(t)) '/' char(run)]) % load correct file

                % Left hemisphere timeseries mean
                if isfield(subjectdata, 'lh_data') == 1
                    lh_mat = zeros(NumberOfParcels/2,size(subjectdata.lh_data,2));
                    fprintf( 'Performing parcellation of lh for ' + string(times(t)) + ' - ' + run + '\n');
                    for i = 1:(NumberOfParcels/2)
                        for y = 1:length(lh_labels)
                            parc_lh(y,i) = lh_labels(y,1)==i;
                        end
                        parcel_verts =  subjectdata.lh_data(parc_lh(:,i),:);
                        Parcels(i,:) = mean(parcel_verts);
                    end
                else
                    fprintf( 'lh timeseries mean for ' + string(times(t)) + ' - ' + subjlist{z} + "doesn't exist \n");
                    continue
                end

                % Right hemisphere timeseries mean
                if isfield(subjectdata, 'rh_data') == 1
                    rh_mat = zeros(NumberOfParcels/2,size(subjectdata.rh_data,2));
                    fprintf( 'Performing parcellation of rh for ' + string(times(t)) + ' - ' + run + '\n');
                    for i = (NumberOfParcels/2)+1:NumberOfParcels
                        for y = 1:length(rh_labels)
                            parc_rh(y,i-(NumberOfParcels/2)) = rh_labels(y,1)==i;
                        end
                        parcel_verts =  subjectdata.rh_data(parc_rh(:,i-(NumberOfParcels/2)),:);
                        Parcels(i,:) = mean(parcel_verts);
                    end
                else
                    fprintf( 'rh timeseries mean for ' + string(times(t)) + ' - ' + subjlist{z} + "doesn't exist \n");
                    continue
                end
                
                
                % Lets make some output directories
                if ~exist(outpath, 'dir')
                   mkdir(outpath)
                end
                if ~exist([outpath '/' char(subjlist{z})], 'dir')
                   mkdir([outpath '/' char(subjlist{z})])
                end
                if ~exist([outpath '/' char(subjlist{z}) '/' char(times(t))], 'dir')
                   mkdir([outpath '/' char(subjlist{z}) '/' char(times(t))])
                end
                
                % And now save the data
                if contains(run,char(subjlist{z})) == 1
                    save(sprintf('%s/%s/%s/%s', outpath, char(subjlist{z}), char(times(t)), char(run)), 'Parcels')
                else
                    save(sprintf('%s/%s/%s/%s_%s.mat', outpath, char(subjlist{z}), char(times(t)), char(subjlist{z}), char(run)), 'Parcels')
                end
                clear Parcels
            end
        end 
    end
end