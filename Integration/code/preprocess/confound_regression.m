function confound_regression(subject, inpath, outpath, confounds, times, conditions, format)

%  Function to regress out confound or other timeseries of interest (eg. task onsets)
%  
%  This script assumes that confound timeseries have been formatted previously to the
%  format
%
SubjectList = subject.SubjectList; 
NumberOfSubjects = subject.NumberOfSubjects;

% Load task times
loader = load(sprintf('../regressors/%s', confounds));
regressor = loader.(char(fieldnames(loader)));

% Tweak some features of the regressor
opts = statset('glmfit');
opts.MaxIter = 10000;

% Regress out confounds (e.g. task effects)
for t = 1:length(times)
    for z = 1:NumberOfSubjects 
        SubjectName = char(SubjectList(z));
        for c = 1:length(conditions)
            task = extractAfter(conditions(c),'-');
            load([inpath,'/', SubjectName, '/', char(times(t,3)), '/', SubjectName, ...
                '_', char(conditions(c)) '.mat']);

            if contains(format,'gii') == 1
                    fprintf( 'Regressing out confounds for ' + string(times(t)) + ' - ' + SubjectName + ' - ' + task + '\n');
                    for p = 1:length(Parcels)
                            controls = regressor.(sprintf('%s_%s_conv',char(task),char(times(t,3))))(z,:);
                            [b,dev,stats] = glmfit(controls(1:length(Parcels(p,:))), Parcels(p,:));
                            Parcels(p,:) = stats.resid';
                    end
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
                    save([outpath,'/', SubjectName, '/', char(times(t,3)), '/', SubjectName, '_', char(conditions(c)) '.mat'], ...
                         'Parcels');
                    
            elseif contains(format,'nii') == 1
                rois = fieldnames(subjectdata); rois(contains(rois,'SubjectName')) = [];
                for r = 1:length(rois)
                    if isfield(subjectdata,(char(rois(r)))) == 1
                        controls = regressor.(sprintf('%s_%s_conv',char(task),char(times(t,3))))(z,:);
                        [b,dev,stats] = glmfit(controls(1:length(subjectdata.(char(rois(r))))), subjectdata.(char(rois(r))));
                        subjectdata.(char(rois(r))) = stats.resid';
                    end
                end
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
                save([outpath,'/', SubjectName, '/', char(times(t,3)), '/', SubjectName, '_', char(conditions(c)) '.mat'], ...
                     'subjectdata');
                if exist([inpath '_cat'], 'dir')
                        if c < 2
                            subjectdata_cat = subjectdata;
                        elseif c < length(conditions)
                            subjectdata_cat = [subjectdata_cat subjectdata];
                        else
                            subjectdata_cat = [subjectdata_cat subjectdata]; 
                        end
                        subjectdata = subjectdata_cat;
                        % Lets make some output directories
                        if ~exist(sprintf('%s/',[outpath '_cat']), 'dir')
                           mkdir(sprintf('%s/',[outpath '_cat']))
                        end
                        if ~exist(sprintf('%s/%s',[outpath '_cat'],SubjectName), 'dir')
                           mkdir(sprintf('%s/%s',[outpath '_cat'],SubjectName))
                        end
                        if ~exist(sprintf('%s/%s/%s',[outpath '_cat'],SubjectName,string(times{t,3})), 'dir')
                           mkdir(sprintf('%s/%s/%s',[outpath '_cat'],SubjectName,string(times{t,3})))
                        end
                        save([[outpath '_cat'],'/', SubjectName, '/', char(times(t,3)), '/', SubjectName, '_', char(conditions(c)) '.mat'], ...
                         'subjectdata');
                end
            else
                print('')
            end
        end
    end
end