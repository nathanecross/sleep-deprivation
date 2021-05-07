% ~ IMPORT DATA ~
% 
% Script to load in fmri data in matlab and save to .mat files for further analyses
%
% times:        Define study-specific session & run combinations, as per named in BIDS directory
%               A MxN cell:
%                   - Column #1: BIDS session name (e.g. 'ses-01')
%                   - Column #2: BIDS run name (e.g. 'run-01')
%                   - Column #3: Study run name (e.g. 'Baseline')
% keyword:      Set to specific sequence - this is case sensitive 
%                 - to concatenate sequences, set keyword1 = 'task-All', and set conditions
% subjectlist   Either 'all' (to read all files in directory), or a cell array of subject ids (e.g. {'sub-02', {'sub-02'}, ...}

%% General options
options.times = [{'ses-01','','Control'} 
                 {'ses-02','run-01','PreNap'}
                 {'ses-02','run-02', 'PostNap'}]; %[{'ses-02','','nap'}];
options.keyword = {'task-ANT', 'task-Nback', 'task-PVT'}; %{'run-01_nap','run-02_nap'};  %{'task-CROSS'};
options.subjectlist = 'all'; %{'sub-XX'}

% Import surface?
options.import.surf.run = 'OFF';    
options.import.surf.inpath = '../../../nifti_bids/derivatives/xcpEngine/output_bids/clean_nogsr_surfsmooth6/';
options.import.surf.template = 'fsaverage5';
options.import.surf.outpath = '../data/Subjectdata_gii';
options.import.surf.concat = 'ON';
options.import.surf.concatname = 'task-All';

% Parcellate surface?
options.parcel.run = 'OFF';
options.parcel.template = {'Schaeffer',400};
options.parcel.outpath = '../data/Parcels';
options.network.template = {'Yeo',17};

% Import volume data? e.g. for subcortical structures
options.import.vol.run = 'ON';     
options.import.vol.inpath = '../../../nifti_bids/derivatives/xcpEngine/output_bids/clean_nogsr_volsmooth6/';
options.import.vol.outpath = '../data/Subjectdata_nii';
options.import.vol.roi = {'amygdala_lh',...
                          'amygdala_rh',...
                          'caudate_lh',...
                          'caudate_rh',...
                          'hippocampus_lh',...
                          'hippocampus_rh',...
                          'pallidum_lh',...
                          'pallidum_rh',...
                          'putamen_lh',...
                          'putamen_rh',...
                          'thalamus_lh',...
                          'thalamus_rh'};
options.import.vol.roipath = '../labels/MNIatlas/subc/func_masks/';
options.import.vol.concat = 'ON';
options.import.vol.concatname = 'task-All';

% Regress confounds of interest?
options.conf.confounds = 'Tasktimes.mat';
options.conf.surf.run = 'ON';
options.conf.surf.inpath = options.parcel.outpath; 
options.conf.surf.outpath = '../data/Parcels_regr';
options.conf.vol.run = 'OFF';
options.conf.vol.inpath = options.import.vol.outpath;
options.conf.vol.outpath = '../data/Subjectdata_nii_regr';



%   ~ RUN ~
addpath(genpath(pwd))
import_pipeline(options)



