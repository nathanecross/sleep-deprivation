%   
%                                                         ,aaaaaa.
%                                                      ,aaabbbbbbaaa. 
%                    _ _         _                  ,aaabbbccccccbbbaaa.
%   __ _ _ _ __ _ __| (_)___ _ _| |_ ___          ,aabbbcccddddddcccbbbaa.  
%  / _` | '_/ _` / _` | / -_) ' \  _(_-<         ,aabbcccddd'    'dddcccbbaa.
%  \__, |_| \__,_\__,_|_\___|_||_\__/__/      ,aabbccddd'          'dddccbbaa.
%  |___/                                     ,abbccdd'         z      'ddccbba.
%                                           ,abccdd'      _ _ z         'ddccba.
%                                          ,abcdd'     '( -l- )'       |  'ddcba.
%                                          abccd'     /\ \ O / ___ ----|   'dccba
%                       
%         
% Nathan Cross feat. Casey Paquola
%   2020
%
%% 1. Setup
clc
%   i. Set key variables
keyword1 = 'task-All'; %<-------------------------------- Set to specific task - this is case sensitive (for all tasks: keyword1 = 'task-All')
keyword2 = '.gii';  %<-------------------------------- Set to format in which you would like to import (e.g. volumetric = '.nii'; surface = '.gii')
times = [{'ses-01','','Control'} %<------------------------¬
         {'ses-02','run-01','PreNap'} %<-------------------¬
         {'ses-02','run-02'}, 'PostNap']; % <-------------- Set up study specific session & run combinations (LEAVE AS IS FOR SDEP)
load('../SubjectName.mat');
tasks = {'ANT', 'Nback', 'PVT'}; % <--- list + order of tasks to concatenate - Needed ONLY IF keyword1 = 'task-All'

%   ii. Specify Schaefer Parcellation and Yeo Network templates 
fs = 'fsaverage5'; %<-------------------------------------- (Choose number of vertices as desired)
nparc = '400'; %<------------------------------------------ (Choose number of parcellations as desired)
nnetw = '17'; %<------------------------------------------- (Choose number of networks as desired)

%  iii. Load in the Schaeffer parcellations
load(sprintf('../labels/%s/ShfParcels/ShfLabels%s_%s.mat',fs,nparc,nnetw)) 
NumberOfParcels = max(rh_labels_Shf); % Number of Parcels used (saved for future)

%  iv. Load in the Yeo network templates
Yeo_7Clusters_ref = load('../labels/fsaverage5/YeoNetworks/1000subjects_clusters007_ref.mat'); % 7 Networks
Yeo_17Clusters_ref = load('../labels/fsaverage5/YeoNetworks/1000subjects_clusters017_ref.mat'); % 17 Networks

%	v. Map Schaeffer parcellations to the Yeo Networks 
for i = 1:length(lh_labels_Shf)
    lh_labels_Shf_Yeo(i,1) = lh_labels_Shf(i);
    lh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.lh_labels(i);
end
for i = 1:length(rh_labels_Shf)
    rh_labels_Shf_Yeo(i,1) = rh_labels_Shf(i);
    rh_labels_Shf_Yeo(i,2) = Yeo_17Clusters_ref.rh_labels(i);
end

% vi. add dependencies
addpath('../dependencies/customcolormap');
addpath('../dependencies/fireice');
addpath('../dependencies/cbrewer');
addpath('../dependencies/DNorm2');

%% 2. Create group-level state gradients
load(sprintf('../data/ConMat%s_%s.mat',nparc, keyword1));
embedding = group_embedding(ConMat, times);
%save('embedding.mat','embedding')

%% 3. Create subject-level gradients
indi_embed = individual_embeddings(ConMat, embedding);
%save('gradients/indi_embed.mat','indi_embed')

%% 3. Load embeddings and correlation matrices
load('gradients/embedding.mat');
load('gradients/indi_embed.mat');
load('times/Tasktimes.mat');
load('labels/Yeo17_Shf400.mat');
%load('data/Schaeffer400_task-All.mat');

