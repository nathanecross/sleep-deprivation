
%% LOAD IN STAGES AND PLOT SEGMENTS
home = pwd;
cd('../../')
rootpath = pwd;
cd(home);
filepath = [rootpath '/data/sleep/Sleep_stage_times.xlsx'];
Stages = readtable(filepath);
Stages = Stages(table2array(Stages(:,6))>359,:);
Stage2 = Stages(strcmp(table2cell(Stages(:,3)),"'NREM2'"),:);
Stage3 = Stages(strcmp(table2cell(Stages(:,3)),"'NREM3'"),:);
load('../../data/subject.mat'); SubjectName = subject.SubjectList;

%% PLOT SEGMENTS OF N2
for i = 1:height(Stage2)
    if i == 10
            plot([table2array(Stage2(i,4))+2130 table2array(Stage2(i,5))+2130],[table2array(Stage2(i,2)) table2array(Stage2(i,2))]); hold on;
            ylim([0 21])
    elseif i == 13
            plot([table2array(Stage2(i,4))+1670 table2array(Stage2(i,5))+1670],[table2array(Stage2(i,2)) table2array(Stage2(i,2))]); hold on;
            ylim([0 21])
    elseif i == 18
            plot([table2array(Stage2(i,4))+1770 table2array(Stage2(i,5))+1770],[table2array(Stage2(i,2)) (table2array(Stage2(i,2)))]); hold on;
            ylim([0 21])
    else
        plot([table2array(Stage2(i,4)) table2array(Stage2(i,5))],[table2array(Stage2(i,2)) (table2array(Stage2(i,2)))]); hold on;
        ylim([0 21])
    end
    yticks([0:21])
    yticklabels([{''}; SubjectName(:); {''}])
end
% PLOT SEGMENTS OF N3 
figure;
for i = 1:height(Stage3)
    if i == 10
            plot([table2array(Stage3(i,4))+2130 table2array(Stage3(i,5))+2130],[table2array(Stage3(i,2)) table2array(Stage3(i,2))]); hold on;
            ylim([0 21])
    elseif i == 13
            plot([table2array(Stage3(i,4))+1670 table2array(Stage3(i,5))+1670],[table2array(Stage3(i,2)) table2array(Stage3(i,2))]); hold on;
            ylim([0 21])
    elseif i == 18
            plot([table2array(Stage3(i,4))+1770 table2array(Stage3(i,5))+1770],[table2array(Stage3(i,2)) table2array(Stage3(i,2))]); hold on;
            ylim([0 21])
    else
        plot([table2array(Stage3(i,4)) table2array(Stage3(i,5))],[table2array(Stage3(i,2)) table2array(Stage3(i,2))]); hold on;
        ylim([0 21])
    end
    yticks([0:21])
    yticklabels([{''}; SubjectName(:); {''}])
end

%% Setup 
surf_dir = '../../data/Parcels_cat/';
vol_dir = '../../data/Subjectdata_nii_cat/';
keyword = 'nap';
load('../../data/subject.mat');
SubjectName = subject.SubjectList;
NumberOfSubjects = length(SubjectName);
Stages = readtable('../../data/sleep/Sleep_times_onechunk.xlsx');
%% Extract 26 min segments of NREM sleep

for z = 1:NumberOfSubjects
    sprintf('%s',char(SubjectName{z}))
    skip=false;
    % Surface data
    try
        load([surf_dir char(SubjectName{z}) '/nap/' ...
                  char(SubjectName{z}) '_' keyword '.mat']);
    catch
        sprintf('No data for %s, skipping.',char(SubjectName{z}));
        skip=true;
    end
    
    if isnan(Stages.Start(z))==1
       skip=true; 
    end
    
    if skip==false
        % Truncate data
        Parcels = Parcels(:,(Stages.Start(z)):(Stages.End(z)));

        % Make output directory
        if ~exist(sprintf('%s/',[surf_dir char(SubjectName{z}) '/nap_trunc/']), 'dir')
            mkdir(sprintf('%s/',[surf_dir char(SubjectName{z}) '/nap_trunc/']))
        end

        % Save truncated data
        save([surf_dir char(SubjectName{z}) '/nap_trunc/' ...
                  char(SubjectName{z}) '_' keyword '.mat'],'Parcels');
    else
        continue
    end
    % Volume data
    try
        load([vol_dir char(SubjectName{z}) '/nap/' ...
                  char(SubjectName{z}) '_' keyword '.mat']);
        vols = fieldnames(subjectdata);
    catch
        sprintf('No data for %s, skipping.',char(SubjectName{z}));
        skip=true;
    end     
    
        % Truncate data
        for v = 1:length(vols)
            if contains(vols(v), 'SubjectName') == 0
                subjectdata.(char(vols(v))) = subjectdata.(char(vols(v)))(:,(Stages.Start(z)):(Stages.End(z)));
            end
        end
        
    if skip==false
        % Make output directory
        if ~exist(sprintf('%s/',[vol_dir char(SubjectName{z}) '/nap_trunc/']), 'dir')
            mkdir(sprintf('%s/',[vol_dir char(SubjectName{z}) '/nap_trunc/']))
        end

        % Save truncated data
        save([vol_dir char(SubjectName{z}) '/nap_trunc/' ...
                  char(SubjectName{z}) '_' keyword '.mat'],'subjectdata');
    else
        continue
    end
                
end

%% Extract 7 min segments of N2 sleep
surf_dir = '../../data/Parcels100/';
vol_dir = '../../data/Subjectdata_nii_cat/';
in_keyword = 'nap';
out_keyword = 'nap_N2';
Subs = unique(Stage2.Var2);
for i = 1:length(Subs)
    
    sprintf('Running Subject %s',char(SubjectName{Subs(i)}))
    
    % find stage 2 segments          
    S2 = Stage2(Stage2.Var2==Subs(i),:); 
    
    %load data
    load([surf_dir char(SubjectName{Subs(i)}) '/nap/' ...
                  char(S2.Var1(1)) '_' in_keyword '.mat']);
    
    
    
    % take first segment >120 TRs
    for row = 1: height(S2)
       dur =  table2array(S2(row,8)) - table2array(S2(row,7)) + 1;
       if dur > 120
           break
       end
    end
    
    start = table2array(S2(row,7));
    finish = table2array(S2(row,7))+120;
    % Truncate data
    Parcels = Parcels(:,start:finish);
    size(Parcels,2)
    % Make output directory
    if ~exist(sprintf('%s/',[surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/']), 'dir')
        mkdir(sprintf('%s/',[surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/']))
    end

    % Save truncated data
    save([surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/' ...
              char(S2.Var1(1)) '_' out_keyword '.mat'],'Parcels');
    
end

% Extract 7 min segments of N3 sleep
Subs = unique(Stage3.Var2);
out_keyword = 'nap_N3';
for i = 1:length(Subs)
    
    sprintf('Running Subject %s',char(SubjectName{Subs(i)}))
    
    % find stage 2 segments          
    S3 = Stage3(Stage3.Var2==Subs(i),:); 
    
    %load data
    load([surf_dir char(SubjectName{Subs(i)}) '/nap/' ...
                  char(S3.Var1(1)) '_' in_keyword '.mat']);
    
    
    
    % take first segment >120 TRs
    for row = 1: height(S3)
       dur =  table2array(S3(row,8)) - table2array(S3(row,7)) + 1;
       if dur > 120
           break
       end
    end
    
    start = table2array(S3(row,7));
    finish = table2array(S3(row,7))+120;
    
    % Truncate data
    Parcels = Parcels(:,start:finish);
    
    % Make output directory
    if ~exist(sprintf('%s/',[surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/']), 'dir')
        mkdir(sprintf('%s/',[surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/']))
    end

    % Save truncated data
    save([surf_dir char(SubjectName{Subs(i)}) '/nap_trunc/' ...
              char(S3.Var1(1)) '_' out_keyword '.mat'],'Parcels');
    
end
