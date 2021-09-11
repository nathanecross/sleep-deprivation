function [Fmat, Pmat, Diff, Tmat, PostPmat] = model_rm(sub_dir, subjects, contrasts, covariates)
 
    NumberOfSubjects = length(subjects);    

    for c = 1:length(contrasts)
        for z = 1:NumberOfSubjects
            load( [sub_dir char(subjects{z}) '/' char(contrasts(c)) '/' char(subjects{z}) '_task-All.mat']) 
            holdmat(c,z,:,:) = ConMat.corr_z;            
        end
    end
    NumberOfParcels = length(holdmat(1,1,:,1));
    % Pre-allocate matrices  
    Fmat=NaN(NumberOfParcels);
    Pmat=NaN(NumberOfParcels);
    for i=1:NumberOfParcels
        for j = 1:NumberOfParcels
            %sprintf('%s,%s',char(string(i)),char(string(j)))
            if sum(isnan(holdmat(1,:,i,j))) > 1
                continue
            else
                a = array2table(squeeze(holdmat(:,:,i,j))');
                if nargin > 2 % check for covariates
                    for c = 1:length(covariates(:,1))
                        a.(char(covariates(c))) = cell2mat(covariates(c,2));
                        if c == 1
                            covs = (char(covariates(c)));
                        else
                            covs = [covs '+' (char(covariates(c)))];
                        end
                    end
                    a.Properties.VariableNames = [contrasts covariates(:,1)']; 
                    Time = (1:length(contrasts))';                
                    rm = fitrm(a, sprintf('%s-%s ~ %s',contrasts{1},contrasts{end}, covs), 'WithinDesign', Time); 
                    ranovatbl = ranova(rm);
                    Fmat(i,j) = table2array(ranovatbl('(Intercept):Time','F'));
                    Pmat(i,j) = table2array(ranovatbl('(Intercept):Time','pValue'));
                else
                    a.Properties.VariableNames = contrasts; 
                    Time = (1:length(contrasts))';
                    rm = fitrm(a, sprintf('%s-%s ~ 1',contrasts{1},contrasts{end}), 'WithinDesign', Time); 
                    ranovatbl = ranova(rm);
                    Fmat(i,j) = table2array(ranovatbl('(Intercept):Time','F'));
                    Pmat(i,j) = table2array(ranovatbl('(Intercept):Time','pValue'));
                end
            end
        end
    end

    
    
    % POST HOC t-TESTS %
    
    % Create contrast combinations
    for c = 2:length(contrasts)
        con{c-1,1} = c-1;
        con{c-1,2} = c;
        con{c-1,3} = [char(contrasts{c-1}) '_' char(contrasts{c})];
    end
    % one last one
    con{length(contrasts),1} = 1;
    con{length(contrasts),2} = length(contrasts);
    con{length(contrasts),3} = [char(contrasts{1}) '_' char(contrasts{end})];
    
    % Run through contrasts
    for k = 1:length(con(:,1))
        % Load data
        for z = 1:NumberOfSubjects
            load( [sub_dir char(subjects{z}) '/' char(contrasts(cell2mat(con(k,1)))) '/' char(subjects{z}) '_task-All.mat']) 
            holdmat1(z,:,:) = ConMat.corr_z;
            load( [sub_dir char(subjects{z}) '/' char(contrasts(cell2mat(con(k,2)))) '/' char(subjects{z}) '_task-All.mat']) 
            holdmat2(z,:,:) = ConMat.corr_z;            
        end
        
        NumberOfParcels = length(holdmat1(1,:,1));
        
        % Pre-allocate matrices
        Diff=squeeze(NaN((length(con(:,1))),NumberOfParcels,NumberOfParcels));
        Tmat=squeeze(NaN((length(con(:,1))),NumberOfParcels,NumberOfParcels));
        PostPmat=squeeze(NaN((length(con(:,1))),NumberOfParcels,NumberOfParcels));
        
        % Compare contrasts
        for i=1:NumberOfParcels
            for j = 1:NumberOfParcels
                Diff(k,i,j) = mean(squeeze(holdmat2(:,i,j))-squeeze(holdmat1(:,i,j)));
                if sum(isnan(holdmat1(:,i,j))) > 1
                    continue
                else
                    [~,PostPmat(k,i,j),~,stats] = ttest(squeeze(holdmat2(:,i,j)), squeeze(holdmat1(:,i,j)));
                    Tmat(k,i,j) = stats.tstat;
                end
            end 
        end
    end
end