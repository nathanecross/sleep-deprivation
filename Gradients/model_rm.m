function [Fmat, Pmat] = model_rm(ConMat,nparc)
    subs = struct2cell(ConMat.times);  
    SubjectName = fieldnames(subs{1});
    subs = struct2cell(subs{1});  
    NumberOfSubjects = length(subs);
    NumberOfParcels = str2num(nparc);
    times = fieldnames(ConMat.times);
    holdmat=zeros(NumberOfSubjects, 3);
    Fmat=NaN(NumberOfParcels);
    Pmat=NaN(NumberOfParcels);
    for i=1:NumberOfParcels
        for j = 1:NumberOfParcels
            for z=1:NumberOfSubjects
                if isfield(ConMat.times.(char(times(1))),(SubjectName{z})) == 1 && isfield(ConMat.times.(char(times(2))),(SubjectName{z})) == 1 && isfield(ConMat.times.(char(times(3))),(SubjectName{z})) == 1 
                    holdmat(z,1) = ConMat.times.(char(times(1))).(SubjectName{z}).ShfCorr_z(i,j); 
                    holdmat(z,2) = ConMat.times.(char(times(2))).(SubjectName{z}).ShfCorr_z(i,j);
                    holdmat(z,3) = ConMat.times.(char(times(3))).(SubjectName{z}).ShfCorr_z(i,j);
                else
                    continue
                end
            end
            if isnan(ConMat.times.(char(times(1))).(SubjectName{z}).ShfCorr_z(i,j)) == 1 || isnan(ConMat.times.(char(times(2))).(SubjectName{z}).ShfCorr_z(i,j)) == 1
                continue
            else
            a = array2table(holdmat);
    %       a.sex = [1;2;1;1;2;2;2;1;2;2;2;2;2;2;1;1;2;2;1];  % <--- add sex + age as covariates. 1=M, 2=F
    %       a.age = [20;22;19;21;19;21;22;21;20;30;22;18;18;21;22;22;19;24;21];
            a.Properties.VariableNames = {'Session1','Session2','Session3'}; 
            Time = [1;2;3];
            rm = fitrm(a, 'Session1-Session3 ~ 1', 'WithinDesign', Time); 
            ranovatbl = ranova(rm);
            Fmat(i,j) = table2array(ranovatbl('(Intercept):Time','F'));
            Pmat(i,j) = table2array(ranovatbl('(Intercept):Time','pValue'));
            end
        end 
    end

end
