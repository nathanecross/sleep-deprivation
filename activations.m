clc;

% Set variables 
keyword = 'Nback';
tasklength = 192;  % # of TRs in fMRI run

% Set your contrasts for conditions
C = [-1 1 0]; %e.g. this is for 
%% Obtain active areas in response to task onset
load(sprintf('data/Schaeffer400_task-%s.mat',keyword));
load(sprintf('times/%s/blocktimes.mat',keyword));
mat = zeros(400,tasklength);
SubjectName = fieldnames(Schaeffer.times.Control);
NumberOfSubjects = length(SubjectName);
times = fieldnames(Schaeffer.times);
Tval = zeros(length(times),NumberOfSubjects,400);
for t = 1:length(times)
    for z = 1:NumberOfSubjects
        
        % import and arrange task regressors
        if strcmp(keyword, 'ANT') == 1
           cov1 = Tasktimes.(char(times(t))).DoubleCuetimes_conv(z,:)'; cov1=cov1(3:tasklength,1);
           cov2 = Tasktimes.(char(times(t))).ValidCuetimes_conv(z,:)'; cov2=cov2(3:tasklength,1);
           cov3 = Tasktimes.(char(times(t))).NoCuetimes_conv(z,:)'; cov3=cov3(3:tasklength,1);
           cov4 = Tasktimes.(char(times(t))).InvalidCuetimes_conv(z,:)'; cov4=cov4(3:tasklength,1);
           X = [cov1 cov2 cov3 cov4 zeros(length(cov1),1)];
        elseif strcmp(keyword, 'Nback') == 1
           cov0 = Tasktimes.(char(sprintf("%s",char(times(t))))).ZeroBack_conv(z,:)';  cov0=[cov0(3:tasklength,1)];
           cov1 = Tasktimes.(char(sprintf("%s",char(times(t))))).OneBack_conv(z,:)';  cov1=[cov1(3:tasklength,1)];
           cov2 = Tasktimes.(char(sprintf("%s",char(times(t))))).TwoBack_conv(z,:)';  cov2=[cov2(3:tasklength,1)];
           X = [cov0 cov2 zeros(length(cov1),1)];
        elseif strcmp(keyword, 'PVT') == 1
           cov1 = Tasktimes.(char(sprintf("%s",char(times(t))))).Blocktimes_conv(z,:)';  cov1=[cov1(3:tasklength,1)];
           X = [cov1 zeros(length(cov1),1)];
        end

        % design matrix, df
        X = X;
        df = tasklength - rank(X);
        if t == 1 && z == 1
            Xs = spm_DesMtx('sca', X);
            figure
            colormap('gray')
            image((Xs+1)*32)
            title(sprintf('Design matrix for %s',keyword))
        end
                  
        % the analysis, giving slopes (=means) 
        for i = 1:400
            mat(i,:) = Schaeffer.times.(char(times(t))).(SubjectName{z}).ShfNetMean.(['p' char(string(i))]);
            Y = mat(i,3:tasklength)';
            B = pinv(X)*Y;

            % t statistic and significance test
            RSS   = sum((Y - X*B).^2);
            MRSS  = RSS / df;
            SE    = sqrt(MRSS*(C*pinv(X'*X)*C'));
            T(i)  = C*B./SE;
            P(i)  = 1-spm_Tcdf(T(i), df); % upper tail p
        end
        
        Tval(t,z,:) = T; %first level t-test
        Pval(t,z,:) = P;
    end

end


%% Second level analysis and significance testing
% Significant activatons at baseline
BL = squeeze(Tval(1,:,:)); 
[h,p,ci,stats] = ttest(BL); tBL = squeeze(stats.tstat); [sigBL,~,~,adjp] = fdr_bh(p);
table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
BLsig.(char(keyword)) = [num2cell(table.Var1) table.Var2 num2cell(tBL)' num2cell(adjp')];
BLsig.(char(keyword)) = BLsig.(char(keyword))(sigBL,:);

% Significant changes in activations across sessions+runs
PRE = squeeze(Tval(2,:,:));
POST = squeeze(Tval(3,:,:));
clear stats table
for p = 1:400 % for each node
    in_table = table([1:size(BL,1)]', BL(:,p), ...
        PRE(:,p), POST(:,p), ...
        'VariableNames', {'subject', 'state1', 'state2', 'state3'});

    rm = fitrm(in_table, 'state1-state3~1');
    ranovatbl = ranova(rm);
    stats(p,:) = table2array(ranovatbl(1,:));
end
rois = fdr_bh(stats(:,5)); % FDR corrected across all cortex at BL
table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
parcels_stats = [num2cell(table.Var1) table.Var2 num2cell(stats(:,4:5))];
uncorr_rois = parcels_stats(stats(:,5)<0.05,:); % Uncorrected changes in activations

% FDR corrected across activations at BL
nonsigBL = ~sigBL;
all = [1:400; stats(:,5)'];
keep = all(:,logical(sigBL));
holdd = all(:,nonsigBL); holdd(2,:) = 1;
[h, crit_p, adj_ci_cvrg, keep(2,:)]=fdr_bh(keep(2,:));
p = array2table([keep holdd]');
p = sortrows(p);
p = table2array(p);
vec = stats(:,4); 
table = readtable('labels/Schaefer2018_400Parcels_17Networks_order.txt');
Sig.(char(keyword)) = [num2cell(table.Var1) table.Var2 num2cell(vec) num2cell(p(:,2))];
Sig.(char(keyword)) = Sig.(char(keyword))(p(:,2)<0.05,:);


Activations.BL = BL; Activations.PRE = PRE; Activations.POST = POST;
save(sprintf('data/activations_%s.mat',keyword), 'Activations');

