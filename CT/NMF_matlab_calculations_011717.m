load('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/demoBehavData/n427forSVRmatlab20170117.mat')

% format data from table to matrix: 1: remove bblid; 2: convert to array;

C1 = table2array(n427forSVRmatlab20170117(:,[2:27,29:35,37:end]));

%BMI = table2array(n427forSVRmatlab20170117(:,36));

C1(:,41)=(C1(:,2)-mean(C1(:,2))).^2; % add squared age term
 
% variables in C: 1: medu; 2: age; 3: sex; 4: T1 Rating; 5: k; 6: logk; 7:
% overall accuracy; 8: F1_Exec_Comp_Res_Accuracy; 9: F2_Social Cog; 10: F3 Memory accuracy; 
% 11:29: NMF Comps; 30: BAS Dr; 31: BAS Fs; 32: BAS RR; 33: BIS; 34: total
% substances tried; 35: tobacco; 36: alcohol; 37: marijuana; 38: average
% sub frequency; 39: ENV score; 40: race
% is Comp17 to be excluded

% pearson's correlations
[r p]=corr(C1(:,6),C1(:,11),'Type','Pearson')


% partial correlations with components, partialling out age, ageSq, and sex
whichCors=[11:29];

corCoefs=NaN; 
% col1 = partial corr; col2=partial corr p; col3 = raw correlation

for  k=1:length(whichCors)
    k
    x = horzcat(C1(:,6),C1(:,whichCors(k)),C1(:,[2 41 3]));
    
    [rho p] = partialcorr(x);
    corCoefs(k,1)=rho(1,2);
    corCoefs(k,2)=p(1,2);
    corCoefs(k,3)=corr(x(:,1),x(:,2));
end

% partial correlations with cognitive variables
[rho p] = partialcorr(C1(:,[6 7 2 41 3])) % Overall accuracy
[rho p] = partialcorr(C1(:,[6 8 2 41 3])) % ECRA
[rho p] = partialcorr(C1(:,[6 9 2 41 3])) % social cognition
[rho p] = partialcorr(C1(:,[6 10 2 41 3])) % Memory accuracy

% partial correlations with other variables
[rho p] = partialcorr(C1(:,[6 1 2 41 3]),'Type','Pearson','rows','complete') % medu
[rho p] = partialcorr(C1(:,[6 4 2 41 3]),'Type','Spearman','rows','complete') % T1 rating

[rho p] = partialcorr(C1(:,[6 35 2 41 3]),'Type','Spearman') % tobacco
[rho p] = partialcorr(C1(:,[6 36 2 41 3]),'Type','Spearman') % alcohol
[rho p] = partialcorr(C1(:,[6 37 2 41 3]),'Type','Spearman') % marijuana
[rho p] = partialcorr(C1(:,[6 34 2 41 3]),'Type','Spearman') % num substances tried

[rho p] = partialcorr(C1(:,[6 30 2 41 3]),'Type','Spearman') % BAS drive
[rho p] = partialcorr(C1(:,[6 31 2 41 3]),'Type','Spearman') % BAS fun seeking
[rho p] = partialcorr(C1(:,[6 32 2 41 3]),'Type','Spearman') % BAS reward responsiveness
[rho p] = partialcorr(C1(:,[6 33 2 41 3]),'Type','Spearman') % BIS 

[rho p] = partialcorr(C1(:,[6 1 2 41 3]),'Type','Pearson','Rows','Complete') % mother's education
