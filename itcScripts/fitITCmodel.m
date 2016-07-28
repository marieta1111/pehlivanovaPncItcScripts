% adapting script from summer 2014 to run on terminal
cd('/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/itcScripts')

load('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/ITCdata.mat')
%load('/Users/marieta/Dropbox/PNC_project/originalData/matlabData/RISKdata.mat')

load('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/itemOrderITC.mat')
load('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/itemOrderRisk.mat')

%% ITC

k_new=nan(size(itcC,2),1);
rsq_new=nan(size(itcC,2),1);
corr1_new=nan(size(itcC,2),1);
medRT_new=nan(size(itcC,2),1);
pctPred_new=nan(size(itcC,2),1);
noise_new=nan(size(itcC,2),1);

tic
for i=1:size(itcC,2)
i
    if ~isnan(itcC(1,i))

    %[hyperbolic] = ITCanalysis(itcC(:,i),itemOrderITC(:,2),itemOrderITC(:,3),...
     %       itemOrderITC(:,4),itemOrderITC(:,5),itcRT(:,i),ids(i))

     [hyperbolic] = KableLab_ITC_hyp(itcC(:,i),itemOrderITC(:,2),...
            itemOrderITC(:,4),itemOrderITC(:,5),itcRT(:,i))

            
    k_new(i)=hyperbolic.k;
    rsq_new(i)=hyperbolic.r2;
    corr1_new(i)=hyperbolic.RTandSubjValueCorr;
    medRT_new(i)=hyperbolic.medianRT;
    pctPred_new(i)=hyperbolic.percentPredicted;
    noise_new(i)=hyperbolic.noise;
    
    end
    
end
toc

save('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/newScript.mat', ...
    'k_new','rsq_new','noise_new','medRT_new','pctPred_new','corr1_new')

plot(k,k_new,'b*','MarkerFaceColor','b')
xlabel('k from script 2014')
ylabel('k from script 2016')

%% create a list with all immediate or all delayed choices, 5/4/15
allGood=find(strcmp(error,'NA'))
allImm=find(strcmp(error,'allimmediate'));
allDelayed=find(strcmp(error,'alldelayed'));

% comparing medRT 
[p h]=ranksum(medRT(allGood),medRT(allImm))
[p h]=ranksum(medRT(allGood),medRT(allDelayed))
[p h]=ranksum(medRT(allDelayed),medRT(allImm))

% saving data 
% numbers didn't get exported properly; manually corrected

csvwrite('/Users/marieta/Dropbox/PNC_project/csvData/DD_all_delayed_bblids_20150504.csv', sort(bblid(allDelayed)))
csvwrite('/Users/marieta/Dropbox/PNC_project/csvData/DD_all_immediate_bblids_20150504.csv', sort(bblid(allImm)))

%% risk

a=nan(size(riskC,2),1);
rsqR=nan(size(riskC,2),1);
errorR=cell(size(riskC,2),1);

for j= [1:size(riskC,2)]
    
j     
 if ~isnan(riskC(1,j))
    [output] = riskutilityfit(itemOrderRisk(:,2),itemOrderRisk(:,3),riskC(:,j),riskRT(:,j),ids(j))

    a(j)=output.a;
    rsqR(j)=output.r2;
    errorR{j}=output.errorcode;
 end

end

% using new script to get updated R-squared values
for j= [1:size(riskC,2)]
    
j     
 if ~isnan(riskC(1,j))
    [output] = agnosticrisk(itemOrderRisk(:,2),itemOrderRisk(:,3),riskC(:,j),riskRT(:,j))

    a(j)=output.a;
    rsqR(j)=output.r2;
    errorR{j}=output.errorcode;
 end

end
