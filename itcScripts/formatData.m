% change path to reflect new name of folder
load('/Users/marieta/Dropbox/PNC_project/originalData/matlabData/origItcRiskData.mat')

ids=str2num(cell2mat(textdata(:,1)));
%bblid=str2num(cell2mat(textdata(:,2)));

%% ITC

choicesITC=textdata(:,52:2:118);
rtsITC=textdata(:,53:2:119);

itcC=nan(size(choicesITC,2),size(data,1));
itcRT=nan(size(choicesITC,2),size(data,1));


for i=1:size(data,1)
    for j=1:size(choicesITC,2)
       
        if ~isempty(str2num(cell2mat(choicesITC(i,j))))
            itcC(j,i)=str2num(cell2mat(choicesITC(i,j)));
        end
        if ~isempty(str2num(cell2mat(rtsITC(i,j))))
            itcRT(j,i)=str2num(cell2mat(rtsITC(i,j)));
        end
    end    
end

[row col]=find(isnan(itcC));
itcC(:,unique(col))=[];
itcRT(:,unique(col))=[];
itcC=itcC-1;

idsITC=ids;
idsITC(unique(col))=[];

%% risk

riskC=data(:,1:2:81)'; % choice 
riskRT=data(:,2:2:82)'; % reaction time

%{
[row col]=find(isnan(riskC));
riskC(:,unique(col))=[];
riskC=2-riskC;
riskRT(:,unique(col))=[];

idsR=ids;
idsR(unique(col))=[];
%}

%% BIS BAS
BASdrive=zeros(size(textdata,1),1);
BASfs=zeros(size(textdata,1),1);
BASrr=zeros(size(textdata,1),1);
BIS=zeros(size(textdata,1),1);

for i=1:size(textdata,1)

% BAS drive
qs=[6 11 13 21];
 if ~strcmp(textdata(i,qs(1)),'')
     for j=1:length(qs)
            BASdrive(i)=BASdrive(i)+str2num(cell2mat(textdata(i,qs(j)))); 
            %BASdrive(i)
     end
 else 
     BASdrive(i)=nan;
 end
 
% BAS fun seeking
qs=[8 12 16 20];
 if ~strcmp(textdata(i,qs(1)),'')
     for j=1:length(qs)
            BASfs(i)=BASfs(i)+str2num(cell2mat(textdata(i,qs(j)))); 
            %BASdrive(i)
     end
 else 
     BASfs(i)=nan;
 end

% BAS reward responsiveness
qs=[7 9 15 18 23];
 if ~strcmp(textdata(i,qs(1)),'')
     for j=1:length(qs)
            BASrr(i)=BASrr(i)+str2num(cell2mat(textdata(i,qs(j)))); 
            %BASdrive(i)
     end
 else 
     BASrr(i)=nan;
 end

% BIS
qs=[5 10 14 17 19 22 24];
 if ~strcmp(textdata(i,qs(1)),'')
     
     BIS(i)=3-str2num(cell2mat(textdata(i,qs(1)))); 
     for j=2:length(qs)
         BIS(i)=BIS(i)+str2num(cell2mat(textdata(i,qs(j)))); 
     end
 else 
     BIS(i)=nan;
 end

i 
j
end

%% BSSS
BSSS=nan(size(textdata,1),1);
for i=1:size(textdata,1)
    
    if ~strcmp(textdata(i,26),'')
    BSSS(i)= str2num(cell2mat(textdata(i,26))) + str2num(cell2mat(textdata(i,27))) + ...
         str2num(cell2mat(textdata(i,28))) + str2num(cell2mat(textdata(i,29)));   
    end

end





























