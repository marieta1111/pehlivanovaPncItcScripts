% read the data

mask_nii = load_untouch_nii('/cbica/home/pehlivam/CT/n427_ctImages_2mm/n1359_antsCT_mask.nii.gz');
mask = mask_nii.img ;

data_lst = '/cbica/home/pehlivam/CT/n427_ctImages_2mm/n427_ants2mmCtPaths.csv';
fid=fopen(data_lst,'r');
datafullpath = textscan(fid,'%s %d\n');
fclose(fid);

data.y = datafullpath{1,2} ;
datafullpath = datafullpath{1,1} ;

count = numel(datafullpath);
info = load_untouch_header_only(datafullpath{1});

data.X = zeros(sum(mask(:)>0),count);
for i=1:count
    disp(i/count)
    nii = load_untouch_nii(datafullpath{i});    
    data.X(:,i) = nii.img(mask>0) ;
end

% load results and calculate reconstruction error
resultsPath='/cbica/home/pehlivam/CT/NMF/sge_job_output/'
numBases=[2:2:30];

RecError=zeros(length(numBases),1);
for b=1:length(numBases)
    disp(b/length(numBases))
    load([resultsPath 'NumBases' num2str(numBases(b)) '/OPNMF/ResultsExtractBases.mat'])  
    Est = B(mask>0,:)*C ;
    RecError(b) = norm(data.X-Est,'fro') ;    
    clear B C
end

% make figure
% 1) reconstruction error
figure;plot(numBases,RecError,'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Reconstruction error','fontsize',12)
set(gca,'fontsize',12)

% 2) gradient of reconstruction error
gr=diff(RecError)
figure;plot(numBases(2:end),gr,'r','LineWidth',3)
xlabel('Number of networks','fontsize',12)
ylabel('Gradient of reconstruction error','fontsize',12)
set(gca,'fontsize',12)
xlim([3 31])
set(gca,'XTick',[4:2:30])
hold on
%plot(20,gr(9),'ob','MarkerSize',14,'LineWidth',3.5);
%line([20 20],[gr(1) gr(9)-.50],'LineWidth',3)
plot(20,gr(9),'ob','MarkerSize',14,'LineWidth',3.5)

,'MarkerFaceColor','w')


