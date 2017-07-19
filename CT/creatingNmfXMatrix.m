% read in nifti files used for NMF parcellation and create X matrix for
% visualization

cd('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm')

%% read in text file with 427 image names
filename = '/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/antsCT/ctImages_2mm/listOf427Ct2mmImages.txt';
delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
listOf427Ct2mmImages = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% mask 
mask = load_nifti('n1359_antsCT_mask.nii.gz');
mask_vector = reshape(mask.vol,[prod(size(mask.vol)) 1]);
keep_voxels = find(mask_vector==1);  % voxels in mask 

% sample voxels from mask
howMany = 10000;
blah=randsample(128155,howMany);

% preallocate X matrix
X=nan(howMany, length(listOf427Ct2mmImages));

% read in nifti files 
for j=1:size(listOf427Ct2mmImages,1)
    nii_file = load_nifti(listOf427Ct2mmImages{j});
    test = reshape(nii_file.vol,[prod(size(nii_file.vol)) 1]);
    test_masked = test(keep_voxels);
    sampled_test_masked = test_masked(blah);
    X(:,j) = sampled_test_masked; 
end

imagesc(X)
axis off

set(gca,'position',[0 0 3 3],'units','inches')

% formattting B matrix
load('/data/joy/BBL/projects/pehlivanovaPncItc/CT/NMF/NMF20results/ResultsExtractBases.mat')
B_masked = B(keep_voxels,:);
sampled_B_masked = B_masked(blah,:);

imagesc(sampled_B_masked)
set(gca,'position',[0 0 0.75 3],'units','inches')

axis off
colormap(flipud(bone))