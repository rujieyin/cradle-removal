function MCAdecomposition_cluster(row_ind, result_path, job_id)

addpath(genpath('~/matlab/MCALabWithUtilities'));
addpath(genpath('~/matlab/cradleremoval'));
addpath('~/matlab/CurveLab-2.1.3/fdct_wrapping_matlab/');

row_ind = str2num(row_ind);
load('/home/grad/rachel/Documents/MATLAB/MCAdecomp_Ghent_8/Ghent_8.mat','subimg');%%
subpatch = subimg(row_ind,:);
texture = cell(length(subpatch),1);
for col_ind = 1 : length(subpatch)
    img = subpatch{col_ind};
    %% check img, display = 0; non-display of results
    MCAdecomp;
    texture{col_ind} = parts(:,:,2);
end
path = [result_path 'texture_row' job_id '/'];
if ~exist(path,'dir')
    mkdir(path);
end
% if ~exist([path 'texture.mat'], 'file')
    save([path 'texture.mat'],'texture');
% end
disp('result saved');