function MCAdecomposition_cluster( ngrid, result_path, load_mat, job_id )

addpath(genpath('~/matlab/MCALabWithUtilities'));
addpath(genpath('~/matlab/cradleremoval'));
addpath('~/matlab/CurveLab-2.1.3/fdct_wrapping_matlab/');

% row_ind = str2num(row_ind);
load(load_mat,'subimg');%%
npatch = length(subimg(:))
ngrid = str2num(ngrid);
job_id = str2num(job_id);
jobsize = ceil(npatch/ngrid)
jobsize*(job_id-1)+1
jobind = jobsize*(job_id-1)+1 : min(jobsize*job_id, npatch) 
subpatch = subimg(jobind);
% subpatch = subimg(row_ind,:);
texture = cell(length(subpatch),1);
for col_ind = 1 : length(subpatch)
    img = subpatch{col_ind};
    %% check img, display = 0; non-display of results
    MCAdecomp;
    texture{col_ind} = parts(:,:,2);
end
path = [result_path 'texture_job_id' num2str(job_id) '/'];
if ~exist(path,'dir')
    mkdir(path);
end
% if ~exist([path 'texture.mat'], 'file')
    save([path 'texture.mat'],'texture');
% end
disp('result saved');