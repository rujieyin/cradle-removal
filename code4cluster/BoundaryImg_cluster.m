function BoundaryImg_cluster(result_path,job_id)

addpath('/home/grad/rachel/matlab/cradleremoval/texture_removal/')
addpath(genpath('/home/grad/rachel/matlab/FFST'));
load('/home/grad/rachel/matlab/cradleremoval/block19a.mat','subpatch_new');
ind = str2num(job_id);
subpatch = subpatch_new((3*ind-2):3*ind);
load(['/home/grad/rachel/Documents/MATLAB/Preparation/result/job' job_id '/output.mat'],'opt');
bdImg = cell(3,1);
lowfreqImg = cell(3,1);

if mod(ind,1) == 0
    opt{1}.BDdemean = 1;
    opt{2}.BDdemean = 1;
    opt{3}.BDdemean = 1;
end

for i = 1:3
    [bdImg{i},lowfreqImg{i}] = BoundaryImg(subpatch{i},opt{i});
end
path = [result_path 'result/job' job_id '/'];
if ~exist(path,'dir')
    mkdir(path);
end
save([path 'bdImg.mat'],'bdImg','lowfreqImg');
