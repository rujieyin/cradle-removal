function blocktexturefeature(result_path,job_id)

ind = str2num(job_id);
path = [result_path 'result/job' job_id];
if ~exist(path,'dir')
    mkdir(path);
end
cd /home/grad/rachel/matlab/cradleremoval/
load('block19a.mat','texture');
load('block19a.mat','mask');

output = struct();
output.opt = cell(3,1);
output.featureST = cell(3,1);
output.Y = cell(3,1);

texture = texture((3*(ind-1)+1):3*ind);
mask = mask((3*(ind-1)+1):3*ind);
for i = 1:length(texture)
opt = struct();
opt.direction = 'horizontal';
    img = texture{i};
% img = texture19a{i,j};
opt.imgsize = size(img);
% opt.transform = 'curvelet';
% opt.curveletisreal = 1;
% cw = fdct_wrapping(img,opt.curveletisreal,1,7);
% mask = ones(opt.imgsize);
% [cwnew,imgnew,opt] = angleselection2(cw,mask,opt);
% opt.mask_h = blockmask{i,j}.h;
% opt.mask_v = blockmask{i,j}.v;
% [featureCW,opt] = getCWfeature(cwnew,opt);
% % [region,distscore,trans] = featuresimilarityregion(feature,opt,img);

opt.transform = 'shearlet';
opt.mask_v = mask{i}.mask_v;
opt.mask_h = mask{i}.mask_h;
addpath(genpath('/home/grad/rachel/matlab/FFST'));
addpath('/home/grad/rachel/matlab/cradleremoval/texture_removal/');
if max(size(img)) <= 512
    ST = shearletTransformSpect(img);
    ST = real(ST);
end
[j,k,cone] = arrayfun(@(x)shearletScaleShear(x), 1:size(ST,3));
[STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
% imgnew = inverseShearletTransformSpect(STnew);
[featureST,opt] = getSTfeature(STnew,j,opt);

    Y = [featureST.cradle;featureST.free];
originImg = invfeatureST(Y,opt,1);
opt.originImg = originImg;

output.opt{i} = opt;
output.featureST{i} = featureST;
output.Y{i} = Y;
end
if ~exist(path,'dir')
    mkdir(path);
end
if ~exist([path '/output.mat'], 'file')
    save([path '/output.mat'],'-struct','output');
end