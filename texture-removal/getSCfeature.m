function [feature,opt] = getSCfeature(S,opt)

% this function extract feature vectors from Scattering coefficients,
% coefficients in the same layer are combined, different layers are left
% separate.

if nargin < 2
    opt = struct();
end

if ~isfield(opt,'featuredim')% along which dim feature vectors are stored
    opt.featuredim = 'row';
end

% whether retain 0th order scattering coefficients
if isfield(opt,'basecoeff') & opt.basecoeff
    feature = cellfun(@(x)x.signal, S,'UniformOutput',0);
else
    feature = cellfun(@(x)x.signal,S(2:end),'UniformOutput',0);
end

opt.level = length(feature);
opt.featuresize = cell(opt.level,1);
for i = 1:opt.level
    opt.featuresize{i} = size(feature{i}{1});
    if strcmp(opt.featuredim,'row')
        feature{i} = cellfun(@(x)reshape(x,[],1),feature{i},'UniformOutput',0);
        feature{i} = cell2mat(feature{i});
    else
        feature{i} = cellfun(@(x)reshape(x,1,[]),feature{i},'UniformOutput',0);
        feature{i} = cell2mat(feature{i});
    end
end

if isfield(opt,'mask_h')
    opt.featuremask_h = arrayfun(@(i)imresize(opt.mask_h,'OutputSize',opt.featuresize{i},'Method','bilinear'),1:2,'UniformOutput',0);
    if strcmp(opt.featuredim,'row')
        opt.featuremask_h = cellfun(@(x)reshape(x,[],1),opt.featuremask_h,'UniformOutput',0);
    else
        opt.featuremask_h = cellfun(@(x)reshape(x,1,[]),opt.featuremask_h,'UniformOutput',0);
    end
end

if isfield(opt,'mask_v')
    opt.featuremask_v = arrayfun(@(i)imresize(opt.mask_v,'OutputSize',opt.featuresize{i},'Method','bilinear'),1:2,'UniformOutput',0);
    if strcmp(opt.featuredim,'row')
        opt.featuremask_v = cellfun(@(x)reshape(x,[],1),opt.featuremask_v,'UniformOutput',0);
    else
        opt.featuremask_v = cellfun(@(x)reshape(x,1,[]),opt.featuremask_v,'UniformOutput',0);
    end
end

