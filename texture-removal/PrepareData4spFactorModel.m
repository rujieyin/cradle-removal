function [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt)

% this function split the whole dataset Y into labeled and unlabled
% datasets as an input for the unsuperised model

% -- get original mask -- %
switch opt.direction
    case 'vertical'
        mask0 = opt.mask_v;
        downsamplesize = opt.vdownsamplesize;
        if isfield(opt,'bdmask_v')
            bdmask = opt.bdmask_v > .3;
        else
            bdmask = zeros(size(mask0));
        end
    case 'horizontal'
        mask0 = opt.mask_h;
        downsamplesize = opt.hdownsamplesize;
        if isfield(opt,'bdmask_h')
            bdmask = opt.bdmask_h > .3;
        else
            bdmask = zeros(size(mask0));
        end
end
resizemask = @(x)logical(imresize(x,downsamplesize,'bilinear'));
mask0 = resizemask(mask0);
bdmask = resizemask(bdmask);

if isfield(opt,'method')
    method = opt.method;
else
    method = 'supervised';
end

switch method
    case 'supervised'
            mcradle = mask0 & (~bdmask);
            mnoncradle = (~mask0) & (~bdmask);
    case 'unsupervised'        
        se = strel('square',40);
        mcradle = imerode(mask0,se);
        mnoncradle = imerode(~mask0,se);
        mcradle = logical(mcradle);
        mnoncradle = logical(mnoncradle);
end

switch opt.direction
    case 'vertical'
        Indcradle = mcradle(opt.vfeatureind.cradle);
        Indnoncradle = mnoncradle(opt.vfeatureind.free);
        featureind = opt.vfeatureind;
    case 'horizontal'
        Indcradle = mcradle(opt.hfeatureind.cradle);
        Indnoncradle = mnoncradle(opt.hfeatureind.free);
        featureind = opt.hfeatureind;
end
featureind.cradle = featureind.cradle(:)';
featureind.free = featureind.free(:)';
    Z = Y([~Indcradle(:); ~Indnoncradle(:)],:);
    Y = Y([Indcradle(:); Indnoncradle(:)],:);
    mask = logical([ones(sum(Indcradle),1);zeros(sum(Indnoncradle),1)]);
    opt.factormodel_setup = struct();
    opt.factormodel_setup.IndY = [featureind.cradle(Indcradle),featureind.free(Indnoncradle)];
    opt.factormodel_setup.Ymask = mask;
    opt.factormodel_setup.Zmask = logical([ones(sum(~Indcradle),1);zeros(sum(~Indnoncradle),1)]);
    opt.factormodel_setup.IndZ = [featureind.cradle(~Indcradle),featureind.free(~Indnoncradle)];
    loc = struct();
    loc.Y = opt.factormodel_setup.IndY;
    loc.Z = opt.factormodel_setup.IndZ;
%     loc = [opt.factormodel_setup.IndY,opt.factormodel_setup.IndZ];

end


