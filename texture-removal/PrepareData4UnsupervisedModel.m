function [Y,mask,Z,opt] = PrepareData4UnsupervisedModel(Y,opt)

% this function split the whole dataset Y into labeled and unlabled
% datasets as an input for the unsuperised model

% -- get original mask -- %
switch opt.direction
    case 'vertical'
mask0 = opt.mask_v;
mask0 = imresize(mask0,opt.vdownsamplesize,'bilinear');
    case 'horizontal'
        mask0 = opt.mask_h;
        mask0 = imresize(mask0,opt.hdownsamplesize,'bilinear');
end


se = strel('square',40);
mcradle = imerode(mask0,se);
mnoncradle = imerode(~mask0,se);
switch opt.direction
    case 'vertical'
Indcradle = mcradle(opt.vfeatureind.cradle);
Indnoncradle = mnoncradle(opt.vfeatureind.free);
    case 'horizontal'
        Indcradle = mcradle(opt.hfeatureind.cradle);
        Indnoncradle = mnoncradle(opt.hfeatureind.free);
end
Indcradle = logical(Indcradle);
Indnoncradle = logical(Indnoncradle);
Z = Y([~Indcradle(:); ~Indnoncradle(:)],:);
Y = Y([Indcradle(:); Indnoncradle(:)],:);
mask = logical([ones(sum(Indcradle),1);zeros(sum(Indnoncradle),1)]);
opt.unsupervised_setup = struct();
opt.unsupervised_setup.IndY = [opt.vfeatureind.cradle(Indcradle),opt.vfeatureind.free(Indnoncradle)];
opt.unsupervised_setup.Zmask = logical([ones(sum(~Indcradle),1);zeros(sum(~Indnoncradle),1)]);
opt.unsupervised_setup.IndZ = [opt.vfeatureind.cradle(~Indcradle),opt.vfeatureind.free(~Indnoncradle)];
end


