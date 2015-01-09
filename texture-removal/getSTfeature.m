function [feature,opt,ST,img] = getSTfeature(ST,j,opt)

% this function extract feature vector from shearlet coefficients  STnew

% -- exclude low frequency feature -- %
% for vertical pattern (support in horizontal direction in spectrum),
% downsample rate is 8 - by - 2, for horizontal pattern (support in
% vertical dirction in spectrum) downsample rate is 2 - by - 8, which is
% checked for size 512 - by - 512, guaranteeing exact reconstruction with a
% modification of scale 16.

k = opt.ST.k;
cone = opt.ST.cone;

switch opt.direction
    case 'horizontal'
        if islogical(opt.hangleind)
            if isfield(opt,'lowfreqfeature') && opt.lowfreqfeature
                opt.hangleind = find(opt.hangleind & (j>0));
            else
                opt.hangleind = find(opt.hangleind & (j > 1));
            end
        end
        ind2split = find((j == max(j)) & (k== 0) & (cone == 'v'));
        opt.hangleind(opt.hangleind == ind2split) = [];
        st = ST(1:2:end,1:8:end,opt.hangleind);
        opt.hangleind = [opt.hangleind, ind2split,ind2split];
        st = cat(3,st,ST(1:2:end,1:8:end,ind2split),ST(2:2:end,1:8:end,ind2split));
        mask = imresize(opt.mask_h,size(st(:,:,1)),'bilinear');
        mask = logical(mask(:));
        opt.hdownsamplesize = size(st(:,:,1));
    case 'vertical'
        if islogical(opt.vangleind)
            if isfield(opt,'lowfreqfeature') && opt.lowfreqfeature
                opt.vangleind = find(opt.vangleind & (j>0));
            else
                opt.vangleind = find(opt.vangleind & (j>1));
            end
        end
        ind2split = find((j == max(j)) & (k == 0) & (cone == 'h'));
        opt.vangleind(opt.vangleind == ind2split) = [];
        st = ST(1:8:end,1:2:end,opt.vangleind);
        opt.vangleind = [opt.vangleind,ind2split,ind2split];
        st = cat(3,st,ST(1:8:end,1:2:end,ind2split),ST(1:8:end,2:2:end,ind2split));
        mask = imresize(opt.mask_v,size(st(:,:,1)),'bilinear');
        mask = logical(mask(:));
        opt.vdownsamplesize = size(st(:,:,1));
end
if nargout > 2
    if strcmp(opt.direction,'vertical')
        ST(:,:,setdiff(1:size(ST,3),opt.vangleind)) = 0;
    else
        ST(:,:,setdiff(1:size(ST,3),opt.hangleind)) = 0;
    end
    img = inverseShearletTransformSpect(ST);
end
st = reshape(st,[],size(st,3));
energy = sum(st.^2,2);
[energy, ind] = sort(energy,'descend');
tmp = cumsum(energy);
cutind = length(ind);
% cutind = find(tmp/tmp(end) > .95,1);
indmask = zeros(size(ind));
indmask(ind(1:cutind)) = 1;
feature = struct();
switch opt.direction
    case 'horizontal'
        opt.hfeatureind.cradle = find(mask & indmask);
        opt.hfeatureind.free = find(~mask & indmask);
    case 'vertical'
        opt.vfeatureind.cradle = find(mask & indmask);
        opt.vfeatureind.free = find(~mask & indmask);
end
feature.cradle = st(mask & indmask,:);
feature.free = st(~mask & indmask,:);

