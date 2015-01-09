i = 4;
j = 7;
mask_v = blockmask{i,j}.v;
mask_h = blockmask{i,j}.h;
imgtmp = subpatch{i,j};

opt = struct();
opt.direction = 'horizontal';
% mask_v = zeros(size(img));
% mask_v(:, ) = 1;
% mask_h = zeros(size(img));
% mask_h( ,:) = 1;
opt.imgsize = size(imgtmp);
opt.mask_v = mask_v;
opt.mask_h = mask_h;
texturecw = fdct_wrapping(texture,1,1,7);
% non-prob direct masking model
[hcradlecw,tmp1,opt] = angleselection2(texturecw,mask_h,opt);
opt.direction = 'vertical';
[vcradlecw,tmp2,opt] = angleselection2(texturecw,mask_v,opt);
texturenew = tmp1+tmp2;
% compute probability of texture within cradling region
opt.direction = 'horizontal';
hcw = angleselection2(texturecw,ones(opt.imgsize),opt);
[hfeature,opt] = getCWfeature(hcw,opt);
opt.direction = 'vertical';
vcw = angleselection2(texturecw,ones(opt.imgsize),opt);
[vfeature,opt] = getCWfeature(vcw,opt);