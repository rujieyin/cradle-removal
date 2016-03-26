function [subimgnew,vinfo,I] = PostRmCradleProfile(img,vinfo,i1,j1,i2,j2,I)
% this function remove the single vcradle (i1,j1) after the first round of
% RmRegularVerticalCradleSegmentation using the profile of vcradle (i2,j2)
% Input: vinfo: output of RmRegularVerticalCradleSegmentation
%        i1,j1: index of vcradle to be removed
%        i2,j2: index of vcradle whose profile will be used
%        img:   the image where hcradle has been successfully removed
%        I:     the output image of RmRegularVerticalCradleSegmentation
% Output: imgnew: the subimg where vcradle is removed

info = vinfo{i1,j1};
rowind = info.rowind;
colind = info.colind;
% the target subimg
subimg = img(rowind(1):rowind(2),colind(1):colind(2));
angle = info.angle/180*pi;
% mask 
mask = vinfo{i1,j1}.cradleimg> 1;

% the profile to be used
p1profile = vinfo{i2,j2}.p1profile;
C = vinfo{i2,j2}.C;
width = vinfo{i2,j2}.width;

% == main function == %
[subimgnew, info] = RemoveAttenuationProfile(rot90(subimg,-1),info.midpoint,angle,...
    width, p1profile, C);
subimgnew = rot90(subimgnew,1);

subimgnew = subimg.*(~mask) + subimgnew.*mask;

% update vinfo
vinfo{i1,j1}.C = info.C;
vinfo{i1,j1}.p1profile = info.p1profile;
vinfo{i1,j1}.width = info.width;

if nargin > 6
    I(rowind(1):rowind(2),colind(1):colind(2) ) = subimgnew;
else
    I = [];
end
