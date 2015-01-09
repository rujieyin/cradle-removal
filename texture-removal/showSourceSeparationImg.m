function [] = showSourceSeparationImg(noncradleImg,cradleImg,residualImg,opt,mask,loc1,loc2)
% this function show the result of noncradleImg,cradleImg and residualImg
% Input:  mask:  label as an input of sparse factor model
%         loc1:  location of labeled data
%         loc2:  location of unlabeled data

if length(mask) ~= length(loc1)
    disp('The mask and loc1 do not match in length!')
end

M = max(abs([noncradleImg(:);cradleImg(:);residualImg(:)]));
% noncradleImg = noncradleImg/3/M + .5;
% cradleImg = cradleImg/3/M + .5;
% residualImg = residualImg/3/M + .5;
% 
w = size(noncradleImg,1);
l = size(noncradleImg,2);
switch opt.direction
    case 'vertical'
w0 = opt.vdownsamplesize(1);
l0 = opt.vdownsamplesize(2);
    case 'horizontal'
        w0 = opt.hdownsamplesize(1);
        l0 = opt.hdownsamplesize(2);
end

I = zeros(w0,l0,3);
I = reshape(I,[],3);
ncind = loc1(~mask);
cind = loc1(mask);
I(ncind,1) = 1;
I(cind,2) = 1;
I(loc2,3) = 1;
I = reshape(I,w0,l0,3);
I1 = zeros(w,l,3);
I1(:,:,1) = imresize(I(:,:,1),[w,l],'bilinear');
I1(:,:,2) = imresize(I(:,:,2),[w,l],'bilinear');
I1(:,:,3) = imresize(I(:,:,3),[w,l],'bilinear');
I2 = cat(2,I1,ones(w,10,3),I1,ones(w,10,3),I1);
alpha = [abs(noncradleImg)/M, ones(w,10),abs(cradleImg)/M, ones(w,10),abs(residualImg)/M];
% I2 = cat(2,I1.*repmat(noncradleImg,1,1,3),ones(w,10,3),I1.*repmat(cradleImg,1,1,3),ones(w,10,3),I1.*repmat(residualImg,1,1,3));
figure('name','source separation','number','off');
imagesc([noncradleImg,zeros(w,10),cradleImg,zeros(w,10),residualImg]);
axis image; axis off; colormap gray
hold on
h = imagesc(I2);axis off;axis image;
set(h, 'AlphaData',.5*alpha);
title('noncradle -- cradle -- residual, \rho = , unsupervied')
hold off
