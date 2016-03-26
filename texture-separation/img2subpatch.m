function [subpatch,orow,ocol,rowind,colind] = img2subpatch(img,patchsize,opt)

% to resemble patches to image, use block2img(subpatch,rowind,colind)
% edited by Rachel Yin March 19, 2015. remove the default option of padding
% image

imgsiz = size(img);
if nargin < 2
    patchsize = 512;
elseif length(patchsize) == 1
        patchsize = ones(2,1)*patchsize;
end
% if nargin < 3
%     opt = 'add';
% end

nrow= ceil(imgsiz(1)/patchsize(1));
if nrow == 1
    orow = 0;
else
    orow = (nrow*patchsize(1) - imgsiz(1))/(nrow-1);
end

% if strcmp(opt,'add')
%     orow = floor(orow);
%     img = [img;repmat(img(end,:),patchsize(1)*nrow - (nrow-1)*orow - imgsiz(1),1)];
% else
orow = ceil(orow);
% img = img(1: (patchsize(1)*nrow - (nrow-1)*orow),:);
% end
rowind = 1:(patchsize(1)-orow):(1+(nrow-1)*(patchsize(1)-orow));
rowind(end) = imgsiz(1)+1-patchsize(1); % fix the last row index
% imgsiz = size(img);

ncol = ceil(imgsiz(2)/patchsize(2));
if ncol == 1
    ocol = 0;
else
    ocol = (ncol*patchsize(2) - imgsiz(2))/(ncol-1);
end

% if strcmp(opt,'add')
%     ocol = floor(ocol);
%     img = [img repmat(img(:,end),1,patchsize(2)*ncol - (ncol-1)*ocol - imgsiz(2))];
% else
    ocol = ceil(ocol);
%     img = img(:,1:(patchsize(2)*ncol - (ncol-1)*ocol));
% end
colind = 1:(patchsize(2)-ocol):(1+(ncol-1)*(patchsize(2)-ocol));
colind(end) = imgsiz(2)+1-patchsize(2);
% imgsiz = size(img);

subpatch = cell(nrow,ncol);
for i = 1:nrow
    for j = 1:ncol
        subpatch{i,j} = img(rowind(i):(rowind(i) + patchsize(1)-1),colind(j):(colind(j)+patchsize(2)-1));
    end
end

