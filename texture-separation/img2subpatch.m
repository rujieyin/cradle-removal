function [subpatch,orow,ocol,rowind,colind] = img2subpatch(img,patchsize,opt)

imgsiz = size(img);
if nargin < 2
    patchsize = 512;
elseif length(patchsize) == 1
        patchsize = ones(2,1)*patchsize;
end
if nargin < 3
    opt = 'add';
end

nrow= ceil(imgsiz(1)/patchsize(1));
orow = (nrow*patchsize(1) - imgsiz(1))/(nrow-1);

if strcmp(opt,'add')
    orow = floor(orow);
    img = [img;repmat(img(end,:),patchsize(1)*nrow - (nrow-1)*orow - imgsiz(1),1)];
else
orow = ceil(orow);
img = img(1: (patchsize(1)*nrow - (nrow-1)*orow),:);
end
rowind = 1:(patchsize(1)-orow):(1+(nrow-1)*(patchsize(1)-orow));
imgsiz = size(img);

ncol = ceil(imgsiz(2)/patchsize(2));
ocol = (ncol*patchsize(2) - imgsiz(2))/(ncol-1);

if strcmp(opt,'add')
    ocol = floor(ocol);
    img = [img repmat(img(:,end),1,patchsize(2)*ncol - (ncol-1)*ocol - imgsiz(2))];
else
    ocol = ceil(ocol);
    img = img(:,1:(patchsize(2)*ncol - (ncol-1)*ocol));
end
colind = 1:(patchsize(2)-ocol):(1+(ncol-1)*(patchsize(2)-ocol));
imgsiz = size(img);

subpatch = cell(nrow,ncol);
for i = 1:nrow
    for j = 1:ncol
        subpatch{i,j} = img(rowind(i):(rowind(i) + patchsize(1)-1),colind(j):(colind(j)+patchsize(2)-1));
    end
end

