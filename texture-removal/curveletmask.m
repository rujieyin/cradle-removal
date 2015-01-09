function cwmask = curveletmask(mask,cw)
% this function compute the mask of curvelet coefficients based on the
% given mask of original image size

% Input:   mask:    size(mask) = size(img)
%            cw:    curvelet transform coefficients
%                   or
%                   a vector of curvelet transform parameters [~,isreal,
%                   finest, nscale] 
% Output:   cwmask:  mask of curvelet coefficients, which has the same
%                    structure of cw obtained by fdct_wrapping

siz = size(mask);
if isnumeric(cw)
    cw = fdct_wrapping(ones(siz),cw(1),cw(2),cw(3));
end

cwmask = cell(size(cw));
for i = 1:length(cwmask)
    cwmask{i} = cell(size(cw{i}));
    for j = 1:length(cwmask{i})
        cwmask{i}{j} = imresize(mask,'OutputSize',size(cw{i}{j}),'Method','bilinear');
    end
end
