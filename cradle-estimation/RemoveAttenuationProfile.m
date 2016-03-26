function [imgnew,info] = RemoveAttenuationProfile(img, midpoint, angle, width,p1profile,C)
% this function remove the horizontal cradle attenuation using a pre-computed
% cross-section profile of cradle
% inpute: the linear attenuation model used is
%   y(:) = p1(:) .* x(:) + p2(:), p2 = C*(1-p1)
%       p1profile: a matrix where each column is a single profile
%       C:         a vector where each entry is correspond to a coloumn in
%                   p1profile

% compute the coordinates (row indices)
len = size(img,2);

if size(p1profile,1) < size(p1profile,2)
p1profile = p1profile';
end

x = (1:len) - ceil(len/2);% column
y1 = midpoint(1) + round(cos(angle)*x);% row

% rescale the profile if the input profile is not consist with cradle img
p1length = size(p1profile,1); % same as p2profile
p1length = p1length - width*2; % the edge distance of the input profile
d = abs(midpoint(1) - midpoint(2)) + 1; % the edge distance of current cradle
if d ~=  p1length
    width = round(width/p1length*d);
    newlength = d + width*2;
    p1profile = imresize(p1profile,[newlength,size(p1profile,2)]);
end
% compute p2
p2profile = (1-p1profile)*diag(C);

% compute mean
if size(p1profile,2) > 1
    p1profile = mean(p1profile,2);
    p2profile = mean(p2profile,2);
end

info = struct();
info.p1profile = p1profile;
info.C = C;
info.width = width;

% padding
p1profile = [p1profile; ones(size(img,1)-size(p1profile,1),1)];
p2profile = [p2profile; zeros(size(img,1)-size(p2profile,1),1)];


transformSlice = @(x,y)(circshift(p1profile,y-width-1).*x + circshift(p2profile,y-width-1));
imgnew = num2cell(img,1);
imgnew = cellfun(@(x,y)transformSlice(x,y),imgnew,num2cell(y1),'UniformOutput',0);
imgnew = cell2mat(imgnew);
