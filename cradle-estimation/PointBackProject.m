function [upline,downline,l,r] = PointBackProject(b,angle,siz,range)

% this function compute a point with bin index b back project along angle
% in radon transform onto an image of certain size. At each column, the
% point is split into half and rounded to upper and lower index.
% Notice here we only consider the horizontal back projection, i.e. angle
% is around 90 degree

% Input:  b:              bin index
%         angle:          angle of projection in radon transform
%         siz:            size of the image
%         range:          (optional) specific range of columns to compute
%                         the back projection
% Output: upline:         index of upper line of the back projection 
%         downline:       index of lower line of the back projection
%         l:              index of the starting column of vector upline &
%                         downline, if the vectors are empty, L = 0
%         r:              index of ending column

b = b(:)';
if nargin < 4
    range = [1,siz(2)];
end
L = siz(1)+1;
M = (siz(2)+1)/2;
theta = angle/180*pi;

if abs(angle - 90) < .1
    upline = ones(siz(2),1)*ceil(L/2-b);
    downline = ones(siz(2),1)*floor(L/2-b);
    l = range(1);
    r = range(2);
else
    
    if angle < 90
        l = max([ceil(tan(theta)*(.5+1/2/sin(theta)+b/sin(theta)-L/2)+M),range(1)]);
        r = min([floor(tan(theta)*(L/2-1/2/sin(theta)+b/sin(theta))+M),range(2)]);
    else
        r = min([floor(tan(theta)*(.5+1/2/sin(theta)+b/sin(theta)-L/2)+M),range(2)]);
        l = max([ceil(tan(theta)*(L/2-1/2/sin(theta)+b/sin(theta))+M),range(1)]);
    end
    if l > M | r < M
        l = 0;
        r = 0;
        upline = [];
        downline = [];
        return;
    end
    x = L/2-b/sin(theta)+[l-M:r-M]*cot(theta);
    upline = ceil(x - 1/2/sin(theta));
    downline = floor(x+1/2/sin(theta));
end
