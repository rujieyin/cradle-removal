function [st,opt] = shearlet_angleselection(ST,j,k,cone,opt)
% this function select only horizontal or vertical directions of shearlets

maxj = max(j);
angle = zeros(1,length(j));
if isfield(opt,'maxangle')
    maxangle = opt.maxangle;
else
    maxangle = 5;
end

for i = 1 : maxj
    angle(j == i) = (abs(k(j == i)) < 2^i/45*maxangle);
%     angle(j == i) = (abs(k(j == i)) < 2^i/2);
end

switch opt.direction
    case 'vertical'
        angle = angle & (cone == 'h'); % periodic in horizontal
        opt.vangleind = angle;
    case 'horizontal'
        angle = angle & (cone == 'v'); % periodic in veritcal
        opt.hangleind = angle;
end

opt.ST.j = j;
opt.ST.k = k;
opt.ST.cone = cone;

st = zeros(size(ST));
st(:,:,angle) = ST(:,:,angle);
