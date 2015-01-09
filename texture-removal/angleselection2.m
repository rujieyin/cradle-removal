function [cwnew,imgnew,opt] = angleselection2(cw,mask,opt)

% Input:  mask:    it could be either a matrix of the same size as img, or
%                  it could be a structure which is the same as curvelet
%                  coefficients cw

if ~iscell(mask)
    if islogical(mask)
        mask = double(mask);
    end
    cwmask = curveletmask(mask,cw);
else
    cwmask = mask;
end

switch opt.direction
    case 'horizontal'
        angleind = cellfun(@(x)[1:floor(length(x)/4),floor(length(x)/2)+1:floor(length(x)/4)*3],cw,'UniformOutput',0);
        opt.hangleind = angleind;
    case 'vertical'
        angleind = cellfun(@(x)[floor(length(x)/4)+1:floor(length(x)/2),floor(length(x)/4)*3+1:floor(length(x)/2)*2],cw,'UniformOutput',0);
        opt.vangleind = angleind;
end

cwnew = cell(size(cw));

for i = 1:length(angleind)
    cwnew{i} = cellfun(@(x)zeros(size(x)),cw{i},'UniformOutput',0);
    cwnew{i}(angleind{i}) = cw{i}(angleind{i});
end

cwnew = ApplyMask2CW(cwmask,cwnew);
% if nargin == 3
%     if isfield(opt,'direction') 
%         cwnew = angleselection(cwnew,opt);
%     end
% else
%     opt = struct();
% end
score = cellfun(@(x)sum(abs(x(:))),cwnew{end});
[~,ind] = sort(score,'descend');
energy = cellfun(@(x)sum(abs(x(:)).^2),cwnew{end});
energy = energy(ind);
energy = cumsum(energy);
energy = energy./energy(end);
i = 1;
while energy(i) < .9
    i = i+1;
end
ind = ind(1:i);
ind = sort(ind,'ascend');
nangles = cellfun(@(x)length(x),cwnew);
angleind = cell(length(nangles),1);
angleind{1} = [];
for i = 2:length(angleind)
    angleind{i} = ceil(ind/nangles(end)*nangles(i));
    angleind{i} = unique(angleind{i});
end
tmp = cell(size(cw));
% for i = 1:length(cw)
%     tmp{i} = cell(size(cw{i}));
%     tmp{i} = cellfun(@(x)zeros(size(x)),cw{i},'UniformOutput',0);
%     tmp{i}(angleind{i}) = cw{i}(angleind{i});
% end
% cwnew = ApplyMask2CW(cwmask,tmp);
for i = 1:length(cw)
    tmp{i} = cell(size(cwnew{i}));
    tmp{i} = cellfun(@(x)zeros(size(x)),cwnew{i},'UniformOutput',0);
    tmp{i}(angleind{i}) = cwnew{i}(angleind{i});
end

if isfield(opt,'imgsize')
    imgnew = ifdct_wrapping(cwnew,opt.curveletisreal,opt.imgsize(1),opt.imgsize(2));
    if isfield(opt,'display') && opt.display
        figure;imshow(imgnew,[]);
    end
end

if isfield(opt,'direction')
    if strcmp(opt.direction,'horizontal')
        opt.hangleind = angleind;
    else
        opt.vangleind = angleind;
    end
else
    opt.angleind = angleind;
end