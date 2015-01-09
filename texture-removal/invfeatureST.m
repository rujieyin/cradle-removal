function img = invfeatureST(st,opt,relocate,loc)

% this function compute the inverse of shearlet feature vector to image
% domain
% Input:   st:  Npts - by - p ST feature vector, either is concatenation of
%               [cradle feature; non-cradle feature] , or Npts is the size of
%               prod(opt.downsamplesize), where relocation is not needed
%               OR Nsample - by - 1 cell, each cell is a ST feature vector.
%    relocate:  logical, whether relocation is needed


% -- average over samples -- %
if iscell(st)
    siz = size(st{1});
    st = cellfun(@(x)x(:),st,'UniformOutput',0);
    st = cell2mat(st(:)');
    st = mean(st,2);
    st = reshape(st,siz);
end


switch opt.direction
    case 'vertical'
        ind = [opt.vfeatureind.cradle; opt.vfeatureind.free];
        siz = opt.vdownsamplesize;
        angle = opt.vangleind;
    case 'horizontal'
        ind = [opt.hfeatureind.cradle; opt.hfeatureind.free];
        siz = opt.hdownsamplesize;
        angle = opt.hangleind;
end

l = length(opt.ST.j);
ST = zeros(siz(1)*siz(2),length(angle));
if nargin < 3 || relocate
    if nargin > 3
        ind = loc;
    end
    % -- if the index is wrong, return all zero image -- %
    if length(ind) ~= size(st,1)
        img = zeros(size(opt.imgsize));
        return;
    end
    ST(ind,:) = st;
else
    ST = st;
end
ST = reshape(ST,siz(1),siz(2),[]);
STnew = zeros(opt.imgsize(1),opt.imgsize(2),l);
downsamplerate = opt.imgsize./siz;
n = round(downsamplerate(1));
m = round(downsamplerate(2));
if angle(end) == angle(end-1)
    STnew(1:n:end,1:m:end,angle(1:end-2)) = ST(:,:,1:end-2);
    if strcmp(opt.direction,'vertical')
        STnew(1:n:end,1:m:end,angle(end)) = ST(:,:,end-1)/2;
        STnew(1:n:end,(1+m/2):m:end,angle(end)) = ST(:,:,end)/2;
    else
        STnew(1:n:end,1:m:end,angle(end)) = ST(:,:,end-1)/2;
        STnew((1+n/2):n:end,1:m:end,angle(end)) = ST(:,:,end)/2;
    end
else
    STnew(1:n:end,1:m:end,angle) = ST;
end
img = inverseShearletTransformSpect(STnew);
img = img*m*n; % rescale img

