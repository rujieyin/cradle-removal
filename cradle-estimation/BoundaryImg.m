function [bdImg,lowfreqImg,mask] = BoundaryImg(img,opt)

% this function returns boundary image of the cradling

if ~isfield(opt,'bdwidth')
    opt.bdwidth = 10;
end
if ~isfield(opt,'transform')
    opt.transform = 'shearlet';
end
if ~isfield(opt,'rmBDhighfreq')
    opt.rmBDhighfreq = 1;
end
if ~isfield(opt,'BDdemean')
    opt.BDdemean = 0;
end

if isfield(opt,'bdwidthIn')
    widthIn = opt.bdwidthIn;
else
    widthIn = 4;
end
if isfield(opt,'bdwidthOut')
    widthOut = opt.bdwidthOut;
else
    widthOut = 20;
end
if ~isfield(opt,'splitBD')
    opt.splitBD = 1;
end

filterOut = cdf('normal',1:widthOut,widthOut/2,widthOut/6);
filterIn = 1 - cdf('normal',1:widthIn,widthIn/2,widthIn/6);
filter = [filterOut,1,filterIn];

% if isfield(opt,'bdwidth')
%     bdwidth = opt.bdwidth; % boundary width
% else
%     bdwidth = 40;
% end
% 
% filter = cdf('normal',1:bdwidth,bdwidth/2,bdwidth/6);
% filter = [filter,1,1-filter];
switch opt.direction
    case 'vertical'
        mask = opt.mask_v;
        mask = double(mask);
        [mask,middleind]= horizontalboundarymask(rot90(mask,-1),widthIn,widthOut,filter);
        mask = rot90(mask,1);
        
        %         mask = imfilter(mask,filter,'replicate','full');
        %         mask = mask(:,round(bdwidth/2):(round(bdwidth/2)+size(img,2)));
        %         mask = (mask~=0);
    case 'horizontal'
        mask = opt.mask_h;
        mask = double(mask);
        [mask,middleind] = horizontalboundarymask(mask,widthIn,widthOut,filter);
        %         mask = imfilter(mask,filter','replicate','full');
        %         mask = mask(round(bdwidth/2):(round(bdwidth/2)+size(img,1) - 1),:);
        %         mask = abs(mask)./max(abs(mask(:)));
end

if isempty(middleind)
    bdImg = zeros(size(img));
    lowfreqImg = zeros(size(img));
    return;
end

switch opt.transform
    case 'shearlet'
        ST = shearletTransformSpect(img);
        %     case 'curevelet'
        %         CW = fdct_wrapping(img,1,1,7);
end

[j,k,cone] = arrayfun(@(x)shearletScaleShear(x),1:size(ST,3));
sizST = size(ST);
bdST = zeros(sizST);
if strcmp(opt.direction,'horizontal')
    directmask = (cone == 'v');
else
    directmask = (cone == 'h');
end
ind = (j < 2) & (k == 0) & directmask;
bdST(:,:,ind) = ST(:,:,ind);
if opt.rmBDhighfreq
    if opt.splitBD       
        bdST(:,:,(j>1)&directmask) = threshold_highfreq(ST,middleind,opt.direction);
    else
    bdST(:,:,(j>1)&directmask) = threshold_highfreq(ST);
    end
end
lowfreqImg = inverseShearletTransformSpect(bdST);
if  opt.BDdemean
    lowfreqImg = lowfreqImg - median(lowfreqImg(mask(:) > .5));
end

bdImg = lowfreqImg.*mask;
% tmpbd = lowfreqImg(mask>0);
% tmpbd = tmpbd - mean(tmpbd(:));
% bdImg = zeros(size(img));
% bdImg(mask>0) = tmpbd.*mask(mask>0);

    function [bdmask,middleRind] = horizontalboundarymask(mask,widthIn,widthOut,filter)
        bdwidth = max(widthIn,widthOut);
        mask = double(mask > 0);
        mask = [zeros(1,size(mask,2));mask(1:end-1,:)-mask(2:end,:)];
        mask = abs(mask);
        % exclude wrong signal within boundary
        for i = 1:size(mask,2)
            indup = find(mask(:,i)>0,1,'first');
            inddown = find(mask(:,i)>0,1,'last');
            tmp = zeros(size(mask,1),1);%tmp = zeros(1,size(mask,2));
            tmp(indup) = mask(indup,i);
            tmp(inddown) = mask(inddown,i);
            mask(:,i) = tmp;
        end
        if sum(mask(:)) == 0
            bdmask = zeros(size(mask));
            middleRind = [];
            return;
        end
        bdmask = [zeros(bdwidth,size(mask,2));mask;zeros(bdwidth,size(mask,2))];
        filter = filter(:);
        bdind = find(bdmask > 0);
        [rowind,colind] = ind2sub(size(bdmask),bdind);
        middleRind = round(mean(rowind));
        bdmask = zeros(size(bdmask));
        for ii = 1:length(rowind)
            if rowind(ii) < middleRind
            bdmask(rowind(ii)-widthOut:rowind(ii)+widthIn,colind(ii)) = bdmask(rowind(ii)-widthOut:rowind(ii)+widthIn,colind(ii)) +filter;
            else
            bdmask(rowind(ii)-widthIn:rowind(ii)+widthOut,colind(ii)) = bdmask(rowind(ii)-widthIn:rowind(ii)+widthOut,colind(ii)) + flipud(filter);
             end
       end
        bdmask = bdmask(bdwidth+1:end-bdwidth,:);
    end

    function stnew = threshold_highfreq(ST,middleind,dir)
        stnew = zeros(size(ST));
        %         if strcmp(direction,'horizontal')
        %             directmask = (cone == 'v');
        %         else
        %             directmask = (cone == 'h');
        %         end
        J = max(j);
        for i = 2:J
            indmask = (j == i) & directmask;
            if nargin == 1
            stnew(:,:,indmask) = threshold_st(ST(:,:,indmask));
            else
                switch dir
                    case 'horizontal'
                        stnew(1:middleind,:,indmask) = threshold_st(ST(1:middleind,:,indmask));
                        stnew(middleind:end,:,indmask) = threshold_st(ST(middleind:end,:,indmask));
                    case 'vertical'
                        stnew(:,1:middleind,indmask) = threshold_st(ST(:,1:middleind,indmask));
                        stnew(:,middleind:end,indmask) = threshold_st(ST(:,middleind:end,indmask));
                end
            end
        end
        stnew = stnew(:,:,(j>1)&directmask);
    end
    
    function stnew = threshold_st(st)
        th = max(abs(st(:)))/2;
        stnew = st.*(abs(st) > th);
    end
end



