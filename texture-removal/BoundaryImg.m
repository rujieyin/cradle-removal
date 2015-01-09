function [bdImg,lowfreqImg,mask] = BoundaryImg(img,opt)

% this function returns boundary image of the cradling

if ~isfield(opt,'bdwidth')
    opt.bdwidth = 10;
end
if ~isfield(opt,'transform')
    opt.transform = 'shearlet';
end
if ~isfield(opt,'BDdemean')
    opt.BDdemean = 0;
end

if isfield(opt,'bdwidthIn')
    widthIn = opt.bdwidthIn;
else
    widthIn = 6;
end
if isfield(opt,'bdwidthOut')
    widthOut = opt.bdwidthOut;
else
    if isfield(opt,'cutwidth')
       widthOut = max(opt.cutwidth)*2;
    else
    widthOut = 10;
    end
end
if ~isfield(opt,'splitBD')
    opt.splitBD = 0;
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
        [maskin,maskout]= horizontalboundarymask(rot90(mask,-1),widthIn,widthOut,filter);
        mask = rot90(maskin+maskout,1); % mask = rot90(maskin,1);
        
        %         mask = imfilter(mask,filter,'replicate','full');
        %         mask = mask(:,round(bdwidth/2):(round(bdwidth/2)+size(img,2)));
        %         mask = (mask~=0);
    case 'horizontal'
        mask = opt.mask_h;
        mask = double(mask);
        [maskin,maskout] = horizontalboundarymask(mask,widthIn,widthOut,filter);
        mask_cut = maskout.*opt.mask_v;
        mask_cut = imfilter(mask_cut,normpdf(1:2*widthOut,widthOut,widthOut/2),'replicate');
        mask = maskin;
        %         mask = imfilter(mask,filter','replicate','full');
        %         mask = mask(round(bdwidth/2):(round(bdwidth/2)+size(img,1) - 1),:);
        %         mask = abs(mask)./max(abs(mask(:)));
end

if sum(mask(:)) == 0
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
% if opt.rmBDhighfreq
%     if opt.splitBD
%         bdST(:,:,(j>1)&directmask) = threshold_highfreq(ST,middleind,opt.direction);
%     else
    bdST(:,:,(j>1)&directmask) = threshold_highfreq(ST);
%     end
% end
lowfreqImg = inverseShearletTransformSpect(bdST);
if  opt.BDdemean
    lowfreqImg = lowfreqImg - median(lowfreqImg(mask(:) > .5));
end
bdImg = lowfreqImg.*mask;
if strcmp(opt.direction,'horizontal')
    bdST = zeros(sizST);
    bdST(:,:,ind) = ST(:,:,ind);
    ind2 = (j>2) & directmask & (abs(k) <= 1);
    ind2 = ind2 | ((j==2) & directmask & (abs(k) == 0));
    bdST(:,:,ind2) = ST(:,:,ind2);
%     bdST(:,:,(j>1)&directmask) = threshold_highfreq(ST,middleind,opt.direction);
    highfreqImg = inverseShearletTransformSpect(bdST);
    bdImg = bdImg + highfreqImg.*mask_cut;
end
% tmpbd = lowfreqImg(mask>0);
% tmpbd = tmpbd - mean(tmpbd(:));
% bdImg = zeros(size(img));
% bdImg(mask>0) = tmpbd.*mask(mask>0);

    function [bdmaskIn,bdmaskOut] = horizontalboundarymask(mask,widthIn,widthOut,filter)
        bdwidth = max(widthIn,widthOut);
        filter = filter(:);
        mask = double(mask > 0);
        mask = [zeros(1,size(mask,2));mask(1:end-1,:)-mask(2:end,:)];
        filterin = [zeros(widthOut,1);.5;filter(end-widthIn+1:end)];
        filterout = [filter(1:widthOut);.5;zeros(widthIn,1)];
        if sum(abs(mask(:))) == 0
            bdmaskIn = zeros(size(mask));
            bdmaskOut = zeros(size(mask));
            return;
        end
        tmpmaskIn = imfilter(double(mask > 0), filterin,'full','replicate');
        tmpmaskIn = tmpmaskIn(widthIn+1:widthIn+size(mask,1),:);
        bdmaskIn = imfilter(double(mask < 0), flipud(filterin),'full','replicate');
        bdmaskIn = bdmaskIn(widthOut+1:widthOut+size(mask,1),:) + tmpmaskIn;
        tmpmaskOut = imfilter(double(mask > 0), filterout,'full','replicate');
        tmpmaskOut = tmpmaskOut(widthIn+1:widthIn+size(mask,1),:);
        bdmaskOut = imfilter(double(mask < 0), flipud(filterout),'full','replicate');
        bdmaskOut = bdmaskOut(widthOut+1:widthOut+size(mask,1),:) + tmpmaskOut;
    end

    function stnew = threshold_highfreq(ST,middleind,dir)
        stnew = zeros(size(ST));
        %         if strcmp(direction,'horizontal')
        %             directmask = (cone == 'v');
        %         else
        %             directmask = (cone == 'h');
        %         end
        J = max(j);
        for i = 2:(J-1)
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
        stnew(:,:,(j==J) & directmask) = ST(:,:,(j == J) & directmask);
        stnew = stnew(:,:,(j>1)&directmask);
    end

    function stnew = threshold_st(st)
        th = max(abs(st(:)))/2;
        stnew = st.*(abs(st) > th);
    end
end



