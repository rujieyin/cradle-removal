function [subimg,subind] = subpatch4HorBdryRm(img,hinfo,vinfo)
% this function get subimgs of size no larger than 512-by-512, which
% containe horizontal boundary of cradling need to be removed

% Input: vinfo:    connected vertical segmentation info

nh = length(hinfo);
L = min(512,min(size(img))); % subimg size
D = 400; % biggest gap b/t upper and lower cut

if isfield(vinfo{1},'upcutind')
    % function of getting upper and lower cut index of xth horizontal
    % cradle
    cutind1 = @(x)min(cellfun(@(y)y.downcutind(x),vinfo));
    cutind2 = @(x)max(cellfun(@(y)y.upcutind(x+1),vinfo));
    cut = 1;
else
    cut = 0;
end
% function of getting upper and lower bdry index of xth horizontal cradle
bd1 = @(x)min(hinfo{x}.lp(1),hinfo{x}.rp(1));
bd2 = @(x)max(hinfo{x}.ld(1),hinfo{x}.rd(1));

subimg = cell(2*nh,1);
subind = cell(2*nh,1);

for i = 1:nh
    if cut
        ind1 = cutind1(i);
        ind2 = cutind2(i);
    else
        ind1 = bd1(i);
        ind2 = bd2(i);
    end
    width = ind2 - ind1;
    center = round((ind1+ind2)/2);
    if width < D
        ind1 = center - L/2;
        subimg{2*i-1} = img(ind1:(ind1 + L -1),:);
        [subimg{2*i - 1},~,~,~,subind{2*i-1}] = img2subpatch(subimg{2*i-1},L);
        ncol = length(subind{2*i-1});
        subind{2*i-1} = [ind1*ones(1,ncol);subind{2*i-1}];
        subimg{2*i} = [];
    else
        ind1 = max(1,ind1 - L/2);
        ind2 = min(ind2 - L/2,size(img,2)-L+1);
        
        subimg{2*i-1} = img(ind1:(ind1 + L-1),:);
        [subimg{2*i - 1},~,~,~,subind{2*i-1}] = img2subpatch(subimg{2*i-1},L);
        ncol = length(subind{2*i-1});
        subind{2*i-1} = [ind1*ones(1,ncol);subind{2*i-1}];
        
        subimg{2*i} = img(ind2:(ind2 + L-1),:);
        [subimg{2*i},~,~,~,subind{2*i}] = img2subpatch(subimg{2*i},L);
        subind{2*i} = [ind2*ones(1,ncol);subind{2*i}];
    end
end        

subimg = subimg(~cellfun(@isempty,subimg));
subind = subind(~cellfun(@isempty,subind));
