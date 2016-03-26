function [subimg,subind] = subpatch4VerBdryRm(img,hinfo,vinfo)
% similar function to subpatch4HorBdryRm

% Input: vinfo:  connected vertical segmentation info

nv = length(vinfo);
L = min(512,min(size(img))); % subimg size
D = 400;

bd1 = @(x)min(vinfo{x}.lp(2),vinfo{x}.ld(2));
bd2 = @(x)max(vinfo{x}.rp(2),vinfo{x}.rd(2));

subimg = cell(1,2*nv);
subind = cell(1,2*nv);

for i = 1:nv
    ind1 = bd1(i);
    ind2 = bd2(i);
    width = ind2 - ind1;
    if width < D
        ind1 = max(1, round((ind1+ind2)/2 - L/2));
        subimg{2*i-1} = img(:,ind1:(ind1 + L -1));
        [subimg{2*i-1},~,~,subind{2*i-1}] = img2subpatch(subimg{2*i-1},L);
        nrow = length(subind{2*i-1});
        subind{2*i-1} = [subind{2*i-1};ind1*ones(1,nrow)];
        subimg{2*i} = [];
    else
        ind1 = max(1,ind1 - L/2);
        ind2 = min(ind2 - L/2,size(img,2)-L+1);
        
        subimg{2*i-1} = img(:,ind1:(ind1+L-1));
        [subimg{2*i-1},~,~,subind{2*i-1}] = img2subpatch(subimg{2*i-1},L);
        nrow = length(subind{2*i-1});
        subind{2*i-1} = [subind{2*i-1};ind1*ones(1,nrow)];
        
        subimg{2*i} = img(:,ind2:(ind2 + L-1));
        [subimg{2*i},~,~,subind{2*i}] = img2subpatch(subimg{2*i},L);
        subind{2*i} = [subind{2*i};ind2*ones(1,nrow)];
    end
end

subimg = subimg(~cellfun(@isempty,subimg));
subind = subind(~cellfun(@isempty,subind));