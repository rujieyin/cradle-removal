function [region,distscore,trans] = featuresimilarityregion(feature,opt,img)

if strcmp(opt.direction,'horizontal')
    cradleftr = feature.hcradle;
    noncradftr = feature.nonhcradle;
%     cradleftr = feature.hcradle_normalized;
%     noncradftr = feature.nonhcradle_normalized;
    cradleind = opt.hfeatureind.cradle;
    noncradind = opt.hfeatureind.free;
    siz = opt.hdownsamplesize;
else
    cradleftr = feature.vcradle;
    noncradftr = feature.nonvcradle;
%     cradleftr = feature.vcradle_normalized;
%     noncradftr = feature.nonvcradle_normalized;
    cradleind = opt.vfeatureind.cradle;
    noncradind = opt.vfeatureind.free;
    siz = opt.vdownsamplesize;
end

% dist_thresh = zeros(length(cradleind),1);
% dist2noncrad = zeros(length(cradleind),length(noncradind));%each row is the dist 2 all noncrad pts
% normcradleftr = sqrt(sum(cradleftr.^2,1));
normnoncradftr = sqrt(sum(noncradftr.^2,1));
distscore = zeros(length(cradleind),length(noncradind));
for i = 1:length(cradleind)
    distscore(i,:) = sum((noncradftr - repmat(cradleftr(:,i),1,length(noncradind))).^2,1)./normnoncradftr/norm(cradleftr(:,i));
    distscore(i,:) = exp(-distscore(i,:));%abs(cradleftr(:,i)'*noncradftr);
%     neighind = findneighbour(cradleind(i),cradleind,siz);
%     if ~isempty(neighind)
%         %squared distance threshold
%         dist_thresh(i) = max(sum((cradleftr(:,neighind) - repmat(cradleftr(:,i),1,length(neighind))).^2,1)./normcradleftr(neighind));
%         dist_thresh(i) = dist_thresh(i)/norm(cradleftr(:,i));
%         dist2noncrad(i,:) = sum((noncradftr - repmat(cradleftr(:,i),1,length(noncradind))).^2,1)./normnoncradftr;
%         distscore(i,:) = dist_thresh(i)./dist2noncrad(i,:);
%     else
%         distscore(i,:) = -1;
%     end
end

maxdistscore = max(distscore,[],2);
maxdistscore = (maxdistscore - min(maxdistscore))/(1- min(maxdistscore));
tmp = zeros(siz);
tmp(cradleind) = maxdistscore;
figure;imagesc(tmp);axis off;title('distance score');

candidate = find(maxdistscore > .5);
accepted = [];
rejected = [];
region = cell(1);
trans = cell(1);
nregion = 0;
fullsiz = size(img);
if nargin > 2
    img = (img - min(img(:)))/(max(img(:)) -  min(img(:)));
    img = repmat(img,[1,1,3]);
%     img = repmat(imresize(img,siz),[1,1,3]);
    tmpimg = zeros(fullsiz);
    tmpimg = repmat(tmpimg,[1,1,3]);
    tmp = zeros(siz);
    tmp(cradleind) = 1;
    tmpimg(:,:,1) = imresize(tmp,fullsiz,'method','nearest');
    tmp = zeros(siz);
    tmp(noncradind) = 1;
    tmpimg(:,:,3) = imresize(tmp,fullsiz,'method','nearest');
    alpha = zeros(siz);
end
figure;image(img);axis off;hold on
h = image(tmpimg);
set(h,'AlphaData',imresize(alpha,fullsiz));
drawnow;
while ~isempty(candidate)
    tmpalpha = alpha;
    tmpbeta = zeros(size(alpha));
    [~,center] = max(maxdistscore(candidate));
    center = candidate(center);
    localvisited = center;
    localaccepted = center;
    activebd = center;
    [~,cocenter] = max(distscore(center,:));
    [m,n] = ind2sub(siz,cradleind(center));
    [com,con] = ind2sub(siz,noncradind(cocenter));
    if nargin > 2
        tmpalpha(m,n) = maxdistscore(center)*.3; 
        tmpbeta(com,con) = tmpalpha(m,n);
    end
    deltasub = [com,con] - [m,n];
    while ~isempty(activebd)
        for k = 1:length(activebd)
            tmpcenter = activebd(k);
            newactivebd = [];
            neighind = findneighbour(cradleind(tmpcenter),cradleind(candidate),siz);
            neighind = candidate(neighind);
            neighind = setdiff(neighind,localvisited);
            if ~isempty(neighind)
                localvisited = [localvisited, neighind(:)'];
                for i = 1:length(neighind)
                    tmpneighind = neighind(i);
                    [m,n] = ind2sub(siz,cradleind(tmpneighind));
                    co = [m,n] + deltasub;
                    com = co(1);
                    con = co(2);
                    if com <= siz(1) && con <= siz(2) && com >=1 && con >= 1
                        coneighind = sub2ind(siz,com,con);
                        [~,coneighind] = intersect(noncradind,coneighind);
                        if distscore(tmpneighind,coneighind)/maxdistscore(tmpneighind) > .7
                            newactivebd = [newactivebd, tmpneighind];
                            if nargin > 2
                                tmpalpha(m,n) = distscore(tmpneighind,coneighind)*.3;
                                tmpbeta(com,con) = tmpalpha(m,n);
%                                 tmpimg(m,n,1) = 1;
%                                 tmpimg(m,n,2:3) = 0;
%                                 tmpimg(com,con,3) = 1;
%                                 tmpimg(com,con,1:2) = 0;
%                                 scoremap(m,n,1) = distscore(tmpneighind,coneighind);
%                                 image([tmpimg,scoremap]);axis off;drawnow;
                            end
                        else
                            rejected = [rejected,tmpneighind];
                        end
                    else
                        rejected = [rejected, tmpneighind];
                    end
                end
            end
        end
        newactivebd = unique(newactivebd);
        activebd = newactivebd;
        localaccepted = [localaccepted,activebd];
    end
    if length(localaccepted) > 4
        if nargin >2
            alpha = tmpalpha;
            set(h,'AlphaData',imresize(alpha+tmpbeta,fullsiz,'method','nearest'));drawnow;
        end
        nregion = nregion + 1;
        region{nregion} = localaccepted;
        trans{nregion} = deltasub;
        accepted = [accepted,localaccepted];
    else
        rejected = [rejected,center];
    end
    
    candidate = setdiff(candidate,union(rejected,accepted));
    
end
        



function neighind = findneighbour(ind,indset,siz)
[mm,nn] = ind2sub(siz,ind);
if mm < siz(1) && mm > 1
    neighind = [mm-1,nn;mm+1,nn];
else
    neighind = [];
end
if nn < siz(2) && nn > 1
    neighind = [neighind;mm,nn-1;mm,nn+1];
end
if ~isempty(neighind)
neighind = sub2ind(siz,neighind(:,1),neighind(:,2));
[~,~,neighind] = intersect(neighind,indset);
end
end

end
