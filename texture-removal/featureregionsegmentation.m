function [region] = featureregionsegmentation(opt,img,cw)
if strcmp(opt.direction,'horizontal')
    ind = opt.hangleind;
    mask = opt.mask_h;
else
    ind = opt.vangleind;
    mask = opt.mask_v;
end
fullsiz = size(img);
if nargin > 2
    img = (img - min(img(:)))/(max(img(:)) -  min(img(:)));
    img = repmat(img,[1,1,3]);
else
    img = zeros([opt.imgsize,3]);
end
alpha = zeros(fullsiz);

region = struct('active',cell(1),'fixed',cell(1));
nregion = 0;

% red for cradling part; blue for noncradling part
tmpimg = zeros(fullsiz);
tmpimg = repmat(tmpimg,[1,1,3]);
tmp = mask;
tmpimg(:,:,1) = tmp;
tmp = 1 - mask;
tmpimg(:,:,3) = tmp;

scalesiz = zeros(length(ind),2);
for i = 1:length(ind)
    if ~isempty(ind{i})
        scalesiz(i,1) = min(arrayfun(@(x)size(cw{i}{x},1),ind{i}));
        scalesiz(i,2) = min(arrayfun(@(x)size(cw{i}{x},2),ind{i}));
    end
end


for i = 5:length(ind)
    if ~isempty(ind{i})
        % compute candidates
        siz = scalesiz(i,:);        
        energy = arrayfun(@(x)reshape(imresize(cw{i}{x},siz),[],1),ind{i},'UniformOutput',0);
        energy = reshape(cell2mat(energy(:)'),[siz,length(energy)]);
        energy = sum(energy.^2,3);
        [sortenergy,energyind] = sort(reshape(energy,[],1),'descend');
        ncandidate = find(cumsum(sortenergy)/sum(sortenergy) > .7,1);
        energyind = energyind(1:ncandidate);
        if ~isempty(region.fixed)
            % rescale region.fixed, region.active
            region.fixed = cellfun(@(x)finescaleind(x,i-1,i),region.fixed,'UniformOutput',0);
            fixed = cell2mat(region.fixed(:)');
               region.active = cellfun(@(x)finescaleind(x,i-1,i),region.active,'UniformOutput',0);
         % remove candidates in fixed region
            energyind = setdiff(energyind,fixed);
        else
            fixed = [];
        end
        for j = 1:nregion
            [region.active{j},rmind] = intersect(energyind,region.active{j});
            flag = arrayfun(@dominatescale,region.active{j});
            region.fixed{j} = [region.fixed{j},region.active{j}(flag)];
            region.active{j} = region.active{j}(~flag);            
            energyind(rmind) = [];
        end
        for j = 1:length(energyind)
            addnewregion(energyind(j));
        end

% %         tmpactive = cell(nregion,1);
%         for j = 1:ncandidate
%             [flag,regionind] = newregion(energyind(j));
%             if flag == 1
%                 addnewregion;
%             elseif flag == 0
%                 refineoldregion(regionind,energyind(j));
%             end
%         end
%         region.active(1:length(tmpactive)) = tmpactive;
        rmemptyregion;
       displayregion; 
    end
end



    function addnewregion(currentind)
        % add the current block as a new region
        nregion = nregion + 1;
%         tmpind = zeros(siz);
%         tmpind(energyind(j)) = 1;
%         tmpind = imresize(tmpind,fullsiz,'method','nearest');
%         tmpind = find(tmpind > 0);
%         tmpind = tmpind(:)';
        if dominatescale(currentind)
            region.fixed{nregion} = currentind;
            region.active{nregion} = [];
        else
            region.active{nregion} = currentind;
            region.fixed{nregion} = [];
        end
    end

%     function refineoldregion(regionind,currentind)
%         % apply to blocks that contained in some region's active part
%         liftind = finescaleind(currentind)';
%         if dominatescale(currentind)
%             region.fixed{regionind} = [region.fixed{regionind},liftind];
%         else
%             tmpactive{regionind} = [tmpactive{regionind},liftind];
%         end
%     end


%     function [flag,regionind] = newregion(ind)
%         % decide whether the block creates a new region, contained in old
%         % region's fixed part ( then just ignore it ) or contained in old
%         % region's active part
%         [m,n] = ind2sub(siz,ind);
%         m1 = floor(m/siz(1)*fullsiz(1));
%         n1 = floor(n/siz(2)*fullsiz(2));
%         m1 = max([1,min([m1,fullsiz(1)])]);
%         n1 = max([1,min([n1,fullsiz(1)])]);
%         
%         ind1 = sub2ind(fullsiz,m1,n1);
%         % [m1,n1] the subindex of the bottomright corner, only need to
%         % check this pt due to the decrease of size of block, or
%         % equivalently the increase of resolution
%         flag = 1;
%         regionind = 0;
%         for ii = 1:nregion
%             if ~isempty(intersect(region.fixed{ii},ind1))
%                 flag = -1;
%                 regionind = ii;
%                 break;
%             elseif ~isempty(intersect(region.active{ii},ind1))
%                 flag = 0;
%                 regionind = ii;
%             end
%         end
%     end

    function displayregion
        
        fixedmask = zeros(fullsiz);
        fixedmask(finescaleind(fixed,i)) = 1;
        tmpalpha = energy/max(energy(:))/2;
        tmpalpha = imresize(tmpalpha,fullsiz,'method','nearest');
        alpha = alpha.*fixedmask + tmpalpha.*(1-fixedmask);% using the old intensity for parts fixed in previous scale
        alphamask = zeros(fullsiz);
        tmpmaskind = [cell2mat(region.active(:)'),cell2mat(region.fixed(:)')];
        tmpmaskind = arrayfun(@(x)finescaleind(x,i),tmpmaskind,'UniformOutput',0);
        tmpmaskind = cell2mat(tmpmaskind(:)');
        alphamask(tmpmaskind) = 1;% mask for the current regions ( active + fixed )
        figure;image(img);axis off;hold on
        h = image(tmpimg);
        set(h,'AlphaData',alpha.*alphamask);
        drawnow;hold off
    end

    function flag = dominatescale(currentind)
        upsiz = scalesiz(i:end,:);
        [m,n] = ind2sub(siz,currentind);
        currentenergy = 0;
        restenergy = 0;
        
        for jj = 1:size(upsiz,1)
            if jj == 1
                for k = 1:length(ind{i})
                    currentenergy = sum(reshape(cw{i}{ind{i}(k)}(m,n).^2,[],1))+currentenergy;
                end
            else
                m1 = floor(m/siz(1)*upsiz(jj,1));
                m2 = ceil((m-1)/siz(1)*upsiz(jj,1));
                n1 = floor(n/siz(2)*upsiz(jj,2));
                n2 = ceil((n-1)/siz(2)*upsiz(jj,2));
                m1 = max([1,min([m1,fullsiz(1)])]);
                m2 = max([1,min([m2,fullsiz(1)])]);
                n1 = max([1,min([n1,fullsiz(1)])]);
                n2 = max([1,min([n2,fullsiz(1)])]);
                
                for k = 1:length(ind{i+jj-1})
                    restenergy = sum(reshape(cw{i+jj-1}{ind{i+jj-1}(k)}(m2:m1,n2:n1).^2,[],1))/(m1-m2+1)/(n1-n2+1) + restenergy;
                end
            end
        end
        flag = currentenergy > restenergy;
    end

    function liftind = finescaleind(currentind,indsyslevel,liftlevel)
        % get all the ind contained in the currentind rectangle
        [m,n] = ind2sub(siz,currentind);
        if nargin < 3
            mlift = fullsiz(1);
            nlift = fullsiz(2);
        else
            mlift = siz(1);
            nlift = siz(2);
        end
            m0 = scalesiz(indsyslevel,1);
            n0 = scalesiz(indsyslevel,2);
        m1 = floor((m-1)*mlift/m0)+1;
        m2 = max([ceil(m*mlift/m0),m1]);
        n1 = floor((n-1)*nlift/n0)+1;
        n2 = max([ceil(n*nlift/n0),n1]);
        m1 = max([1,min([m1,mlift])]);
        m2 = max([1,min([m2,mlift])]);
        n1 = max([1,min([n1,nlift])]);
        n2 = max([1,min([n2,nlift])]);
        
        liftind = sub2ind([mlift,nlift],repmat(m1:m2,[1,n2-n1+1])',reshape(repmat(n1:n2,[m2-m1+1,1]),[],1));
        liftind = liftind(:)';
    end

    function rmemptyregion
        emptyind = cellfun(@isempty,region.active) & cellfun(@isempty,region.fixed);
        region.active = region.active(~emptyind);
        region.fixed = region.fixed(~emptyind);
        nregion = length(region.active);
    end

end


