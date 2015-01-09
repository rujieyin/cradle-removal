function [imgnew,info] = RmRegularVerticalCradleSegmentation(img,hinfo,verest,opt)
% this function removes vertical cradling by first divide the input img
% into pieces containing only a single vertical crale, and then use
% RmHorizontalCradle on the rotated subimg


% image should not have zeros background unless cradle is aligned with edge
% of the painting panel

% temporary opt input for RmHorizontalCradle
tmpopt = struct();
tmpopt.cornerind = 'single';
tmpopt.s = opt.vs;
if isfield(opt,'smoothbd')
    tmpopt.smoothbd = opt.smoothbd;
end
if ~isfield(opt,'edgecut')
    opt.edgecut = 1;
end
tmphorest = verest;

% obtain subimg indices and angle of each vcradle
[subimgind, subHangle] = get_subimage_ind(size(img),hinfo);
% initialize the output image, info
imgnew = img;
% the hcradle and vcradle cut the full image into (nh+1)-by-(nv+1) subimgs
info = cell(length(hinfo)+1,length(verest)/2);% nh = length(hinfo), nv = length(verest)/2


for i = 1:size(subimgind,1)
    % crop and get subimg
    tmpimg = img(subimgind(i,1):subimgind(i,2),:);
    
    % check if the vertical cradle has a horizontal edge cut, i.e. there is
    % a horizontal slot between the vertical cradle and the horizontal cradle
    if isfield(opt,'edgecut') && opt.edgecut
        if isfield(opt,'edgeth')
            tmpimginfo = get_edgecut_loc(tmpimg,verest,subHangle(i,:),opt.edgeth);
        else
            tmpimginfo = get_edgecut_loc(tmpimg,verest,subHangle(i,:));
        end
        % the upper and lower bound of the vertical cradle
        tmpopt.upcutind = cellfun(@(x)max(x.upcutloc(:,1)),tmpimginfo);
        tmpopt.downcutind = cellfun(@(x)min(x.downcutloc(:,1)),tmpimginfo);
        % compute the mask of left/right cuts assume the image has been
        % rotated
        tmpopt.rightmaskloc = get_maskloc(tmpimginfo,'up',size(tmpimg,1));
        tmpopt.leftmaskloc = get_maskloc(tmpimginfo,'down',size(tmpimg,1));
    end
    
    % load prefixed information, location of edge and intensity difference,
    % this option is for dealing with problematic cradles that require user
    % input 
    if isfield(opt,'VdI')
        tmpopt.HdI = opt.VdI; % intensity difference for nv vertical cradles
    end
        if isfield(opt,'Vangle')
        tmpopt.angle = opt.Vangle(i,:);
    end
    if isfield(opt,'smoothedgecut')
        tmpopt.smoothmaskbd = 1;
    end
    
    if isfield(opt,'edgeloc_v')
        tmpopt.edgeloc = opt.edgeloc_v(i,:);
    elseif isfield(opt,'edgecut') && opt.edgecut
        if isfield(opt,'smoothbd') && opt.smoothbd 
            tmpopt.edgeloc = arrayfun(@(x)tmpimginfo{x}.verest,1:length(tmpimginfo),'UniformOutput',0);
        end
    end
    
    % redo the estimation of vertical location for subimages to improve accuracy
    tmphorest = cradledetect(tmpimg,length(verest)/2,[1,size(tmpimg,2)],0,[1,size(tmpimg,1)],opt);
    tmpverest = [ 1; size(tmpimg,1)-min(tmpopt.downcutind); size(tmpimg,1)-max(tmpopt.upcutind);size(tmpimg,1)];
    
    % ==== RmHorizontalCradle ====%
    [tmpimgnew,tmpinfo] = RmHorizontalCradle(rot90(tmpimg,-1),tmphorest,tmpverest,tmpopt);
   
    if isfield(opt,'edgecut') && opt.edgecut
        for ii = 1:length(tmpinfo)
            
            tmpinfo{ii}.upcutind = tmpopt.upcutind(ii);
            tmpinfo{ii}.downcutind = tmpopt.downcutind(ii);
        end
    end
    info(i,:) = globalind(tmpinfo(:)',[subimgind(i,1),1],size(tmpimg));
    imgnew(subimgind(i,1):subimgind(i,2),:) = rot90(tmpimgnew,1);
end

%% transfer the local info of subimg into global info by adjusting the indices
    function infonew = globalind(info,initialind,sizsubimg)
        infonew = info;
        for j = 1:length(info)
            infonew{j}.colind = info{j}.rowind + initialind(2) - 1;
            infonew{j}.rowind = info{j}.colind + initialind(1) - 1;
            infonew{j}.lp = [initialind(1),info{j}.rp + initialind(2) - 1];
            infonew{j}.rp = [initialind(1),info{j}.rd + initialind(2) - 1];
            infonew{j}.ld = [initialind(1) + sizsubimg(1) - 1,info{j}.lp + initialind(2) - 1];
            infonew{j}.rd = [initialind(1) + sizsubimg(1) - 1,info{j}.ld + initialind(2) - 1];
            infonew{j}.cradleimg = rot90(info{j}.cradleimg,1);
            if isfield(info{j},'upcutind')
                infonew{j}.upcutind = info{j}.upcutind + infonew{j}.rowind(1) - 1;
                infonew{j}.downcutind = info{j}.downcutind + infonew{j}.rowind(1) - 1;
            end
            try
                infonew{j}.edgeloc = info{j}.edgeloc;
            end
        end
    end

%% compute the location of mask of the slot between hcradle and the cut of vcradle
    function maskloc = get_maskloc(imginfo,dir,l)
        
        n = length(imginfo);
        switch dir
            case 'up'
                if isfield(imginfo{1},'upcutloc')
                    maskloc = arrayfun(@(x)[l+1-imginfo{x}.upcutloc(:,1) imginfo{x}.upcutloc(:,2)],1:n,'UniformOutput',0);
                else
                    maskloc = [];
                end
            case 'down'
                if isfield(imginfo{1},'downcutloc')
                    maskloc = arrayfun(@(x)[l+1-imginfo{x}.downcutloc(:,1) imginfo{x}.downcutloc(:,2)], 1:n, 'UniformOutput',0);
                else
                    maskloc = [];
                end
            otherwise
                maskloc = [];
        end
    end
    function [subind, subangle] = get_subimage_ind(imgsize,hinfo)
        nh = length(hinfo);
        subind = zeros(nh+1,2);
        subangle = zeros(nh+1,2);
        subind(1,1) = 1;
        subangle(1,1) = -1;
        subind(end,2) = imgsize(1);
        subangle(end,2) = -1;
        for k = 1:nh
            subind(k,2) = max(hinfo{k}.lp(1),hinfo{k}.rp(1));
            subind(k+1,1) = min(hinfo{k}.ld(1),hinfo{k}.rd(1));
            subangle(k,2) = hinfo{k}.angle;
            subangle(k+1,1) = hinfo{k}.angle;
        end
    end

%% detect the location of edge cut
    function info = get_edgecut_loc(tmpimg,verest,hangle,th)
        tmpimg = tmpimg.*(tmpimg > 0);
        nv = length(verest)/2;
        info = cell(nv,1);
        % set threshold for cut detection
        if nargin < 4
            th = 20;
        end
        th = th*3;
        % compute the gradient of image
        gradmap1 = imfilter(tmpimg(1:50,:),[-ones(5,1);ones(5,1)],'replicate');
        gradmap1 = gradmap1.*(gradmap1>th);
        gradmap2 = imfilter(tmpimg(end-49:end,:),[-ones(5,1);ones(5,1)],'replicate');
        gradmap2 = -gradmap2.*(gradmap2 < -th);
        for j = 1:nv
            % check if there is cut above
            tmpest = verest((2*j-1) : 2*j);
            info{j}.ind = [max(1,tmpest(1) - 30) min(size(tmpimg,2), tmpest(2) + 30)];
            tmpind = info{j}.ind;
            tmpgrad = gradmap1(:,verest(j*2-1):verest(j*2))>th;
            if hangle(1) > 0 && (sum(max(tmpgrad,[],1)) > size(tmpgrad,2)*.7)
                info{j}.upcut = 1;
                tmpgrad = gradmap1(:,tmpind(1):tmpind(2));
                [newverest,cutloc] = get_new_verest(tmpgrad,hangle(1),th);
                info{j}.upverest = newverest + tmpind(1)-1;
                info{j}.upcutloc = [cutloc(:) (tmpind(1):tmpind(2))'];
            else
                info{j}.upcut = 0;
                info{j}.upverest = [];
                info{j}.upcutloc = [ones(tmpind(2)-tmpind(1)+1,1) (tmpind(1):tmpind(2))'];
                
            end
            % check if there is cut below
            tmpgrad = gradmap2(:,verest(j*2-1):verest(j*2)) > th;
            tmpest = verest((2*j-1) : 2*j);
            if hangle(2) > 0 && (sum(max(tmpgrad,[],1)) > size(tmpgrad,2)*.7)
                info{j}.downcut = 1;
                tmpgrad = gradmap2(:,tmpind(1):tmpind(2));
                [newverest,cutloc] = get_new_verest(tmpgrad,hangle(2),th);
                info{j}.downverest = newverest + tmpind(1)-1;
                info{j}.downcutloc = [cutloc(:)+size(tmpimg,1)-50 (tmpind(1):tmpind(2))'];
            else
                info{j}.downcut = 0;
                info{j}.downverest = [];
                info{j}.downcutloc = [size(tmpimg,1)*ones(tmpind(2)-tmpind(1)+1,1) (tmpind(1):tmpind(2))'];
            end
            % renew estimate of verest if cuts exit
            if ~isfield(info{j},'verest') || isempty(info{j}.verest)
                info{j}.verest = verest(2*j-1:2*j);
                info{j}.verest = info{j}.verest(:)';
            end
            info{j}.verest = round(mean(info{j}.verest,1));
        end
    end

%%  refine estimation of vertical cradle location
    function [verest,cut] = get_new_verest(gradmap,angle,th)
        [R,xp] = radon(gradmap,angle);
        [~,binind] = max(R);
        tmpedge = PointBackProject(xp(binind),angle,size(gradmap));
        cut = tmpedge;
        tmpedge = sub2ind(size(gradmap),tmpedge,1:size(gradmap,2));
        tmpl = find(gradmap(tmpedge) > th,1,'first');
        tmpr = find(gradmap(tmpedge) > th,1,'last');
        [tmplr,tmplc] = ind2sub(size(gradmap),tmpedge(tmpl));
        [tmprr,tmprc] = ind2sub(size(gradmap),tmpedge(tmpr));
        verest = [tmplc tmprc];
    end
end