function [imgnew,info] = RmHorizontalCradle(img,horest,verest,opt)
% this function try to find the horizontal cradling edges


if ~isfield(opt,'cornerind')
    opt.cornerind = 'full'; % whether ind of corners in both cords are saved
end

vmask = ones(1,size(img,2));
for i = 1:length(verest)/2
    vmask(verest(2*i-1):verest(2*i)) = 0;
end

nh = length(horest)/2;
info = cell(nh,1);
imgnew = img;

for i = 1 : nh
    if isfield(opt,'ind')
        initialind = opt.ind{i}(1);
        endind = opt.ind{i}(2);
        imgtmp = img(opt.ind{i}(1):opt.ind{i}(2),:);
    else
        if i == 1
            initialind = max(1,horest(2*i-1) - opt.s*4);
        else
            initialind = horest(2*i-2) + opt.s;
        end
        if i == nh
            endind = min(size(img,1),horest(2*i) + opt.s*4);
        else
            endind = horest(2*i+1) - opt.s;
        end
        imgtmp = img(initialind:endind,:);
    end
    mask = ones(size(imgtmp));
    if isfield(opt,'leftmaskloc')
        for j = 1:size(opt.leftmaskloc{i},1)
            mask(opt.leftmaskloc{i}(j,2) - initialind + 1,1:opt.leftmaskloc{i}(j,1)) = 0;
        end
    end
    if isfield(opt,'rightmaskloc')
        for j = 1:size(opt.rightmaskloc{i},1)
            mask(opt.rightmaskloc{i}(j,2) - initialind + 1,opt.rightmaskloc{i}(j):end) = 0;
        end
    end
    if isfield(opt,'edgeloc') || isfield(opt, 'angle') || isfield(opt,'HdI')
        s = 20;
        prefix = struct();
        try
            prefix.edgeloc = opt.edgeloc{i} - initialind + 1;
        end
        try
            prefix.angle = opt.angle(i);
        end
        try
            if length(opt.HdI) == 1
                prefix.dI = opt.HdI;
            else
                prefix.dI = opt.HdI(i);
            end
        end
        [imgnewtmp,infotmp] = rm_single_cradle(imgtmp.*mask,horest(2*i-1:2*i)-initialind + 1,vmask,s,prefix);
    else
        s = opt.s*2;
        [imgnewtmp,infotmp] = rm_single_cradle(imgtmp.*mask,horest(2*i-1:2*i)-initialind + 1,vmask,s);
    end
    %     figure;subplot(2,1,1);imshow(imgtmp,[]);subplot(2,1,2);imshow(imgnewtmp,[]);
    imgnew(initialind:endind,:) = imgnewtmp.*mask + imgtmp.*(1-mask);
    infotmp.rowind = [initialind endind];
    infotmp.colind = [1 size(imgnewtmp,2)];
    infotmp.cradleimg = infotmp.cradleimg.*mask;
    info{i} = globalind(infotmp,initialind,opt.cornerind);
    try
        info{i}.edgeloc = opt.edgeloc{i};
    end
end

    function infonew = globalind(info,initialind,flag)
        infonew = info;
        if strcmp(flag,'full')
            infonew.lp = [info.lp + initialind - 1, info.colind(1)];
            infonew.ld = [info.ld + initialind - 1, info.colind(1)];
            infonew.rp = [info.rp + initialind - 1, info.colind(2)];
            infonew.rd = [info.rd + initialind - 1, info.colind(2)];
        else
            infonew.lp = info.lp + initialind - 1;
            infonew.ld = info.ld + initialind - 1;
            infonew.rp = info.rp + initialind - 1;
            infonew.rd = info.rd + initialind - 1;
        end
    end




    function [imgnew,info] = rm_single_cradle(img,loc,vmask,s,prefix)
        
        if nargin > 4 && isfield(prefix,'edgeloc')
            loc = prefix.edgeloc;
        end
        s = min([(loc(1) - 1)/2,(size(img,1) - loc(2))/2, s]);
        s = round(s);
        if nargin > 4 && isfield(prefix,'angle')
            theta = prefix.angle;
        else
            theta = 85:.1:95;
        end
        
        ind1 = [loc(1) - s,loc(1) + s];
        gradmap = imfilter(img(ind1(1):ind1(2),:),[-ones(5,1);ones(5,1)],'symmetric');
        gradmap = gradmap.*(gradmap>0);
        [R,xp] = radon(gradmap,theta);
        [angle1,R] = get_angle(R,theta);
        if nargin > 4 && isfield(prefix,'edgeloc')
            ind1 = prefix.edgeloc(1);
        else
            [~,peak] = max(R);
            ind1 = loc(1) - xp(peak);
        end
        
        ind2 = [loc(2) - s,loc(2)+s];
        gradmap = imfilter(img(ind2(1):ind2(2),:),[-ones(5,1);ones(5,1)],'symmetric');
        gradmap = -gradmap.*(gradmap < 0);
        [R,xp] = radon(gradmap,theta);
        [angle2,R] = get_angle(R,theta);
        if nargin > 4 && isfield(prefix,'edgeloc')
            ind2 = prefix.edgeloc(2);
        else
            [~,peak] = max(R);
            ind2 = loc(2) - xp(peak);
        end
        
        angle = (angle1+angle2)/2;
        
        subimg = img(ind1-s:ind1+s,:);
        subimg = bsxfun(@times,subimg,vmask);
        %     figure;imshow(subimg,[]);
        [R1,xp1] = radon(subimg,angle);
        %     figure;plot(1:length(R1),R1);
        [edgeshape1,dI1] = get_edges(R1,xp1,angle,'first',s);
        %     figure;plot(1:(2*s-1),R(abs(xp) < s)' - edgeshape*dI);
        subimg = img(ind2-s:ind2+s,:);
        subimg = bsxfun(@times,subimg,vmask);
        [R2,xp2] = radon(subimg,angle);
        %         figure;plot(1:length(R2),R2);
        [edgeshape2,dI2] = get_edges(R2,xp2,angle,'second',s);
        dI = (dI1 + dI2)/2;
        deltaIntense = size(img,2)/sum(vmask);
        %         figure;plot(1:length(R1),R1 - edgeshape1*dI1,1:length(R1),R1);
        %         figure;plot(1:length(R2),R2 - edgeshape2*dI2,1:length(R2),R2);
        bound1 = BackProjection(edgeshape1*deltaIntense/dI1*dI,xp1,angle,[size(subimg,1)*3,size(subimg,2)]);
        % enlarge the size of back projection image to exclude artifacts
        bound1 = bound1((size(subimg,1)+1):size(subimg,1)*2,:);
        bound2 = BackProjection(edgeshape2*deltaIntense/dI2*dI,xp2,angle,[size(subimg,1)*3,size(subimg,2)]);
        bound2 = bound2((size(subimg,1)+1):size(subimg,1)*2,:);
        cradleimg = zeros(size(img));
        cradleimg(ind1-s:ind1+s,:) = bound1;
        cradleimg(ind2-s:ind2+s,:) = bound2;
        tmpimg = arrayfun(@(x)full_bar(cradleimg(ind1-s:ind2+s,x)),1:size(img,2),'UniformOutput',0);
        tmpimg = cell2mat(tmpimg(:)');
        th = max(cradleimg(:))/3;
        info = struct();
        info.lp = find(cradleimg(:,1)>th,1,'first');
        info.ld = find(cradleimg(:,1)>th,1,'last');
        info.rp = find(cradleimg(:,end) > th, 1,'first');
        info.rd = find(cradleimg(:,end) > th, 1, 'last');
        if nargin > 4 && isfield(prefix,'dI')
            tmpimg = tmpimg/max(tmpimg(:))*prefix.dI;
        end
        cradleimg(ind1 - s:ind2+s,:) = tmpimg;
        
%         cradleimg(ind1-s:ind1+s,:) = apply_cradlemask(img(ind1-s:ind1+s,:),tmpimg(1:2*s+1,:));
%         cradleimg(ind2-s:ind2+s,:) = apply_cradlemask(img(ind2-s:ind2+s,:),tmpimg(end-2*s:end,:));
        
%         cradleimg(ind1-s:ind2+s,:) = apply_cradlemask(img(ind1-s:ind2+s,:),tmpimg);
        
        %         figure;imshow(cradleimg,[]);
        imgnew = img-cradleimg;
        figure;imshow([imgnew; img],[]);
        
        info.cradleimg = cradleimg;
        info.angle = angle;
        info.dI = max(cradleimg(:)) - min(cradleimg(:));
    end

    function signal = full_bar(signal)
        L = round(length(signal)/2);
        [pk1,ind1] = max(signal(1:L));
        [pk2,ind2] = max(signal(end-L:end));
        ind2 = ind2 + L;
        signal(ind1:ind2) = linspace(pk1,pk2,ind2-ind1+1);
    end

    function [angle,R] = get_angle(R,theta)
        [~,engind] = max(sum(R.^2,1));
        angle = theta(engind);
        R = R(:,engind);
    end

    function [edgeshape,dI] = get_edges(R,xp,angle,bd,s)
        if strcmp(bd,'second')
            R = flipud(R(:));
        end
        [n,x] = hist(R(abs(xp) < sin(angle/180*pi)*s),round(s/3));
        binsize = x(2)-x(1);
        %         figure;plot(x,n);
        [~,pks] = findpeaks([0 n 0]);
        highintense = x(max(pks)-1);
        lowintense = x(min(pks)-1);
        dI = highintense - lowintense;
        if isfield(opt,'smoothbd') && opt.smoothbd
            usecdf = 1;
            leftshape = get_edgeshape(R(xp > -s & xp <= 0),'l',highintense,lowintense,usecdf);
            rightshape = get_edgeshape(R(xp >= 0 & xp < s),'r',highintense,lowintense,usecdf);
        else
            usecdf = 0;
            leftshape = get_edgeshape(R(xp > -s & xp <= 0),'l',highintense-binsize,lowintense+binsize,usecdf);
            rightshape = get_edgeshape(R(xp >= 0 & xp < s),'r',highintense-binsize,lowintense+binsize,usecdf);
        end
        edgeshape = [leftshape rightshape(2:end)];
        edgeshape = edgeshape - edgeshape(end);
        tmpind = find(edgeshape > 0 ,1,'first');
        edgeshape(1:tmpind-1) = 0;
        tmpind = find(edgeshape < 0,1,'first');
        edgeshape(tmpind:end) = 0;
        dI = max(edgeshape) - min(edgeshape);
        if strcmp(bd,'second')
            edgeshape = fliplr(edgeshape);
        end
        L = length(R)- length(edgeshape);
        edgeshape = [zeros(floor(L/2),1);edgeshape(:);zeros(ceil(L/2),1)];
    end

    function [edgeshape] = get_edgeshape(r,dir,highintense,lowintense,usecdf)
        if strcmp(dir,'l')
            ind1 = find(r > highintense,1,'last');
            ind2 = find(r(1:end-1) - r(2:end) < 0,1,'last');
            ind = min([ind1,ind2]);
            if r(ind) < r(end)
                ind = find(r > r(end),1,'last');
            end
            if ind >= length(r)
                ind = length(r) - 1;
            end
            edgeshape = zeros(1,length(r));
            if nargin > 4 && usecdf
                L = length(r) - ind;
                edgeshape(ind:end) = 1 -cdf('norm',1:(L+1),L+1,L/3);
                edgeshape = edgeshape*(r(ind) - r(end))*2 + 2*r(end) - r(ind);
            else
                edgeshape(ind:end) = r(ind:end);
            end
            %             edgeshape(1:ind-1) = r(ind);
        else
            ind1 = find(r < lowintense,1,'first');
            ind2 = find(r(1:end-1) - r(2:end) < 0,1,'first');
            ind = max([ind1,ind2])-1;
            if ind <= 1
                ind = 2;
            end
            %             edgeshape = edgeshape.*(edgeshape > 0);
            %             tmpind = find(edgeshape == 0, 1 ,'first');
            %             edgeshape(tmpind:end) = 0;
            edgeshape = zeros(1,length(r));
            if nargin > 4 && usecdf
                L = ind;
                edgeshape(1:ind) = 1 -cdf('norm',1:L,1,L/3);
                edgeshape = edgeshape*(r(1)-r(ind))*2 + r(ind);
            else
                edgeshape(1:ind) = r(1:ind);
                edgeshape(ind+1:end) = r(ind);
            end
            %             edgeshape(ind+1:end) = r(ind);
        end
    end

    function cradlenew = apply_cradlemask(img,cradle)
        cradlemask = img > (max(cradle(:)) + min(img(img>0)));
        cradlenew = cradle.*cradlemask;
    end
end
