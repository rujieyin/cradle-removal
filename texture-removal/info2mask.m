function mask = info2mask(info,imgind)
% imgind is a struct or matrix
% imgind.rowind imgind.colind
% or [rowind;colind]

if isstruct(imgind)
    rowind = imgind.rowind;
    colind = imgind.colind;
else
    rowind = imgind(1,:);
    colind = imgind(2,:);
end

imgsize = [rowind(2) - rowind(1) + 1, colind(2) - colind(1) + 1];
mask = zeros(imgsize);

for i = 1:size(info,1)
    for j = 1:size(info,2)
        tmpinfo = info{i,j};
        tmpinfo.cradleimg = fillmask(tmpinfo.cradleimg,tmpinfo);
        if overlap_interval(rowind,tmpinfo.rowind) && overlap_interval(colind,tmpinfo.colind)
            tmprowind = intersect(rowind(1):rowind(2),tmpinfo.rowind(1):tmpinfo.rowind(2));
            tmpcolind = intersect(colind(1):colind(2),tmpinfo.colind(1):tmpinfo.colind(2));
            tmprowind1 = [min(tmprowind) max(tmprowind)] - rowind(1) + 1;
            tmpcolind1 = [min(tmpcolind) max(tmpcolind)] - colind(1) + 1;
            tmprowind2 = [min(tmprowind) max(tmprowind)] - tmpinfo.rowind(1) + 1;
            tmpcolind2 = [min(tmpcolind) max(tmpcolind)] - tmpinfo.colind(1) + 1;
            subcradleimg = tmpinfo.cradleimg(tmprowind2(1):tmprowind2(2),tmpcolind2(1):tmpcolind2(2));
            mask(tmprowind1(1):tmprowind1(2),tmpcolind1(1):tmpcolind1(2)) = subcradleimg > (max(subcradleimg(:))/2);
        end
    end
end

    function mask = fillmask(cradleimg,tmpinfo)
        rind = [max(tmpinfo.lp(1),tmpinfo.rp(1)) min(tmpinfo.ld(1),tmpinfo.rd(1))] - tmpinfo.rowind(1) + 1;
        cind = [max(tmpinfo.lp(2),tmpinfo.ld(2)) min(tmpinfo.rp(2),tmpinfo.rd(2))] - tmpinfo.colind(1) + 1;
        mask = cradleimg > (max(cradleimg(:))/2);
        mask(rind(1):rind(2),cind(1):cind(2)) = 1;
    end

    function flag = overlap_interval(seg1,seg2)
        flag = ~( (max(seg1) < min(seg2)) | (max(seg2) < min(seg1)) );
    end

end