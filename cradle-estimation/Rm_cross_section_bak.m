function [imgnew, crossinfo] = Rm_cross_section(img,hinfo,vinfo,opt)
% this function remove the cross section of horizontal and vertical cradles

if nargin < 4
    opt = struct();
    opt.smoothbd = 0;
end

if size(vinfo,1) > 1
    vinfo = connect_vinfo(vinfo);
end

nv = length(vinfo);
nh = length(hinfo);
crossinfo = cell(nh,nv);
imgnew = img;

for ih = 1 : nh
    for iv = 1 : nv
        info = struct();
        hrowind = hinfo{ih}.rowind;
        hcolind = hinfo{ih}.colind;
        vrowind = vinfo{iv}.rowind;
        vcolind = vinfo{iv}.colind;
        info.rowind = [max(hrowind(1),vrowind(1)) min(hrowind(2),vrowind(2))];
        info.colind = [max(hcolind(1),vcolind(1)) min(hcolind(2),vcolind(2))];
        info.crossimg = get_crossimg(hinfo{ih},vinfo{iv},info);
        subimg = img(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2));
        if isfield(vinfo{iv},'upcutind')
            d1 = min(abs(vinfo{iv}.downcutind - max(hinfo{ih}.lp(1),hinfo{ih}.rp(1))));%edge cut width
            d2 = min(abs(vinfo{iv}.upcutind - min(hinfo{ih}.ld(1),hinfo{ih}.rd(1))));
            d = max([d1,d2,30]);
            se2 = strel('disk',2*d);
        else
            se2 = strel('disk',60);
        end
        se = strel('disk',30);
        bdin = info.crossimg - imerode(double(info.crossimg),se);
        bdin = bdin > max(bdin(:))*.3;
        bdout = logical(imdilate(double(info.crossimg),se2) - imdilate(double(info.crossimg),se));
        bdout = bdout > max(bdout(:))*.3;
        if isfield(opt,'crossIntense')
            if isscalar(opt.crossIntense)
            info.deltaIntense = opt.crossIntense;
            else
                info.deltaIntense = opt.crossIntense(ih,iv);
            end
        else
            info.deltaIntense = mean(subimg(bdin)) - mean(subimg(bdout));
        end
        imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimg - info.crossimg*info.deltaIntense;
        crossinfo{ih,iv} = info;
        try
            crossinfo{ih,iv}.cutwidth = [d1 d2];
        catch
            crossinfo{ih,iv}.cutwidth = [0,0];
        end
    end
end

    function crossimg = get_crossimg(hinfo,vinfo,info)
        rind = info.rowind;
        rind1 = hinfo.rowind;
        rind2 = vinfo.rowind;
        cind1 = hinfo.colind;
        cind2 = vinfo.colind;
        cind = info.colind;
        imgsiz = [info.rowind(2)-info.rowind(1)+1 info.colind(2)-info.colind(1)+1];
        crossimg = hinfo.cradleimg((rind(1) - rind1(1) + 1):(rind(2) - rind1(1) + 1),(cind(1) - cind1(1) + 1):(cind(2) - cind1(1) + 1));
        crossimg = crossimg.*vinfo.cradleimg((rind(1) - rind2(1) + 1):(rind(2) - rind2(1) + 1),(cind(1) - cind2(1) + 1):(cind(2) - cind2(1) + 1));
        crossimg = crossimg/max(crossimg(:));
    end

end