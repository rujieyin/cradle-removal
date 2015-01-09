function [imgnew, crossinfo] = Rm_cross_section(img,hinfo,vinfo,opt)
% this function remove the cross section of horizontal and vertical cradles

if ~isfield(opt,'smoothbd')
    opt.smoothbd = 0;
end

if isfield(opt,'model')
    model = opt.model;
else
    model = 'additive';
end

if strcmp(model,'additive') && size(vinfo,1) > 1
    vinfo = connect_vinfo(vinfo);
end

nv = length(vinfo);
nh = length(hinfo);
crossinfo = cell(nh,nv);
imgnew = img;

for ih = 1 : nh
    for iv = 1 : nv
        info = struct();
        % get the indices of subimg containing the cross-section
        if isfield(vinfo{iv},'upcutind')
            % use the indices of edge_cut
            info.rowind = [vinfo{iv}.downcutind(ih) - vinfo{iv}.width(ih)...
                vinfo{iv}.upcutind(ih+1) + vinfo{iv}.width(ih+1)];
            info.cutind = [vinfo{iv}.width(ih) vinfo{iv}.upcutind(ih+1)-info.rowind(1)];
        else            
            hrowind = hinfo{ih}.rowind;
            info.rowind = [hrowind(1)-opt.hs hrowind(2)+opt.hs];
        end
        info.colind = vinfo{iv}.colind;
        
        subimg = img(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2));
        
        switch model
            case 'additive'
                % doesn't touch the slot b/t hcradle and vcradle
                % find the shape of the crossimg, with maximum intensity 1
                info.crossimg = get_crossimg(hinfo{ih},vinfo{iv},info);
                
                % estimate the boundaries of cross-section
                interpolate = @(x,y)round((x*(nv-iv) + y*iv)/nv);
                
                bdrowind(1) = interpolate(hinfo{ih}.lp(1), hinfo{ih}.rp(1));
                bdrowind(2) = interpolate(hinfo{ih}.ld(1), hinfo{ih}.rd(1));
                bdrowind = bdrowind - info.rowind(1) + 1;
                
                vbd = any(info.crossimg);
                bdcolind(1) = find(vbd,1);
                bdcolind(2) = find(vbd,1,'last');
                
                % obtain the rectangular ring inside the boundary
                bdin = zeros(size(subimg));
                bdin([bdrowind(1):bdrowind(1)+opt.hs bdrowind(2)-opt.hs:bdrowind(2)], bdcolind(1):bdcolind(2)) = 1;
                bdin(bdrowind(1):bdrowind(2), [bdcolind(1):bdcolind(1)+opt.vs bdcolind(2)-opt.vs:bdcolind(2)]) = 1;
                bdin = logical(bdin);
                
                % obtain the rectangular ring outside the boundary, if
                % there is a slot, then the slot region is discarded
                bdout = zeros(size(subimg));
                vs = min([opt.vs, bdcolind(1)+1, size(subimg,2)-bdcolind(2)]);
                if isfield(info,'cutind')
                    bdout([1:info.cutind(1) info.cutind(2):end],bdcolind(1)-vs:bdcolind(2)+vs) = 1;
                else
                    bdout([1:bdrowind(1) bdrowind(2):end],bdcolind(1)-vs:bdcolind(2)+vs) = 1;
                end                    
                bdout(1:end, [bdcolind(1)-vs:bdcolind(1) bdcolind(2):bdcolind(2)+vs]) = 1;
                bdout = logical(bdout);
               
                % estimate the intensity difference in/out
                if isfield(opt,'crossIntense')
                    if isscalar(opt.crossIntense)
                        info.deltaIntense = opt.crossIntense;
                    else
                        info.deltaIntense = opt.crossIntense(ih,iv);
                    end
                else
                    info.deltaIntense = median(subimg(bdin)) - median(subimg(bdout));
                end
                % subtract shape (info.crossimg) times intensity difference
                imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimg - info.crossimg*info.deltaIntense;
                crossinfo{ih,iv} = info;
            case 'multiplicative'
        end
    end
end

    function crossimg = get_crossimg(hinfo,vinfo,info)
        % crossimg is the product of vcradleimg and hcradleimg, then
        % normalized
        rind = info.rowind;
        rind1 = hinfo.rowind;
        rind2 = vinfo.rowind;
        cind1 = hinfo.colind;
        cind2 = vinfo.colind;
        cind = info.colind;
        if  rind(1) < rind1(1)
            crossimg = zeros(rind(2)-rind(1)+1,cind(2)-cind(1)+1);
            crossimg(rind1(1)-rind(1)+1:rind1(2)-rind(1)+1,:) = hinfo.cradleimg(:,(cind(1) - cind1(1) + 1):(cind(2) - cind1(1) + 1));
        else
            crossimg = hinfo.cradleimg((rind(1) - rind1(1) + 1):(rind(2) - rind1(1) + 1),(cind(1) - cind1(1) + 1):(cind(2) - cind1(1) + 1));
        end
        crossimg = crossimg.*vinfo.cradleimg((rind(1) - rind2(1) + 1):(rind(2) - rind2(1) + 1),(cind(1) - cind2(1) + 1):(cind(2) - cind2(1) + 1));
        crossimg = crossimg/max(crossimg(:));
    end

end