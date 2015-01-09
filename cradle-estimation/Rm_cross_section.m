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

if strcmp(model,'multiplicative')
    if isfield(opt,'crossmodel')
        crossmodel = opt.crossmodel;
    else
        crossmodel = 'profilefit';
    end
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
        
        %*****************check if vinfo is connected or not!*************%
        
        
        % get the indices of subimg containing the cross-section
        if isfield(vinfo{1,iv},'upcutind')
            if size(vinfo,1) == 1
                % use the indices of edge_cut
                info.rowind = [vinfo{iv}.downcutind(ih)-vinfo{iv}.width(ih)...
                    vinfo{iv}.upcutind(ih+1)+vinfo{iv}.width(ih+1)];
                info.cutind = [vinfo{iv}.width(ih) vinfo{iv}.upcutind(ih+1)-info.rowind(1)];
            else
                info.rowind = [vinfo{ih,iv}.downcutind-vinfo{ih,iv}.width...
                    vinfo{ih+1,iv}.upcutind+vinfo{ih+1,iv}.width];
                info.cutind = [vinfo{ih,iv}.width vinfo{ih+1,iv}.upcutind - info.rowind(1)];
            end
        else
            hrowind = hinfo{ih}.rowind;
            info.rowind = [hrowind(1)-opt.hs hrowind(2)+opt.hs];
        end
        if size(vinfo,1) > 1
            info.colind = vinfo{ih,iv}.colind;
        else
            info.colind = vinfo{iv}.colind;
        end
        
        subimg = img(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2));
        
        switch model
            case 'additive'
                
                % doesn't touch the slot b/t hcradle and vcradle
                % find the shape of the crossimg, with maximum intensity 1
                info.crossimg = get_crossimg(hinfo{ih},vinfo{iv},info);
                
            case 'multiplicative'
                % find the mask used in vcradle removal
                mask = crop_img(info.rowind,info.colind,vinfo{ih,iv}.rowind,...
                    vinfo{ih,iv}.colind,vinfo{ih,iv}.cradleimg);
                mask = mask + crop_img(info.rowind,info.colind,vinfo{ih+1,iv}.rowind,...
                    vinfo{ih+1,iv}.colind,vinfo{ih+1,iv}.cradleimg);
                mask = mask == 0;
                % calculate the input of RemoveAttenuationProfile
                % use the profile of vcradle to deal with low cross
                
                % find midpoint of colind use interpolation, the most
                % robust way !
                m1 = mean(vinfo{ih,iv}.rowind); m2 = mean(vinfo{ih+1,iv}.rowind);
                m = mean(info.rowind);
                mp1 = vinfo{ih,iv}.midpoint + vinfo{ih,iv}.colind(1) - 1;
                mp2 = vinfo{ih+1,iv}.midpoint + vinfo{ih+1,iv}.colind(1) - 1;
                midpoint = (mp1*(m2-m) + mp2*(m-m1))/(m2-m1)-info.colind(1)+1;
                midpoint = round(midpoint);
                angle = (vinfo{ih,iv}.angle + vinfo{ih+1,iv}.angle)/2;
                
                % back up info
                infotmp = info;
                switch crossmodel
                    case 'profilefit'
                        angle = angle/180*pi;
                        % compute transformed image using vcradle profile
                        [subimgnew,info] = RemoveAttenuationProfile(rot90(subimg,-1),midpoint,angle,...
                            vinfo{ih,iv}.width,vinfo{ih,iv}.p1profile,vinfo{ih,iv}.C);
                    case 'linearfit'
                        tmpverest = [1; size(subimg,1)+[-info.cutind(2);-info.cutind(1);0]];
                        % remove cradle by fitting attenuation model directly
                        [subimgnew,info] = cradle_attenuation_fitting(rot90(subimg,-1),midpoint,angle,tmpverest);
                end
                subimgnew = rot90(subimgnew,1);
                % restore info
                info.rowind = infotmp.rowind;
                info.colind = infotmp.colind;
                info.cutind = infotmp.cutind;
                % compute and normalize difference img
                info.crossimg = (subimg-subimgnew).*mask;
                
                if strcmp(crossmodel,'linearfit')
                   imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimgnew.*mask + subimg.*(~mask);
                   crossinfo{ih,iv} = info;
                   continue;
                end
                
                %  % code for scaling the difference image
                %                 info.crossimg = info.crossimg/median( info.crossimg(info.crossimg(:)>0) );
                
                %                 % subtract cradleimage with mask
                %                 imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimg.*(~mask) + subimgnew.*mask;
                %                 crossinfo{ih,iv} = info;
        end
        
        % estimate the boundaries of cross-section
        interpolate = @(x,y)round((x*(nv-iv) + y*iv)/nv);
        
        bdrowind(1) = interpolate(hinfo{ih}.lp(1), hinfo{ih}.rp(1));
        bdrowind(2) = interpolate(hinfo{ih}.ld(1), hinfo{ih}.rd(1));
        bdrowind = bdrowind - info.rowind(1) + 1;
        
        vbd = any(info.crossimg > 0);
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
            avrgin = median(subimg(bdin));
            avrgout = median(subimg(bdout));
            info.deltaIntense = avrgin - avrgout;
        end
        
        switch model
            case 'multiplicative'
                C = vinfo{ih,iv}.C;
                x = avrgin;
                y2 = avrgout;
                y1 = x - median(info.crossimg(bdin));
                width = vinfo{ih,iv}.width;
                p1 = median(vinfo{ih,iv}.p1profile([width+1:2*width end-2*width:end-width]));
                a = 1 - (y1-y2)/(x-C)/p1; % a should be b/t [0,1]
                subimgnew = RemoveAttenuationProfile(rot90(subimg,-1),midpoint,angle,...
                    vinfo{ih,iv}.width,vinfo{ih,iv}.p1profile*a,vinfo{ih,iv}.C);
                subimgnew = rot90(subimgnew,1);
                % discard possible increas of pixel value, due to negative
                % difference around boundary of vcradles
                mask = (subimgnew < subimg) & mask;
                imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimgnew.*mask + subimg.*(~mask);
                info.a = a;
            case 'additive'
                % normalize the crossimg
                info.crossimg = info.crossimg/median(info.crossimg(bdin));
                
                % subtract shape (info.crossimg) times intensity difference
                imgnew(info.rowind(1):info.rowind(2),info.colind(1):info.colind(2)) = subimg - info.crossimg*info.deltaIntense;
                
        end
        
        crossinfo{ih,iv} = info;
        
    end
end
%%
    function cropmask = crop_img( r, c, mr, mc, mask)
        % this function crop the mask with row indices mr1:mr2, column
        % indices mc1:mc2 that is within the region r1:r2, c1:c2
        r1 = r(1); r2 = r(2); c1 = c(1); c2 = c(2);
        mr1 = mr(1); mr2 = mr(2); mc1 = mc(1); mc2 = mc(2);
        cropmask = zeros(r2-r1+1,c2-c1+1);
        C1 = max(c1,mc1); C2 = min(c2,mc2);
        R1 = max(r1,mr1); R2 = min(r2,mr2);
        cropmask( (R1:R2)-r1+1, (C1:C2)-c1+1 ) = mask( (R1:R2)-mr1+1, (C1:C2)-mc1+1 );
    end


%%
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
        %         crossimg = crossimg/max(crossimg(:));
    end

end