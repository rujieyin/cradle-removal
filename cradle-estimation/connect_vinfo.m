function vinfonew = connect_vinfo(vinfo)
% this function connect vertical segments into a full piece

vinfonew = cell(1,size(vinfo,2));
nv = size(vinfo,2);
for j = 1:nv
    info = struct();
    info.lp = vinfo{1,j}.lp;
    info.ld = vinfo{end,j}.ld;
    info.rp = vinfo{1,j}.rp;
    info.rd = vinfo{end,j}.rd;
    info.rowind = [vinfo{1,j}.rowind(1) vinfo{end,j}.rowind(2)];
    info.colind = [vinfo{1,j}.colind(1) vinfo{end,j}.colind(2)];
    try
    info.upcutind = cellfun(@(x)x.upcutind,vinfo(:,j));
    info.downcutind = cellfun(@(x)x.downcutind,vinfo(:,j));
    end
    info.cradleimg = zeros(info.rowind(2) - info.rowind(1)+1,info.colind(2)-info.colind(1)+1);
%     left = linspace(info.lp(2),info.ld(2),info.ld(1)-info.lp(1)+1) - info.colind(1) + 1;
%     left = round(left);
%     right = linspace(info.rp(2),info.rd(2),info.rd(1) - info.rp(1)+1) - info.colind(1) + 1;
%     right = round(right);
%     up = linspace(info.lp(1),info.rp(1),abs(info.lp(2) - info.rp(2))+1) - info.rowind(1) + 1;
%     up = round(up);
%     down = linspace(info.ld(1),info.rd(1),abs(info.ld(2) - info.rd(2))+1) - info.rowind(1) + 1;
%     down = round(down);
%     mask1 = zeros(size(info.cradleimg));
%     ind = (info.lp(1):info.ld(1)) - info.rowind(1) + 1;
%     for i = 1:length(ind)
%         mask1(ind(i),left(i):end) = 1;
%     end
%     mask2 = zeros(size(info.cradleimg));
%     ind = (info.rp(1) : info.rd(1)) - info.rowind(1) + 1;
%     for i = 1:length(ind)
%         mask2(ind(i), 1:right(i)) = 1;
%     end
%     mask3 = zeros(size(info.cradleimg));
%     if info.lp(2) < info.rp(2)
%         ind = (info.lp(2):info.rp(2)) - info.colind(1) + 1;
%     else
%         ind = (info.lp(2):-1:info.rp(2)) - info.colind(1) + 1;
%     end
%     for i = 1:length(ind)
%         mask3(up(i):end,ind(i)) = 1;
%     end
%     mask4 = zeros(size(info.cradleimg));
%     if info.ld(2) < info.rd(2)
%         ind = (info.ld(2):info.rd(2)) - info.colind(1) + 1;
%     else
%         ind = (info.ld(2):-1:info.rd(2)) - info.colind(1) + 1;
%     end
%     for i = 1:length(ind)
%         mask4(1:down(i),ind(i)) = 1;
%     end
%     info.cradleimg = mask1.*mask2.*mask3.*mask4;

    imgsiz = size(info.cradleimg);
    xv = [info.lp(1) info.ld(1) info.rd(1) info.rp(1)] - info.rowind(1) + 1;
    yv = [info.lp(2) info.ld(2) info.rd(2) info.rp(2)] - info.colind(1) + 1;
    info.cradleimg(min(xv(1),xv(4)):max(xv(2),xv(3)),min(yv(1),yv(4)):max(yv(2),yv(3))) = 1;
    info.cradleimg(max(xv(1),xv(4)):min(xv(2),xv(3)),max(yv(1),yv(4)):min(yv(2),yv(3))) = 0;
    [X,Y] = find(info.cradleimg);
%     [Y,X] = meshgrid(1:imgsiz(2),1:imgsiz(1));
    In = inpolygon(X(:),Y(:),xv,yv);
    info.cradleimg(sub2ind(imgsiz,X,Y)) = In;
    info.cradleimg(max(xv(1),xv(4)):min(xv(2),xv(3)),max(yv(1),yv(4)):min(yv(2),yv(3))) = 1;
%     info.cradleimg = reshape(In,imgsiz);
    vinfonew{j} = info;
end