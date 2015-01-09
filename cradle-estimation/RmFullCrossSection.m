function [imgnew,info] = RmFullCrossSection(img,Vinfo,hloc,opt)

info = cell(size(Vinfo,1)-1,size(Vinfo,2));

for i = 1%1:size(info,1)
    for j = 1%1:size(info,2)
        ind1 = Vinfo{i,j}.ind(2) - 2*opt.vs;%max([Vinfo{i,j}.maskld(1),Vinfo{i,j}.maskrd(1)]);
        ind2 = Vinfo{i+1,j}.ind(1) + 2*opt.vs;%min([Vinfo{i+1,j}.masklp(1),Vinfo{i+1,j}.maskrp(1)]);
        ind3 = max([Vinfo{i,j}.ind(3),Vinfo{i+1,j}.ind(3)]);
        ind4 = min([Vinfo{i,j}.ind(4),Vinfo{i+1,j}.ind(4)]);
        bdinfo = struct();
        bdinfo.up = max([Vinfo{i,j}.maskld(1),Vinfo{i,j}.maskrd(1)]) - ind1+1;
        bdinfo.down = min([Vinfo{i+1,j}.masklp(1),Vinfo{i+1,j}.maskrp(1)]) - ind1 +1;
        hangle = hloc{i}.angle;
        vangle = Vinfo{i,j}.angle;
        [subimgnew,info{i,j}] = RmCrossSection(img(ind1:ind2,ind3:ind4),hangle,vangle,bdinfo);    
        info{i,j}.ind = [ind1,ind2,ind3,ind4];
        imgnew(ind1:ind2,ind3:ind4) = subimgnew;
    end
end
