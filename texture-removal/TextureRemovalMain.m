opt.transform = 'shearlet';
opt.direction = 'horizontal';

% imgind = [1,size(img,1);1,size(img,2)];
% opt.mask_h = info2mask(hinfo,imgind);
% opt.mask_v = info2mask(vinfonew,imgind);
% hbdImg = BoundaryImg(img,opt);

% === remove horizontal boundary === %
[subimg,subind] = subpatch4HorBdryRm(img,hinfo,vinfonew);
L = size(subimg{1}{1},1);
subopt.transform = 'shearlet';
subopt.direction = 'horizontal';
hbdImg = cell(size(subimg));
ir = 0;
for i1 = 1:length(subimg)
    if mod(i1,2) == 1
        ir = ir+1;
        local_v1 = cellfun(@(x)[x.ld(2) x.rd(2)],vinfo(ir,:),'UniformOutput',0);
        local_v1 = cell2mat(local_v1(:));
        local_v2 = cellfun(@(x)[x.lp(2) x.rp(2)],vinfo(ir+1,:),'UniformOutput',0);
        local_v2 = cell2mat(local_v2(:));
        local_v = [min(local_v1(:,1),local_v2(:,1)) max(local_v1(:,2),local_v2(:,2))];
        local_cutgap = cellfun(@(x)x.cutwidth,crossinfo(ir,:),'UniformOutput',0);
        local_cutgap = cell2mat(local_cutgap(:));
    end
   
    hbdImg{i1} = cell(size(subimg{i1}));
    for i2 = 1:length(subimg{i1})
        if ~isempty(subimg{i1}{i2})
            subimgind = subind{i1}(:,i2)*ones(1,2)+[0,L-1;0,L-1];
            subopt.mask_h = info2mask(hinfo,subimgind);
            if sum((subimgind(2,1) < local_v(:,2)) & (subimgind(2,2) > local_v(:,1))) == 0
                subopt.mask_v = zeros(size(subimg{i1}{i2}));
                subopt.cutwidth = 0;
            else
                vind = find((subimgind(2,1) < local_v(:,2)) & (subimgind(2,2) > local_v(:,1)));
                subopt.cutwidth = max(local_cutgap(vind,:),[],1);
                subopt.mask_v = info2mask(vinfo,subimgind);
            end
            hbdImg{i1}{i2} = BoundaryImg(subimg{i1}{i2},subopt);
        else
            hbdImg{i1}{i2} = []
        end
    end
end

HbdImg = cell(size(hbdImg));
HbdImg = arrayfun(@(x)block2img(hbdImg{x},subind{x}(1,:),subind{x}(2,:)),1:length(hbdImg),'UniformOutput',0);
imgnew = img;
for i = 1:length(HbdImg)
    ind1 = subind{i}(1,1);
    ind2 = ind1 + size(HbdImg{i},1) - 1;
    imgnew(ind1:ind2,:) = imgnew(ind1:ind2,:) - HbdImg{i};
end

figure;imshow([img,imgnew],[]);
img = imgnew;

% === remove vertical boundary === %
tic
subopt = struct();
subopt.transform = 'shearlet';
[subimg,subind] = subpatch4VerBdryRm(img,hinfo,vinfonew);
subopt.direction = 'vertical';
subopt.bdwidthIn = 15;
subopt.bdwidthOut = 15;
vbdImg = cell(size(subimg));
ic = 0;
for i1 = 1:length(subimg)
    if mod(i1,2) == 1
        ic = ic +1;
    end
    vbdImg{i1} = cell(size(subimg{i1}));
    for i2 = 1:length(subimg{i1})
        if ~isempty(subimg{i1}{i2})
            subimgind = subind{i1}(:,i2)*ones(1,2) + [0,L-1;0,L-1];
            subopt.mask_v = info2mask(vinfonew,subimgind);
            vbdImg{i1}{i2} = BoundaryImg(subimg{i1}{i2},subopt);
        else
            vbdImg{i1}{i2} = [];
        end
    end
end

VbdImg = cell(size(vbdImg));
VbdImg = arrayfun(@(x)block2img(vbdImg{x},subind{x}(1,:),subind{x}(2,:)),1:length(vbdImg),'UniformOutput',0);
imgnew = img;
for i = 1:length(VbdImg)
    ind1 = subind{i}(2,1);
    ind2 = ind1 + size(VbdImg{i},2) - 1;
    imgnew(:,ind1:ind2) = imgnew(:,ind1:ind2) - VbdImg{i};
end
toc
figure;imshow([img,imgnew],[]);
img = imgnew;

clear subimg subind subopt

% opt.direction = 'vertical';
% vbdImg = BoundaryImg(img,opt);
% img = img - hbdImg - vbdImg;

[subimg,orow,ocol,rowind,colind] = img2subpatch(img,512);

MCAdecomp;
tic
img = parts(:,:,2);
ST = shearletTransformSpect(img);
ST = real(ST);
[j,k,cone] = arrayfun(@(x)shearletScaleShear(x), 1:size(ST,3));
opt.imgsize = size(img);
opt.direction = 'vertical';
[STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
[featureST,opt] = getSTfeature(STnew,j,opt);
Y = [featureST.cradle;featureST.free];
originImg = invfeatureST(Y,opt,1);
opt.originImg = originImg;
[Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
opt.matname = '';
opt.method = 'supervised';
opt.updaterho = 1;
toc

opt.updaterho = 1;
opt.mrho = 0.5;
for i = 2%1:4
    for jj = 1%1:2
        img = subpatch{i,jj};
        imgind = [rowind(i),rowind(i)+511;colind(jj),colind(jj)+511];
        opt.mask_v = info2mask(vinfonew,imgind);
        ST = shearletTransformSpect(img);
        ST = real(ST);
        [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
        [featureST,opt] = getSTfeature(STnew,j,opt);
        Y = [featureST.cradle;featureST.free];
        originImg = invfeatureST(Y,opt,1);
        opt.originImg = originImg;
        [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
        spfactcovest_mgploadings(Y,mask,opt,Z,loc);
        load('output_1.mat','*Img');
        vresult{i,jj} = cradleImg + residualImg.*opt.mask_v;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
stage1_result.img{n,3} = img;
stage1_result.hinfo{n,3} = hinfo;
stage1_result.vinfo{n,3} = vinfo;
stage1_result.opt{n,3} = opt;
stage1_result.vinfonew{n,3} = vinfonew;


I = 1;
J = 1;
img = stage1_result.img{I,J};
hinfo = stage1_result.hinfo{I,J};
vinfo = stage1_result.vinfonew{I,J};
opt = stage1_result.opt{I,J};

[subpatch,orow,ocol,rowind,colind] = img2subpatch(img,512);
size(block2img(subpatch,orow,ocol)) - size(img)

vbdImg = cell(size(subpatch,1),2);
for jj = 1:3
    for i = 1:size(subpatch,1)
        img = subpatch{i,jj};
           imgind = [rowind(i) (rowind(i)+size(img,1)-1);colind(jj) (colind(jj)+size(img,2)-1)];
        opt.mask_h = info2mask(hinfo,imgind);
        opt.mask_v = info2mask(vinfonew,imgind);
        opt.direction = 'vertical';
     vbdImg{i,jj} = BoundaryImg(img,opt);
    end
end
vbdImg = block2img(vbdImg,orow,ocol);
figure;imshow(stage1_result.img{I,J} - vbdImg,[])


hbdImg = cell(size(subpatch));
for i = 1:size(subpatch,1)
for jj = 1:size(subpatch,2)
     img = subpatch{i,jj};
           imgind = [rowind(i) (rowind(i)+size(img,1)-1);colind(jj) (colind(jj)+size(img,2)-1)];
        opt.mask_h = info2mask(hinfo,imgind);
        opt.mask_v = info2mask(vinfo,imgind);
        opt.direction = 'horizontal';
     hbdImg{i,jj} = BoundaryImg(img,opt);
end
end
for jj = 1:size(subpatch,2)
    subpatch{i,jj} = subpatch{i,jj} - hbdImg{jj};
end

img = block2img(subpatch,orow, ocol);
        filepath = ['/home/grad/rachel/Documents/MATLAB/18a/subimg_' num2str(I) '_' num2str(J) '/'];
        if (~exist(filepath, 'dir'))
            mkdir(filepath);
        end

cartoon = cell(size(subpatch));
texture = cell(size(subpatch));
for i = 1:size(subpatch,1)
    for j = 1:size(subpatch,2)
        img = subpatch{i,j};
        MCAdecomp;
        cartoon{i,j} = parts(:,:,1);
        texture{i,j} = parts(:,:,2);
    end
end
save([filepath 'result.mat'],'cartoon','texture');

%% == texture separation == %

filepath = '/home/grad/rachel/Documents/MATLAB/MCAdecomp_Ghent_8/';

% vresult = cell(size(texture));
% hresult = cell(size(texture));
% for i = 1:size(hresult,1), for jj = 1:size(hresult,2), hresult{i,jj} = zeros(size(texture{1}));end;end

Npatch = size(texture);
hcradle_sample_size = zeros(Npatch);
hnoncradle_sample_size = zeros(Npatch);
vcradle_sample_size = zeros(Npatch);
vnoncradle_sample_size = zeros(Npatch);

for i = 1:size(texture,1)
    for jj = 1:size(texture,2)
        img = texture{i,jj};
        imgind = [rowind(i) (rowind(i)+size(img,1)-1);colind(jj) (colind(jj)+size(img,2)-1)];
        opt.mask_h = info2mask(hinfo,imgind);
        opt.mask_v = info2mask(vinfo,imgind);
        ST = shearletTransformSpect(img);
        ST = real(ST);
%         [j,k,cone] = arrayfun(@(x)shearletScaleShear(x), 1:size(ST,3));
        opt.imgsize = size(img);
        
        opt.direction = 'vertical';
        opt.lowfreqfeature = 1;
        [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
        [featureST,opt] = getSTfeature(STnew,j,opt);
        vcradle_sample_size(i,jj) = size(featureST.cradle,1);
        vnoncradle_sample_size(i,jj) = size(featureST.free,1);
        save([filepath 'vfeatureST_' num2str(i) '_' num2str(jj) '.mat'],'featureST','opt');
%         Y = [featureST.cradle;featureST.free];
%         originImg = invfeatureST(Y,opt,1);
%         opt.originImg = originImg;
%         [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
%         opt.matname = '/home/grad/rachel/matlab/cradleremoval/output';
%         opt.method = 'supervised';
%         opt.mrho = .5;
%         spfactcovest_mgploadings(Y,mask,opt,Z,loc);
%         load('output_1.mat','*Img');
%         vresult{i,jj} = cradleImg + residualImg.*opt.mask_v;
%         save([path 'vresult.mat'],'vresult');
                opt.direction = 'horizontal';
                opt.lowfreqfeature = 1; % whether to include low frequency feature or not
        [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
        [featureST,opt] = getSTfeature(STnew,j,opt);
        hcradle_sample_size(i,jj) = size(featureST.cradle,1);
        hnoncradle_sample_size(i,jj) = size(featureST.free,1);
%         if sum(opt.mask_h(:)) < length(img(:))*.1
%             featureST = [];
%             hresult{i,jj} = zeros(size(img));
%             continue;
%         end
        save([filepath 'hfeatureST_' num2str(i) '_' num2str(jj) '.mat'],'featureST','opt');
%         
%         Y = [featureST.cradle;featureST.free];
%         originImg = invfeatureST(Y,opt,1);
%         opt.originImg = originImg;
%         [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
%         opt.matname = '/home/grad/rachel/matlab/cradleremoval/output';
%         opt.method = 'supervised';
%         opt.mrho = .8;
%         spfactcovest_mgploadings(Y,mask,opt,Z,loc);
%         load('output_1.mat','*Img');
%         hresult{i,jj} = cradleImg + residualImg.*opt.mask_h;
%         save([path 'hresult.mat'],'hresult');
    end
end

save([filepath 'Ghent_8.mat'],'-append','hcradle_sample_size','hnoncradle_sample_size','vcradle_sample_size','vnoncradle_sample_size');


%% run Input SampleNeighborPatch.m
num_NCsample = 500; % number of non-cradle samples draw per region
num_Csample = 1500; % number of cradle samples draw per region

for i = 2%1:size(texture,1)
    Ync = [];
    Yc = [];
    for jj = 1:size(texture,2)
        try
            load([filepath 'featureST_' num2str(i) '_' num2str(jj) '.mat'])
            Ync = [Ync; featureST.free(randsample(1:size(featureST.free,1),num_NCsample),:)];
            Yc = [Yc; featureST.cradle(randsample(1:size(featureST.cradle,1),num_Csample),:)];
        end
    end
    Y = [Yc; Ync];
    mask = [ones(size(Yc,1),1);zeros(size(Ync,1),1)];
    Z = [];
    opt.mrho = .8;
    spfactcovest_mgploadings(Y,mask,opt,Z);
    for jj = 1:size(texture,2)
        load([filepath 'featureST_' num2str(i) '_' num2str(jj) '.mat'],'featureST','opt');
        Y = [featureST.cradle;featureST.free];
        originImg = invfeatureST(Y,opt,1);
        opt.originImg = originImg;
        [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
        zout = spfactor_post_inference([opt.matname '_1.mat'],Y,opt);
        hresult{i,jj} = zout.cradleImg + (originImg -zout.noncradleImg).*opt.mask_h;
    end
end

cartoon = block2img(cartoon,orow,ocol);
texture = block2img(texture,orow,ocol);
htexture = block2img(hresult,orow,ocol);
vtexture = block2img(vresult,orow,ocol);
result = stage1_result.img{I,J} - htexture - vtexture;
save([filepath 'result.mat'],'cartoon','texture','htexture','vtexture','result');
