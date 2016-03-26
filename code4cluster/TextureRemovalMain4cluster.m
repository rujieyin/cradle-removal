%% set up option, add field to opt

opt.transform = 'shearlet';
opt.direction = 'horizontal';

vinfonew = connect_vinfo(vinfo);
img = Inew;

%% === remove horizontal boundary === % (local machine)
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
        try
        local_cutgap = cellfun(@(x)x.cutwidth,crossinfo(ir,:),'UniformOutput',0);
        local_cutgap = cell2mat(local_cutgap(:));
        catch
            edge = [hinfo{ir}.lp(1)+hinfo{ir}.rp(1),hinfo{ir}.ld(1)+hinfo{ir}.rd(1)]/2;
            cut1 = cellfun(@(x)(x.downcutind + x.rowind(1)),vinfo(ir,:),'UniformOutput',0);
            cut2 = cellfun(@(x)(x.upcutind),vinfo(ir+1,:),'UniformOutput',0);
%             edge = [hinfo{ir}.lp + hinfo{ir}.rp, hinfo{ir}.ld + hinfo{ir}.rd]/2 + hinfo{ir}.rowind(1);
%             cut1 = cellfun(@(x)(x.downcutind + x.rowind),vinfo(ir,:),'UniformOutput',0);
%             cut2 = cellfun(@(x)(x.upcutind + x.rowind), vinfo(ir+1,:),'UniformOutput',0);
            local_cutgap = max([edge(1) - cell2mat(cut1(:)), cell2mat(cut2(:)) - edge(2)],[],2);
        end
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

%% === remove vertical boundary === % (local machine)

subopt = struct();
subopt.transform = 'shearlet';
[subimg,subind] = subpatch4VerBdryRm(img,hinfo,vinfonew);
subopt.direction = 'vertical';
subopt.bdwidthIn = 5;
subopt.bdwidthOut = 5;
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

figure;imshow([img,imgnew],[]);
img = imgnew;

%% clean up the stage of boundary image removal
clear subimg subind subopt imgnew

%% cut image for MCA decompositoin
[subimg,orow,ocol,rowind,colind] = img2subpatch(img,512);


%% == MCA = % (cluster)
 MCAdecomp_cluster_building

% you may want to check if the MCA result is correct before next step

%% == texture separation == % (cluster)

filepath = '/home/grad/rachel/Documents/MATLAB/MCA/';

runGrid; % on  computeShearletDescriptor_cluster

% these paths should be the same as in the 'runGrid.m'
load_mat = '/home/grad/rachel/Documents/MATLAB/MCA/ghent-multiplicative.mat';
result_path = ['/home/grad/rachel/Documents/MATLAB/MCA/ghent-new-model/'];
% colllect the MCA texture result from the result_path
subimg = load(load_mat,'subimg');
texture = cell(size(texture));
%======= write a for loop here to load texture result from file according
%to job_id and assign the result to the cell matrix texture ============%


runGrid; % on computeBayesianSeparation_cluster

%% == assemble the blocks to final result == % (local machine)

img = block2img(subimg,rowind, colind);
