function [feature,opt] = getCWfeature(cw,opt)

% each column is a feature vector

if strcmp(opt.direction,'horizontal')
    if isfield(opt,'hangleind')
        angleind = opt.hangleind;
    else
        angleind = cellfun(@(x)[1:floor(length(x)/4),floor(length(x)/2)+1:floor(length(x)/4)*3],cw,'UniformOutput',0);
        opt.hangleind = angleind;
    end
    mask = opt.mask_h;
    if isfield(opt,'hdownsamplesize')
        downsamplesize = opt.hdownsamplesize;
    else
        if isfield(opt,'htexturesize')
            downsamplesize = round(opt.imgsize(:)./opt.htexturesize(:));
        else
            nangle = length(cw{end});
            downsamplesize = size(cw{end}{round(nangle/8)});
        end
        opt.hdownsamplesize = downsamplesize;
    end
    % choose the coarsest level to be considered by constrain on the ratio
    % of width and height of cw block
    startind = 2;
    while startind < length(cw)
        siz = size(cw{startind}{angleind{startind}(1)});
        if siz(1)/siz(2) > downsamplesize(1)/downsamplesize(2)/3
            break;
        else
            startind = startind + 1;
        end
    end
else
    if isfield(opt,'vangleind')
        angleind = opt.vangleind;
    else
        angleind = cellfun(@(x)[floor(length(x)/4)+1:floor(length(x)/2),floor(length(x)/4)*3+1:floor(length(x)/2)*2],cw,'UniformOutput',0);
        opt.vangleind = angleind;
    end
    mask = opt.mask_v;
    if isfield(opt,'vdownsamplesize')
        downsamplesize = opt.vdownsamplesize;
    else
        if isfield(opt,'vtexturesize')
            downsamplesize = round(opt.imgsize(:)./opt.vtexturesize(:));
        else
            nangle = length(cw{end});
            downsamplesize = size(cw{end}{round(nangle/8*3)});
        end
        opt.vdownsamplesize = downsamplesize;
    end
    % choose the coarsest level
    startind = 2;
    while startind < length(cw)
        siz = size(cw{startind}{angleind{startind}(1)});
        if siz(2)/siz(1) > downsamplesize(2)/downsamplesize(1)/3
            break;
        else
            startind = startind + 1;
        end
    end
end

mask = imresize(mask,downsamplesize,'nearest');
maskcradle = mask > .5;
maskfree = ~(maskcradle);
% maskcradle = zeros(size(mask));
% maskcradle(370:end,:) = 1;
% maskfree = zeros(size(mask));
% maskfree(1:330,:) = 1;
% maskcradle = logical(imresize(maskcradle,downsamplesize));
% maskfree = logical(imresize(maskfree,downsamplesize));
% figure;imagesc([maskfree,maskcradle])

% low frequency excluded as it is zero
feature = cell(length(cw),1);
for i = startind:length(feature)
    feature{i} = cellfun(@(x)reshape(imresize(x,downsamplesize,'nearest'),1,[]),cw{i}(angleind{i}),'UniformOutput',0);
    feature{i} = cell2mat(feature{i}(:));
    if ~opt.curveletisreal | (isfield(opt,'absfeature') && opt.absfeature)
        feature{i} = abs(feature{i});
    end
end
feature = cell2mat(feature(:));

figure;imagesc(reshape(max(feature.*conj(feature),[],1),downsamplesize));

[energy,ind] = sort(sum(feature.*conj(feature),1),'descend');
% figure;imshow(reshape(sum(feature.^2,1),downsamplesize),[])
tmp = cumsum(energy);
cutind = find(tmp/tmp(end) > .95,1);
ind = ind(1:cutind);
tmp = zeros(downsamplesize);
tmp(ind) = energy(1:cutind);
figure;imshow(tmp,[]);title('95% energy retained');
energy = sqrt(energy(1:cutind));
feature = feature(:,ind);
feature1 = feature(:,maskcradle(ind));
feature2 = feature(:,maskfree(ind));
feature = feature./repmat(energy(:)',[size(feature,1),1]);
feature3 = feature(:,maskcradle(ind));
feature4 = feature(:,maskfree(ind));

feature = struct();
if strcmp(opt.direction,'horizontal')
    feature.hcradle = feature1;
    feature.nonhcradle = feature2;
    feature.hcradle_normalized = feature3;
    feature.nonhcradle_normalized = feature4;
    opt.hfeatureind = struct();
    % ind in the downsamplesize frame
    opt.hfeatureind.cradle = ind(maskcradle(ind));
    opt.hfeatureind.free = ind(maskfree(ind));
else
    feature.vcradle = feature1;
    feature.nonvcradle = feature2;
    feature.vcradle_normalized = feature3;
    feature.nonvcradle_normalized = feature4;
    opt.vfeatureind = struct();
    opt.vfeatureind.cradle = ind(maskcradle(ind));
    opt.vfeatureind.free = ind(maskfree(ind));
end
opt.featureangleind = angleind;
opt.featureangleind(1:(startind-1)) = cell(startind-1,1);


    