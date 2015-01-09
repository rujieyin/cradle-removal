function [] = displayfeatureintensity(feature,opt,img,cw)

if strcmp(opt.direction,'horizontal')
    %     cradleftr = feature.hcradle;
    %     noncradftr = feature.nonhcradle;
    % %     cradleftr = feature.hcradle_normalized;
    % %     noncradftr = feature.nonhcradle_normalized;
    %     cradleind = opt.hfeatureind.cradle;
    %     noncradind = opt.hfeatureind.free;
    %     siz = opt.hdownsamplesize;
    ind = opt.hangleind;
    mask = opt.mask_h;
else
    %     cradleftr = feature.vcradle;
    %     noncradftr = feature.nonvcradle;
    % %     cradleftr = feature.vcradle_normalized;
    % %     noncradftr = feature.nonvcradle_normalized;
    %     cradleind = opt.vfeatureind.cradle;
    %     noncradind = opt.vfeatureind.free;
    %     siz = opt.vdownsamplesize;
    ind = opt.vangleind;
    mask = opt.mask_v;
end
fullsiz = size(img);
if nargin > 2
    img = (img - min(img(:)))/(max(img(:)) -  min(img(:)));
    img = repmat(img,[1,1,3]);
    
    for i = 1:length(ind)
        if ~isempty(ind{i})
            siz = zeros(1,2);
            siz(1) = min(arrayfun(@(x)size(cw{i}{x},1),ind{i}));
            siz(2) = min(arrayfun(@(x)size(cw{i}{x},2),ind{i}));
            
            %     img = repmat(imresize(img,siz),[1,1,3]);
            tmpimg = zeros(fullsiz);
            tmpimg = repmat(tmpimg,[1,1,3]);
            %     tmp = zeros(siz);
            %     tmp(cradleind) = 1;
            tmp = mask;
            tmpimg(:,:,1) = tmp;%imresize(tmp,fullsiz,'method','nearest');
            %     tmp = zeros(siz);
            %     tmp(noncradind) = 1;
            tmp = 1 - mask;
            tmpimg(:,:,3) = tmp;%imresize(tmp,fullsiz,'method','nearest');
            alpha = arrayfun(@(x)reshape(imresize(cw{i}{x},siz),[],1),ind{i},'UniformOutput',0);
            alpha = reshape(cell2mat(alpha(:)'),[siz,length(alpha)]);
            alpha = sum(alpha.^2,3);
            alpha = alpha/max(alpha(:))/2;
            % alpha = zeros(siz);
            %     alpha(cradleind) = sum(cradleftr.^2,1)/6e3;
            %     alpha(noncradind) = sum(noncradftr.^2,1)/6e3;
            figure;image(img);axis off;hold on
            h = image(tmpimg);
            set(h,'AlphaData',imresize(alpha,fullsiz,'method','nearest'));
            drawnow;hold off
        end
    end
end


% cut = length(feature.nonhcradle);
% ind1 = ind(ind<=cut);
% ind2 = ind(ind>cut)-cut;
%
% if strcmp(opt.direction,'horizontal')
%     downsamplesize = opt.hdownsamplesize;
%     if nargin < 4
%         img = zeros(downsamplesize);
%     else
%         img = imresize(img,downsamplesize);
%         img = (img - min(img(:)))/(max(img(:))-min(img(:)));
%     end
%     img = repmat(img,[1,1,3]);
%     locind = opt.hfeatureind;
% else
%     downsamplesize = opt.vdownsamplesize;
%     if nargin < 4
%         img = zeros(downsamplesize);
%     else
%         img = imresize(img,downsamplesize);
%         img = (img - min(img(:)))/(max(img(:))-min(img(:)));
%     end
%     img = repmat(img,[1,1,3]);
%     locind = opt.vangleind;
% end
% for i = 1:length(ind1)
%     [m,n] = ind2sub(downsamplesize,locind.free(ind1(i)));
%     img(m,n,:) = [1,0,0];
% end
%
% for i = 1:length(ind2)
%     [m,n] = ind2sub(downsamplesize,locind.cradle(ind2(i)));
%     img(m,n,:) = [0,0,1];
% end
%
% figure;image(img);axis off