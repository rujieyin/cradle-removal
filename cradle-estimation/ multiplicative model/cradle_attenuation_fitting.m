function [imgnew,info,p,in,out] = cradle_attenuation_fitting(img, midpoint ,angle, verest,estimation)
% this function fit the attenuation caused by cradling based on a
% multiplicative model

% Input:    img: an (sub-)image containing only one HORIZONTAL cradle
%           midpoint: row indices of middle points of cradle boundaries
%           angle: the angle of the cradle
%           verest: estimated location of any vertical cradle, need to be
%                   excluded in the estimation
%           
% Output:   imgnew: attenuation recovered image
%           info: the location of corners of the cradle, the edge profiles
%           p: parameters of linear model, p(1) is slope, p(2) is intercept
%           in, out: data points for linear model fitting
%
% The Model: img(panel) = img(cradled)/C0 + C(1-1/C0), where 0 < C0 < 1 is
%                   attenuation factor, and C is the saturation value

%%
display = 0;
if nargin < 5
estimation = 'section';
end

% convert the angle
angle = angle/180*pi;

%% find the boundary segments for estimation

% width of boundary neighbourhood
width = 30;
% length of cradle
len = size(img,2);

% compute the coordinates (row indices)
x = (1:len) - ceil(len/2);% column
y1 = midpoint(1) + round(cos(angle)*x);% row
y2 = midpoint(2) + y1 - midpoint(1);% row

info = struct();
info.lp = y1(1);
info.rp = y1(end);
info.ld = y2(1);
info.rd = y2(end);

%available segments
verest = reshape(verest,2,[]);
nv = size(verest,1); % number of vertical cradles
if nv > 0
    seg = cell(1,nv+1);
    seg{1} = [1 : verest(1)-1];%1:lenseg:(verest(1)-lenseg);%
    for i = 2:nv
        seg{i} = [verest(2,i-1)+1 : verest(1,i)-1];%verest(2,i-1):lenseg:(verest(1,i)-lenseg);%
    end
    seg{nv+1} = [verest(2,nv)+1 : len];%verest(2,nv):lenseg:length;%
else
    seg = {[1:len]};
end
% get segments of data
dataseg = @(d)cellfun(@(x)d(:,x),seg,'UniformOutput',0);

%% compute the average of inside/outside boundary value

% length of segments
lenseg = 40;

% boundary slices
img_out1 = @(x,s)img((y1(x)-s):(y1(x)-1),x);
img_in1 = @(x,s)img((y1(x)):(y1(x)+s),x);
img_in2 = @(x,s)img((y2(x)-s):(y2(x)-1),x);
img_out2 = @(x,s)img((y2(x)):(y2(x)+s),x);

imgout1 = arrayfun(@(x)img_out1(x,width),1:len,'UniformOutput',0);% pixel value above point img(y1(x),x)
imgin1 = arrayfun(@(x)img_in1(x,width),1:len,'UniformOutput',0);% below
imgout2 = arrayfun(@(x)img_out2(x,width),1:len,'UniformOutput',0);% pixel value above point img(y2(x),x)
imgin2 = arrayfun(@(x)img_in2(x,width),1:len,'UniformOutput',0);% below

% smooth data
smfilter = @(x)imfilter(x,ones(1,lenseg)/lenseg,'symmetric');

bdout1 = cellfun(@median,imgout1);
bdout2 = cellfun(@median,imgout2);
out = [dataseg(bdout1) dataseg(bdout2)];% discard vertical cradled part
out = cellfun(smfilter,out,'UniformOutput',0);
out = cell2mat(out(:)');

bdin1 = cellfun(@median, imgin1);
bdin2 = cellfun(@median, imgin2);
in = [dataseg(bdin1) dataseg(bdin2)]; % discard vertical cradled part
in = cellfun(smfilter,in,'UniformOutput',0);
in = cell2mat(in(:)');


%% fitting the linear model
p = polyfit(in,out,1);

if p(1) < 1
    p(1) = 1.01;
    p(2) = median(out)-p(1)*median(in);
end

%if display
% plot the linear model
figure;scatter(in,out,'*');hold on;plot(0:.1:max(img(:))+10,polyval(p,0:.1:max(img(:))+10),'LineWidth',2);hold off;
set(gca, 'fontsize',16);
xlabel('\fontsize{16} attenuated pixel value');
ylabel('\fontsize{16} original pixel value');
title(['\fontsize{16} Linear Model of Attenuation, p(1) = ' num2str(p(1)) ', p(2) = ' num2str(p(2))]);
%end

C = p(2)/(1-p(1));
%C0 = 1/p(1);

%% estimate the profile of edge p(1)

% get boundary slices
bdslice1 = [cell2mat(imgout1);cell2mat(imgin1)];
bdslice1 = dataseg(bdslice1);
bdslice1 = cell2mat(cellfun(smfilter,bdslice1,'UniformOutput',0));

bdslice2 = [cell2mat(imgin2);cell2mat(imgout2)];
bdslice2 = dataseg(bdslice2);
bdslice2 = cell2mat(cellfun(smfilter,bdslice2,'UniformOutput',0));

% transform data
%transform = @(x,base)(C-x)/base; % base = C - img(panel)

C0profile1 = (C-bdslice1)./repmat(C-out(1:size(bdslice1,2)),[size(bdslice1,1),1]);
var1 = var(C0profile1,0,2);
C0profile1 = mean(C0profile1,2);


C0profile2 = (C-bdslice2)./repmat(C-out(end-size(bdslice2,2)+1:end),[size(bdslice2,1),1]);
var2 = var(C0profile2,0,2);
C0profile2 = mean(C0profile2,2);

if display
    figure; errorbar(1:size(var1,1),C0profile1,var1);
    figure; errorbar(1:size(var2,1),C0profile2,var2);
end

switch estimation
    case 'edge'
        
        p1profile = [1./C0profile1;p(1)*ones(midpoint(2)-midpoint(1)-width*2-1,1);1./C0profile2];        
        
    case 'section'
        %% estimate the profile of full cross-section
        img_slice = @(x,s)img((y1(x)+s+1):(y2(x)-s-1),x);
        imgslice = arrayfun(@(x)img_slice(x,width),1:len,'UniformOutput',0);% pixel value above point img(y1(x),x)
        smtrans = (cos(linspace(0,pi,midpoint(2)-midpoint(1)-2*width-1))+1)/2;
        outmedian1 = median(out(1:size(out,2)/2));
        outmedian2 = median(out(size(out,2)/2+1:end));
        out = ones(1,size(smtrans,2))*outmedian2+smtrans*(outmedian1-outmedian2);
        
        imgslice = dataseg(cell2mat(imgslice));
        imgslice = cell2mat(cellfun(smfilter,imgslice,'UniformOutput',0));
        
        C0profile = (C-imgslice)./repmat(C-out(:),[1,size(imgslice,2)]);
        Var = var(C0profile,0,2);
        C0profile = mean(C0profile,2);
        
        if display
            figure; errorbar(1:size(Var,1),C0profile,Var);
            figure; plot(1:size(Var,1),mean(imgslice,2),...
                1:size(Var,1),mean(imgslice,2)+var(imgslice,0,2),...
                1:size(Var,1),mean(imgslice,2)-var(imgslice,0,2));
        end
        
        p1profile = [1./C0profile1 ; 1./C0profile ; 1./C0profile2];        
        
end

% 1d medianfilter to smooth the signal while keep the jump
tmp = medfilt1(p1profile,width);
p1profile(width:end-width) = tmp(width:end-width);

info.p1profile = p1profile;
info.C = C;
info.width = width;
info.midpoint = midpoint;

imgnew = RemoveAttenuationProfile(img,midpoint,angle,width,p1profile,C);


%         p1profile = [p1profile; ones(size(img,1)-size(p1profile,1),1)];
%         p2profile = [p2profile; zeros(size(img,1)-size(p2profile,1),1)];
% 
% figure;plot(1:length(p1profile),p1profile);
% figure;plot(1:length(p2profile),p2profile);
% 
% %% compute the linear transform
% 
% transformSlice = @(x,y)(circshift(p1profile,y-width-1).*x + circshift(p2profile,y-width-1));
% imgnew = num2cell(img,1);
% imgnew = cellfun(@(x,y)transformSlice(x,y),imgnew,num2cell(y1),'UniformOutput',0);
% imgnew = cell2mat(imgnew);

if display
    figure;imshow(imgnew,[]);
end


