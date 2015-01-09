function [imgnew,subimgnew,loc,angle] = RmFullHorizontalCradles(img,horest,verest,opt)
% this function removes all horizontal cradles indicated by the estimated
% position of these cradling

% Input: img:    origianl gray scale image
%       horest:  estimation of location of cradling
%       opt:     structure including options

n = length(horest)/2;
loc = cell(n,1);
angle = zeros(n,1);
subimgnew = cell(n,1);
imgnew = img;

for i = n
    region = zeros(1,4);
    region(1) = max(1,horest(2*i-1)-opt.hs);
    region(2) = min(size(img,1),horest(2*i) + opt.hs);
    region(3) = 1;
    region(4) = size(img,2);
    [loc{i},angle(i),R,flag] = FindHorizontalCradle(img,region,verest,opt);
    loc{i}.region = region;
    [~,subimgnew{i}] = RmHorizontalCradle(img,loc{i},angle(i),horest(2*i-1:2*i),opt,R);
    imgnew(region(1):region(2),loc{i}.left:loc{i}.right) = subimgnew{i};
    figure;subplot(2,1,1);imshow(img(region(1):region(2),loc{i}.left:loc{i}.right),[0,255]);title('raw X-ray image');
    subplot(2,1,2);imshow(subimgnew{i},[0,255]);title('horizontal cradle removed X-ray image');
end
