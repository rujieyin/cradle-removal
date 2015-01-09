function [verest,horest, signal1, signal2] = cradledetect(img,vn,vrange,hn,hrange,opt)

% input: vn:          # of vertical edges in cradle
%        vrange:      a 1 by 2 array indicating the range of position of vertical edges
%        hn:          # of horizontal edges in cradle
%        hrange:      a 1 by 2 array indicating the range of position of
%                     horizontal edges

if nargin > 5
    if isfield(opt,'L')
        L = opt.L; % haarfilter length
    else
L = 5;
    end
    if isfield(opt,'s')
        s = opt.s;
    else
s = 2; % smoothing parameter
    end
end

haarfilt = [ones(1,L),-ones(1,L)];
vert = imfilter(img,haarfilt,'symmetric');
signal1 = smooth(sum(vert,1),s);
[peak,loc] = findpeaks(signal1.*(signal1>0));
% signal1 = smooth(sum(abs(vert),1),s);
% [peak, loc] = findpeaks(signal1);
peak = peak(loc > vrange(1) & loc < vrange(2));
loc = loc(loc > vrange(1) & loc < vrange(2));
[~, ind] = sort(peak,'descend');
loc = loc(ind(1:vn));
verest = loc;
[peak,loc] = findpeaks(-signal1.*(signal1<0));
peak = peak(loc > vrange(1) & loc < vrange(2));
loc = loc(loc > vrange(1) & loc < vrange(2));
[~, ind] = sort(peak,'descend');
loc = loc(ind(1:vn));
verest = [verest;loc];
verest = sort(verest,'ascend');

mask = ones(1,size(img,2));
for i = 1:vn
    mask(verest(2*i-1):verest(2*i)) = 0;
end
mask = logical(mask);
img = img(:,mask);

horiz = imfilter(img,haarfilt','symmetric');
signal2 = smooth(sum(horiz,2),s);
[peak,loc] = findpeaks(signal2.*(signal2>0));
% signal2 = smooth(sum(abs(horiz),2),s);
% [peak, loc] = findpeaks(signal2);
peak = peak(loc > hrange(1) & loc < hrange(2));
loc = loc(loc > hrange(1) & loc < hrange(2));
[~, ind] = sort(peak,'descend');
loc = loc(ind(1:hn));
horest = loc;
[peak,loc] = findpeaks(-signal2.*(signal2<0));
peak = peak(loc > hrange(1) & loc < hrange(2));
loc = loc(loc > hrange(1) & loc < hrange(2));
[~, ind] = sort(peak,'descend');
loc = loc(ind(1:hn));
horest = [horest;loc];
horest = sort(horest,'ascend');
