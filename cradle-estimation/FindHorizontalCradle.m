function [loc, angle, R,flag] = FindHorizontalCradle(Img,region,verest,opt)
% this function find the horizontal cradling in Img within reg
% Here we assume there is only one horizontal cradle in the subimage
% Input:    Img:   graysacle image
%           reg:   4by1 array; region to find the horizontal cradling, if is empty,
%                  then find cradling in the whole image
% Output:   loc:   structure; location of the four corners of the cradling,
%                  .lp: leftup  .ld: leftdown .rp: rightup  .rd: rightdown
%           angle: the angle of the cradling

if isempty(region)
    img = Img;
    region = 1;
else
    img = Img(region(1):region(2),region(3):region(4));
    verest = verest - region(3);
end

% theta = 85:.1:95;
% verest = verest + repmat([-1;1]*opt.vs,[length(verest)/2,1]);
% verest = verest(2:end-1);
% angle = zeros(length(verest)/2,1);

% % for i = 1:length(angle)
% [R,xp] = radon(img,theta);
% 
% % % subtract the region of interest
% % if nargin == 5
% %     subreg = estreg - round(size(img,1)/2);
% %     if ~isempty(region)
% %         subreg(1) = subreg(1) - region(1) + 1;
% %         subreg(2) = subreg(2) - region(1) + 1;
% %     end
% %     Rsub = R(xp > subreg(1) - 10 & xp < subreg(2) + 10,:);
% % else
% %     Rsub = R;
% % end
% Rsub = R;
% 
% % compute /ell_1 norm of the difference
% l1energy = sum(abs(imfilter(Rsub,[1;-1])),1);
% [~,ind] = max(l1energy);
% angle = theta(ind);
% end
% angle = mean(angle);

v1 = (1+verest(1))/2;%v1 = (verest(2)+verest(3))/2;
edge1 = findcenteredges(sum(img(:,1:verest(1)),2),opt);%edge1 = findcenteredges(sum(img(:,verest(2):verest(3)),2));

flag = 0;
if edge1 == 0
    flag = 1;
    i = 1;
    while edge1 == 0
        i = i+1;
        edge1 = findcenteredges(sum(img(:,verest(2*i):verest(2*i+1)),2),opt);
    end
    v1 = (verest(2*i)+verest(2*i+1))/2;
end

v2 = (verest(2)+size(img,2))/2;%v2 = (verest(end-2) + verest(end-1))/2;
edge2 = findcenteredges(sum(img(:,verest(2):end),2),opt);%edge2 = findcenteredges(sum(img(:,verest(end-2):verest(end-1)),2));

if edge2 == 0
    flag = 1;
    i = length(verest)/2;
    while edge2 == 0
        i = i-1;
        edge2 = findcenteredges(sum(img(:,verest(2*i-2):verest(2*i-1)),2),opt);
    end
    v2 = (verest(2*i-2)+verest(2*i-1))/2;
end

if v1 >= v2
    disp('Cannot find the cradle!');
    return;
end

alpha = zeros(2,1);
alpha(1) = atan((v2-v1)/(edge2(1)-edge1(1)))/pi*180;
alpha(2) = atan((v2-v1)/(edge2(2)-edge1(2)))/pi*180;
alpha = 180.*(alpha < 0) + alpha;
angle = sum(alpha)/2;

if ~flag
[R,xp] = radon(img,angle);
else
    [R,xp] = radon(img(:,v1:v2),angle);
end

% estimate the location of the horizonal cradling
loc = struct();
% signal = imfilter(R(:,ind),[ones(4,1);-ones(4,1)]);
% % [peaks,Loc] = findpeaks(signal);
% % [~,I] = sort(peaks,'descend');
% % Loc = Loc(I);
% % Loc = Loc.*(R(Loc,ind) > median(R(Loc,ind))/1.5);
% % Loc = Loc(find(Loc,2));
% Loc(1) = findpositivepeaks(signal);
% Loc(2) = length(singal) + 1 - findpositivepeaks(fliplr(-signal));
Loc = findcenteredges(R,opt);
Loc = xp(Loc);

% estimate the cradle location in the middle of img
middle = round((size(img,1)+1)/2);

if Loc(1) > Loc(2)
    middleup = middle - Loc(1);
    middledown = middle - Loc(2);
else
    middleup = middle - Loc(2);
    middledown = middle - Loc(1);
end
% compute the corner location

if ~flag
diff1 = size(img,2)/2*tan((angle-90)/180*pi);
diff2 = diff1;
else
    diff1 = (v1+v2)/2*tan((angle-90)/180*pi);
    diff2 = (size(img,2)-(v1+v2)/2)*tan((angle-90)/180*pi);
end
loc.lp = middleup + round(diff1) + region(1) - 1;
loc.rp = middleup - round(diff2) + region(1) - 1;
loc.ld = middledown + round(diff1) + region(1) - 1;
loc.rd = middledown - round(diff2) + region(1) - 1;
[upline,downline] = getboundaryline(loc,size(img,2));
intensity = median(reshape(Img(loc.lp:loc.ld,:),[],1));

lmin = min(reshape(Img(upline(1):downline(1),region(3):region(3)+ 20),[],1));
for i = 1:size(img,2)
    %     if sum(Img(upline(i):downline(i),i + region(3)-1) > intensity/2) > (downline(i)-upline(i))/2
    if sum(Img(upline(i):downline(i),i + region(3)-1) > (intensity+lmin)/2) > 10
        break;
    end
end
loc.left = i + region(3) -1;
loc.lp = upline(i);
loc.ld = downline(i);
loc.leftT = (lmin+intensity)/2;

rmin = min(reshape(Img(upline(end):downline(end),region(4)-20:region(4)),[],1));
for i = size(img,2) : -1:1
    %     if sum(Img(upline(i):downline(i),i + region(3)-1) > intensity/2) > (downline(i)-upline(i))/2   
    if sum(Img(upline(i):downline(i), i + region(3) - 1) > (intensity+rmin)/2) > 10
        break;
    end
end
loc.right = i + region(3) -1;
loc.rp = upline(i);
loc.rd =downline(i);
loc.rightT = (rmin+intensity)/2;


%     function ll = findpositivepeaks(sig)
%         % this function find the significant peak that is closest to the
%         % center of the signal on the right half part of the signal
%         m = round(length(sig)/2);
%         sig = sig.*(sig > 0);
%         sig = sig(m:end);
%         [pp,ll] = findpeaks(sig);
%         [~,ii] = sort(pp,'descend');
%         ll = min(ll(ii(1:2))) + m - 1;
% 
%     end

 

end
