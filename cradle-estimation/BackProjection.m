% function back = BackProjection(signal,xp,ind)
% % this function back project signal according to the bin index xp and the
% % assigned index for each pixel
% 
% % Input: signal:     1D signal
% %            xp:     bin index
% %           ind:     index matrix
% 
% m = max(ind(:));
% n = min(ind(:));
% l1 = find(xp == m);
% l2 = find(xp == n);
% back = zeros(size(ind));
% for i = n:m
%     back = back + signal(l2 + i - n)*(ind == i)/sum(reshape(ind == i,[],1));
% end

function back = BackProjection(signal,xp,angle,siz)

% Input:  signal:   1D signal
%         xp:       bin index
%         angle:    projection angle in radon transform
%         siz:      the siz of the image

theta = angle/180*pi;
maxbin = floor((siz(1)-1)/2/sin(theta));
back = zeros(siz);
j = find(xp == -maxbin);
for b = -maxbin : maxbin-1
    [up,down,l,r] = PointBackProject(b,angle,siz);
    delta = signal(j)/((r-l+1)*2);
    j = j + 1;
    for i = l : r
        back(up(i-l+1),i) = back(up(i-l+1),i) + delta;
        back(down(i-l+1),i) = back(down(i-l+1),i) + delta;
    end
end
