   function [edge,mu,sigma] = findcenteredges(sig,opt)
        % this function find the edges of the largest intensity part of the
        % signal
        s = 20;%opt.s;
       [m,I] = max(sig);
       subsig = sig(1:I);
       subsig = [zeros(s,1);subsig(:);zeros(s,1)];
       [peaks,loc] = findpeaks(-subsig);
       loc = loc(loc > s & loc < I+s);
       difference = arrayfun(@(x)max(subsig(x-s:x+s)-subsig(x)),loc);
       [d,m1left] = max(difference);
       m1left = loc(m1left) - s;
       subsig = sig(I:end);
       subsig = [zeros(s,1);subsig(:);zeros(s,1)];
       [peaks,loc] = findpeaks(-subsig);
       loc = loc(loc > s & loc < length(subsig) - s);
       difference = arrayfun(@(x)max(subsig(x-s:x+s)-subsig(x)),loc);
       [d,m2right] = max(difference);
       m2right = loc(m2right)+I-21;
       
       if isempty(m1left) | isempty(m2right)
           edge = 0;
           return;
       end
       
%        [~,m1] = min(abs(sig(m1left:I) - (m+sig(m1left))/2));
%        m1 = m1 + m1left - 1;
%        [~,m2] = min(abs(sig(I:m2right) - (m + sig(m2right))/2));
%        m2 = m2 + I - 1;
%        edge = [m1,m2];

%        % find mu and sigma fitting smooth boundary
       m1right = m1left+1;
       while sig(m1right) > sig(m1right-1) 
           m1right = m1right + 1;
       end
%        m1left = m1-1;
%        while sig(m1left) < sig(m1left+1)
%            m1left = m1left - 1;
%        end
       mu1 = floor((m1left + m1right)/2);
       sigma1 = (m1right-m1left)/6;
       m2left = m2right-1;
       while sig(m2left) > sig(m2left+1)
           m2left = m2left - 1;
       end
%        m2right = m2+1;
%        while sig(m2right) < sig(m2right-1)
%            m2right = m2right + 1;
%        end
       mu2 = ceil((m2left+m2right)/2);
       sigma2 = (m2right-m2left)/6;
       edge = [mu1,mu2];
       mu = [mu1,mu2];
       sigma = [sigma1,sigma2];
    end