  function neighind = get2dneighbor(i,j,m,n,r)
  % this function calculates the 2D neighbor of (i,j) (itself NOT INCLUDED!) inside the full (m,n)
  % within radius r. Typically used for choosing sample region of
  % non-cradle feature vectors for the sparse factor model.
  
        mask = zeros(m,n);
        a1 = max(i-r,1);
        a2 = min(i+r,m);
        b1 = max(j-r,1);
        b2 = min(j+r,n);
        mask(a1:a2,b1:b2) = 1;
        mask(i,j) = 0;
        neighind = find(mask);
    end