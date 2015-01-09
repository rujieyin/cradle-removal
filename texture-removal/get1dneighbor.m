   function neighind = get1dneighbor(i,j,m,n,dir)
   % this function calculate the 1D neighbor of (i,j) inside the full (m,n)
   % for direction 'horizontal' or 'vertical'
        
        switch dir
            case 'vertical'
                if i > 1 && i < m
                    neighind = [i-1, i, i+1;j,j,j];
                else
                    if i == 1
                        neighind = [1,2;j,j];
                    else
                        neighind = [(m-1):m;j,j];
                    end
                end
            case 'horizontal'
                if j > 1 && j < n
                    neighind = [i,i,i;j-1,j,j+1];
                else
                    if j == 1
                        neighind = [i,i;1,2];
                    else
                        neighind = [i,i;n-1:n];
                    end
                end
        end
        neighind = sub2ind([m,n],neighind(1,:),neighind(2,:));
        
    end
