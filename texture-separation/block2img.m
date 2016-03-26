function img = block2img(subpatch,rowinfo,colinfo)

if isscalar(rowinfo) && isscalar(colinfo)
rowoverlap = rowinfo;
coloverlap = colinfo;
siz = size(subpatch);
img = cell(siz(1),1);
for j = 1:siz(1)
    img{j} = block2row(subpatch(j,:),coloverlap);
end
img = block2col(img,rowoverlap);
else
rowind = rowinfo;
colind = colinfo;
img = cell(size(subpatch,1),1);
for j = 1:length(img)
    img{j} = subpatch{j,1};
     for ii = 1:(size(subpatch,2)-1)
            coloverlap = size(subpatch{j,ii},2) + colind(ii) - colind(ii+1);
            img{j} = block2row([img(j) subpatch(j,ii+1)],coloverlap);
end
end
imgnew = img(1);
for j = 1:(length(img)-1)
    rowoverlap = size(img{j},1) + rowind(j) - rowind(j+1);
    imgnew{1} = block2col([imgnew;img(j+1)],rowoverlap);
end
img = imgnew{1};
end





    function rowimg = block2row(subpatch,coloverlap)
        nrow = size(subpatch{1},1);
        filter = cdf('normal',1:coloverlap,coloverlap/2,coloverlap/6);
        Lfilter = filter;
        Rfilter = 1-filter;
        fhandle1 = @(x)cat(2,bsxfun(@times,x(:,1:coloverlap),Lfilter),x(:,coloverlap+1:end));
        fhandle2 = @(x)cat(2,x(:,1:end-coloverlap),bsxfun(@times,x(:,(end-coloverlap+1):end),Rfilter));
        subpatch(2:end) = cellfun(fhandle1,subpatch(2:end),'UniformOutput',0);
        subpatch(1:end-1) = cellfun(fhandle2,subpatch(1:end-1),'UniformOutput',0);
        rowimg = subpatch{1};
        for i = 2:length(subpatch)
            rowimg = [rowimg,zeros(nrow,size(subpatch{i},2)-coloverlap)] + [zeros(nrow,size(rowimg,2) - coloverlap),subpatch{i}];
        end
    end

    function colimg = block2col(subpatch,rowoverlap)
        subpatch = cellfun(@(x)x',subpatch,'UniformOutput',0);
        colimg = block2row(subpatch,rowoverlap);
        colimg = colimg';
    end

end