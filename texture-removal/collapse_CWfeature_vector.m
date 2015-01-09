
function cw = collapse_CWfeature_vector(Y,mask,opt)

% this function collapse curvelet feature vectors Y back to a cell array of
% curvelet coefficients

angleind = opt.featureangleind;
if length(cell2mat(angleind')) ~= size(Y,2)
    disp('Dimension of feature vectors is wrong!');
    cw = [];
    return;
end

% ---- check if Y is cradle part ---- %
if size(Y,1) == length(mask)
    flg = 1;
elseif size(Y,1) == sum(mask)
    flg = 2;
end

siz = opt.vdownsamplesize;      % base rescaled size of original image
ncInd = opt.vfeatureind.free;   % index of non-cradle data points
cInd = opt.vfeatureind.cradle;  % index of cradle data points
L = length(angleind);
cw = fdct_wrapping(zeros(opt.imgsize),opt.curveletisreal,1,L);
p = 0;
for i = 1:L
    for j = 1:length(angleind{i})
        p = p + 1;
        locsiz = size(cw{i}{angleind{i}(j)});
        tmpcw = zeros(siz);
        if flg == 1
            tmpcw(ncInd) = Y(~mask,p);
            tmpcw(cInd) = Y(mask,p);
        else
            tmpcw(cInd) = Y(:,p);
        end
%         cw{i}{angleind{i}(j)} = weightcw(tmpcw,locsiz);
                cw{i}{angleind{i}(j)} = imresize(tmpcw,locsiz,'nearest');
        
    end
end

    function cwnew = weightcw(tmpcw,locsiz)
        rowind = 1:siz(1);
        rowind = (rowind - (1+siz(1))/2)/(1+siz(1))*(1+locsiz(1))+(1+locsiz(1))/2;
        rowind = round(rowind);
        rowind = arrayfun(@(x)find(rowind == x),1:locsiz(1),'UniformOutput',0);
        cwnew = cellfun(@(x)mean(tmpcw(x,:),1),rowind,'UniformOutput',0);
        cwnew = cell2mat(cwnew(:));
        colind = 1:siz(2);
        colind = (colind - (1+siz(2))/2)/(1+siz(2))*(1+locsiz(2))+(1+locsiz(2))/2;
        colind = round(colind);
        colind = arrayfun(@(x)find(colind == x),1:locsiz(2),'UniformOutput',0);
        cwnew = cellfun(@(x)mean(cwnew(:,x),2),colind,'UniformOutput',0);
        cwnew = cell2mat(cwnew(:)');
        %     cwnew = imresize(tmpcw,locsiz,'bilinear');
    end
end
