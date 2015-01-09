function Yimg = invfeatureCW(data,opt)

% this function inverses feature curvelet 1-by-p cell array to image domain

% Input:   data:     N-by-p cell array, each row is a feature curvelet cell
%                    array
%          opt:      struct, including info for inverse curvelet transform

if ~iscell(data)
    display('Wrong data format.')
    return
end

p = size(data,2);
N = size(data,1);

if size(data,1) > 1
    meandata = cell(1,p);
    for i = 1:p
        siz = size(data{1,i});
        tmp = cell2mat(data(:,i)');
        tmp = reshape(tmp,[],N);
        meandata{i} = mean(tmp,2);
        meandata{i} = reshape(meandata{i},siz);
    end
    data = meandata;
end

angleind = opt.featureangleind;

cw = fdct_wrapping(zeros(opt.imgsize),opt.curveletisreal,1,length(angleind));
k = 1;
for i = 1:length(angleind)
    for j = 1:length(angleind{i})
        cw{i}{angleind{i}(j)} = data{k};
        k = k + 1;
    end
end
Yimg = ifdct_wrapping(cw,opt.curveletisreal,opt.imgsize(1),opt.imgsize(2));
