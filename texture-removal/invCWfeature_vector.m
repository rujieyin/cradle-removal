function curveatom = invCWfeature_vector(opt,feature)
% this function compute the inverse curvelet transform of a feature vector

cw = fdct_wrapping(zeros(opt.imgsize),opt.curveletisreal,1,7);
templatesize = ceil(opt.imgsize/2);
% templatesize = ceil(opt.imgsize./size(cw{1}))/3;
% templatesize = 2.^ceil(log2(templatesize));
cw = fdct_wrapping(zeros(templatesize),opt.curveletisreal,1,7);
basesize = size(cw{end}{opt.featureangleind{end}(1)});
loc = ceil(basesize/2);

% [X_rows, X_cols, F_rows, F_cols] = fdct_wrapping_param(cw,templatesize(1),templatesize(2));


% ---  check if the length of feature is correct --- %
if length(feature) ~= length(cell2mat(opt.featureangleind'))
    disp('The length of feature vector is wrong!')
    return
end

x_rows = zeros(length(feature),1);
x_cols = zeros(length(feature),1);

% --- assign value of curvelet coefficients ---%
l = 1;
for i = 1:length(opt.featureangleind)
    if ~isempty(opt.featureangleind{i})
        ind = opt.featureangleind{i};
        for j = 1:length(ind)
%             indlocal = round(loc./templatesize.*size(cw{i}{ind(j)}));
%             cw{i}{ind(j)}(indlocal(1),indlocal(2)) = feature(l);
            indlocal = loc./basesize.*size(cw{i}{ind(j)});
            r1 = indlocal(1)-floor(indlocal(1));
            r2 = indlocal(2)-floor(indlocal(2));
            cw{i}{ind(j)}(floor(indlocal(1)):floor(indlocal(1))+1,floor(indlocal(2)):floor(indlocal(2))+1) = feature(l)*[(1-r1)*(1-r2),(1-r1)*r2;r1*(1-r2),r1*r2];
            l = l + 1;
        end
    end
end

curveatom = ifdct_wrapping(cw,opt.curveletisreal,templatesize(1),templatesize(2));