function InputSampleNeighborPatch(ind, file_path, dir, matname,job_id,result_path)

% %% MCA decomposition
%
addpath(genpath('/home/grad/rachel/matlab/cradleremoval'))
addpath(genpath('/home/grad/rachel/matlab/SparseFactorModel'))
addpath(genpath('/home/grad/rachel/matlab/MCALabWithUtilities'))
rmpath(genpath('/home/grad/rachel/matlab/MCALabWithUtilities/CurveletToolbox'))
addpath(genpath('/home/grad/rachel/matlab/CurveLab-2.1.3/fdct_wrapping_matlab'))
addpath(genpath('/home/grad/rachel/matlab/FFST'))

if nargin > 1
    filepath = file_path;
else
    filepath = '/home/grad/rachel/Documents/MATLAB/MCAdecomp_Ghent_8/';
end

if nargin < 6
    result_path = filepath;
end

if nargin < 4
    matname = 'Ghent_8.mat';
end
% filepath= '/gtmp/DataREU/Prado/Rubens/ground_truth/simulation/';
%
%% calculate the shearlet transform on texture patches

% % ==== use low frequency feature or not ==== %
% lowfreqfeature = struct();
% lowfreqfeature.h = 1;
% lowfreqfeature.v = 0;
%
% load([filepath 'MASTER_sample.mat']);
% opt.imgsize = size(texture{1});
% ST = shearletTransformSpect(zeros(opt.imgsize));
% [j,k,cone] = arrayfun(@(x)shearletScaleShear(x), 1:size(ST,3));
%
% vcradle_sample_size = zeros(size(texture));
% hcradle_sample_size = zeros(size(texture));
% vnoncradle_sample_size = zeros(size(texture));
% hnoncradle_sample_size = zeros(size(texture));
% for i = 1:size(texture,1)
%     for jj = 1:size(texture,2)
%         img = texture{i,jj};
%         imgind = [rowind(i) (rowind(i)+size(img,1)-1);colind(jj) (colind(jj)+size(img,2)-1)];
%         opt.mask_h = info2mask(hinfo,imgind);
%         try
%             opt.mask_v = info2mask(vinfonew,imgind);
%         catch
%             opt.mask_v = info2mask(vinfo,imgind);
%         end
%         ST = shearletTransformSpect(img);
%         ST = real(ST);
%
%         opt.direction = 'vertical';
%         [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
%         opt.lowfreqfeature = lowfreqfeature.v;
%         [featureST,opt] = getSTfeature(STnew,j,opt);
%         save([filepath 'vfeatureST_' num2str(i) '_' num2str(jj) '.mat'],'featureST','opt');
%         vcradle_sample_size(i,jj) = size(featureST.cradle,1);
%         vnoncradle_sample_size(i,jj) = size(featureST.free,1);
%
%         opt.direction = 'horizontal';
%         [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
%         opt.lowfreqfeature = lowfreqfeature.h;
%         [featureST,opt] = getSTfeature(STnew,j,opt);
%         save([filepath 'hfeatureST_' num2str(i) '_' num2str(jj) '.mat'],'featureST','opt');
%         hcradle_sample_size(i,jj) = size(featureST.cradle,1);
%         hnoncradle_sample_size(i,jj) = size(featureST.free,1);
%     end
% end
% save([filepath 'ST_sample_size.mat'],'vcradle_sample_size','hcradle_sample_size','vnoncradle_sample_size','hnoncradle_sample_size');




%% build fator model by random sample on target and neighbor patches
if nargin < 3
    dir = 'horizontal';
end

try
    load([filepath matname],'subimg');
    Npatch = size(subimg);
end

if nargin > 0
    ind = str2num(ind);
else
    ind = [];
end

switch dir
    case 'horizontal'
        if isempty(ind)
            row_ind = 1:Npatch(1)
        else
            row_ind = ind;
        end
        col_ind = 1:Npatch(2);
    case 'vertical'
        row_ind = 1:Npatch(1);
        if isempty(ind)
            col_ind = 1:Npatch(2)
        else
            col_ind = ind
        end
end

num_NCsample = 500; % number of non-cradle samples draw per region
num_Csample = 1500; % number of cradle samples draw per region


r = 1; % radius of neighbor for non-cradle sampling

load([filepath matname],'*sample_size');
switch dir
    case 'horizontal'
        cradle_sample_size = hcradle_sample_size;
        noncradle_sample_size = hnoncradle_sample_size;
        STfeaturepath = [filepath 'hfeatureST_'];
    case 'vertical'
        cradle_sample_size = vcradle_sample_size;
        noncradle_sample_size = vnoncradle_sample_size;
        STfeaturepath = [filepath 'vfeatureST_'];
end

total_sample_size_unit = max(cradle_sample_size(:)+noncradle_sample_size(:));


result = cell(Npatch);
result1 = cell(Npatch);

for i = row_ind(1):row_ind(end)
    for jj = col_ind(1):col_ind(end)
        
        %**************** build sparse factor model *********************%
        if cradle_sample_size(i,jj) > 0
            %=== select non-cradle patches to sample ===%
            NCind = get2dneighbor(i,jj,Npatch(1),Npatch(2),r);
            NCind(noncradle_sample_size(NCind) < total_sample_size_unit/6) = [];
            if length(NCind) > 4
                [~,idx] = sort(noncradle_sample_size(NCind),'descend');
                NCind = NCind(idx(1:4));% for 4 patches containing most non-cradle samples are used
            elseif length(NCind) < 3
                NCind = get2dneighbor(i,jj,Npatch(1),Npatch(2),r+1);
                NCind(noncradle_sample_size(NCind) == 0) = [];
                [~,idx] = sort(noncradle_sample_size(NCind),'descend');
                NCind = NCind(idx(1:4));% for 4 patches containing most non-cradle samples are used
            end
            [row,col] = ind2sub(Npatch,NCind);
            if noncradle_sample_size(i,jj) > 0
                row = [row(:);i]; col = [col(:);jj];
            end
            
            Ync = cell(length(row),1);
            for k = 1:length(Ync)
                load([STfeaturepath num2str(row(k)) '_' num2str(col(k)) '.mat'],'featureST');
                samplesize = noncradle_sample_size(row(k),col(k));
                idx = randsample(samplesize,min(samplesize,num_NCsample));
                Ync{k} = featureST.free(idx,:);
            end
            Ync = cell2mat(Ync);
            
            %=== select cradle patches to sample ===%
            Cind = get1dneighbor(i,jj,Npatch(1),Npatch(2),dir);
            [row,col] = ind2sub(Npatch,Cind);
            Yc = cell(length(row),1);
            for k = 1:length(Yc)
                load([STfeaturepath num2str(row(k)) '_' num2str(col(k)) '.mat'],'featureST');
                samplesize = cradle_sample_size(row(k),col(k));
                idx = randsample(samplesize,min(samplesize,num_Csample));
                Yc{k} = featureST.cradle(idx,:);
            end
            Yc = cell2mat(Yc);
            
            %=== input data for factor model ===%
            Y = [Yc;Ync];
            Ymask = [ones(size(Yc,1),1);zeros(size(Ync,1),1)];
            Z = [];
            
            %=== set parameters for factor model ===%
            load([STfeaturepath num2str(i) '_' num2str(jj) '.mat'],'opt');
            disp(opt.direction)
            opt.matname = [result_path 'factormodel_' num2str(i) '_' num2str(jj) '_' opt.direction];
            opt.mrho = .8;
            opt.updaterho = 1;
            opt.method = 'supervised';
            opt.transform = 'shearlet';
            
            spfactcovest_mgploadings(Y,Ymask,opt,Z);
            save([STfeaturepath num2str(i) '_' num2str(jj) '.mat'],'opt','-append');
        end
        %
        %******************** posterior inference *************************%
        if cradle_sample_size(i,jj) > 0
            load([STfeaturepath num2str(i) '_' num2str(jj) '.mat']);
            Y = [featureST.cradle;featureST.free];
            originImg = invfeatureST(Y,opt,1);
            opt.originImg = originImg;
            [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
            try
                zout = spfactor_post_inference([opt.matname '_1.mat'],Y,opt);
                if strcmp(opt.direction,dir)
                    switch dir
                        case 'horizontal'
                            cradlemask = opt.mask_h;
                        case 'vertical'
                            cradlemask = opt.mask_v;
                    end
                    result{i,jj} = (originImg - zout.noncradleImg).*cradlemask;
                    result1{i,jj} = zout.cradleImg;
                    %                     result2{i,jj} = zout.noncradleImg.*cradlemask;
                end
            end
            % %             result{i,jj} = originImg;
        else
            load([STfeaturepath num2str(i) '_' num2str(jj) '.mat'],'opt');
            result{i,jj} = zeros(opt.imgsize);
            result1{i,jj} = zeros(opt.imgsize);
            %             result2{i,jj} = zeros(opt,imgsize);
        end
        
        
    end
end

switch dir
    case 'horizontal'
        result = result(ind(1):ind(end),:);
        result1 = result1(ind(1):ind(end),:);
        path = [result_path 'texture_separation_row' num2str(ind) '/'];
    case 'vertical'
        result = result(:,ind(1):ind(end));
        result1 = result1(:,ind(1):ind(end));
        path = [result_path 'texture_separation_col' num2str(ind) '/'];
end
if ~exist(path,'dir')
    mkdir(path);
end
% if ~exist([path 'texture.mat'], 'file')
save([path dir 'texture_separation.mat'],'result','result1');
% end
disp('result saved');

% switch dir
%     case 'horizontal'
%         hresult = result;
%     case 'vertical'
%         vresult = result;
% end



