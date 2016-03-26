function computeBayesianSeparation_cluster( ngrid, filepath, result_path, load_mat, job_id ,dir)

addpath(genpath('~/matlab/cradleremoval'));
addpath(genpath('~/matlab/FFST'))

% find target patch
load(load_mat,'texture','hinfo','vinfo','rowind','colind','opt');%%
npatch = length(texture(:));
patch_size = size(texture{1});
ngrid = str2double(ngrid);
job_id = str2double(job_id);
jobsize = ceil(npatch/ngrid);
jobind = jobsize*(job_id-1)+1 : min(jobsize*job_id, npatch)
Npatch = [length(rowind), length(colind)];


%% run Input SampleNeighborPatch.m
num_NCsample = 500; % number of non-cradle samples draw per region
num_Csample = 1500; % number of cradle samples draw per region

for i = 1:length(jobind)
    load([filepath dir 'featureST_' num2str(jobind(i)) '.mat'])
    disp(['load' filepath dir 'featureST_' num2str(jobind(i)) '.mat'])
    % check if there are cradle
    switch dir
        case 'h'
            if hcradle_sample_size == 0
                hresult = zeros(patch_size);
                residual = [];
                originImg = [];
            end
        case 'v'
            if vcradle_sample_size == 0
                vresult = zeros(patch_size);
                residual = [];
                originImg = [];
            end
    end
    
    % do separation if image patch contains cradled region
    if ~(exist('hresult','var') || exist('vresult','var'))
        disp('compute Bayesian separation');
        
        Ync = [];
        Yc = [];
        [r,c] = ind2sub(Npatch,jobind(i));
        
        %=== select non-cradle patches to sample ===%
        NCind = get2dneighbor(r,c,Npatch(1),Npatch(2),1);
        NCind = [NCind;jobind(i)];
        
        Ync = cell(length(NCind),1);
        for k = 1:length(Ync)
            load([filepath dir 'featureST_' num2str(NCind(k)) '.mat']);
            samplesize = size(featureST.free,1);
            if (samplesize > 2*num_NCsample) || k == length(Ync)
                idx = randsample(samplesize,min(samplesize,round(num_NCsample/4)));
                Ync{k} = featureST.free(idx,:);
            else
                Ync{k} = [];
            end
        end
        Ync = cell2mat(Ync);
        %=== select cradle patches to sample ===%
        Cind = get1dneighbor(r,c,Npatch(1),Npatch(2),opt.direction);
        Yc = cell(length(Cind),1);
        for k = 1:length(Yc)
            load([filepath dir 'featureST_' num2str(Cind(k)) '.mat']);
            samplesize = size(featureST.cradle,1);
            idx = randsample(samplesize,min(samplesize,round(num_Csample/2)));
            Yc{k} = featureST.cradle(idx,:);
        end
        Yc = cell2mat(Yc);
        
        % == run Gibbs sampling ==
        Y = [Yc; Ync];
        mask = [ones(size(Yc,1),1);zeros(size(Ync,1),1)];
        clear Yc Ync
        Z = [];
        load([filepath dir 'featureST_' num2str(jobind(i)) '.mat'],'featureST','opt');
        opt.mrho = 1;
        opt.transform = 'shearlet';
        opt.matname = [result_path 'FactorModel_' num2str(jobind(i))];
        spfactcovest_mgploadings(Y,mask,opt,Z);
%         save([filepath dir 'featureST_' num2str(jobind(i)) '.mat'],'opt','-append');
        
        % == posterior inference using the model == %
%         load([filepath dir 'featureST_' num2str(jobind(i)) '.mat'],'featureST');
        Y = [featureST.cradle;featureST.free];
        originImg = invfeatureST(Y,opt,1);
        opt.originImg = originImg;
        [Y,mask,Z,loc,opt] = PrepareData4spFactorModel(Y,opt);
        zout = spfactor_post_inference([opt.matname '_1.mat'],Y,opt);
        switch dir
            case 'h'
                hresult = (originImg -zout.noncradleImg).*opt.mask_h;
                residual = hresult-zout.cradleImg;
                hresult = zout.cradleImg;
            case 'v'
%                vresult = zout.cradleImg + (originImg - zout.noncradleImg).*opt.mask_v;
                vresult = (originImg - zout.noncradleImg).*opt.mask_v;
                residual = vresult-zout.cradleImg;
                vresult = zout.cradleImg;
        end
    end
    % save result
    switch dir
        case 'h'
            save([result_path dir 'Bayesian_result_' num2str(jobind(i)) '.mat'],'hresult','residual','originImg','opt');
            clear hresult residual
        case 'v'
            save([result_path dir 'Bayesian_result_' num2str(jobind(i)) '.mat'],'vresult','residual','originImg','opt');
            clear vresult residual
    end
end
