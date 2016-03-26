function computeShearletDescriptor_cluster( ngrid, result_path, load_mat, job_id )

addpath(genpath('~/matlab/cradleremoval'));
addpath(genpath('~/matlab/FFST'))

% load target texture patch
load(load_mat,'texture','hinfo','vinfo','rowind','colind','opt');%%
npatch = length(texture(:));
ngrid = str2double(ngrid);
job_id = str2double(job_id);
jobsize = ceil(npatch/ngrid);
jobind = jobsize*(job_id-1)+1 : min(jobsize*job_id, npatch) 
texture = texture(jobind);


for i = 1:length(jobind)
        img = texture{i};
        [r,c] = ind2sub([length(rowind),length(colind)],jobind(i));
        imgind = [rowind(r) (rowind(r)+size(img,1)-1);colind(c) (colind(c)+size(img,2)-1)];
        opt.mask_h = info2mask(hinfo,imgind);
        opt.mask_v = info2mask(vinfo,imgind);
        ST = shearletTransformSpect(img);
        ST = real(ST);
        [j,k,cone] = arrayfun(@(x)shearletScaleShear(x), 1:size(ST,3));
        opt.imgsize = size(img);
        optcpy = opt;% make a copy of opt
        
        opt.direction = 'vertical';
        opt.lowfreqfeature = 1;
        [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
        [featureST,opt] = getSTfeature(STnew,j,opt);
        vcradle_sample_size = size(featureST.cradle,1);
        vnoncradle_sample_size = size(featureST.free,1);
        save([result_path 'vfeatureST_' num2str(jobind(i)) '.mat'],'featureST',...
            'opt','vcradle_sample_size','vnoncradle_sample_size');

        opt = optcpy;
        opt.direction = 'horizontal';
        opt.lowfreqfeature = 1; % whether to include low frequency feature or not
        [STnew ,opt] = shearlet_angleselection(ST,j,k,cone,opt);
        [featureST,opt] = getSTfeature(STnew,j,opt);
        hcradle_sample_size = size(featureST.cradle,1);
        hnoncradle_sample_size = size(featureST.free,1);
%         if sum(opt.mask_h(:)) < length(img(:))*.1
%             featureST = [];
%             hresult{i,jj} = zeros(size(img));
%             continue;
%         end
        save([result_path 'hfeatureST_' num2str(jobind(i)) '.mat'],'featureST',...
            'opt','hcradle_sample_size','hnoncradle_sample_size');
%         
end



