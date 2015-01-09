function spfactor_cluster(result_path,job_id)

% load('wksp_unsupervised_model.mat','y','z','opt','mask','loc');
% load('wksp_vanilla_modle.mat','Y','cwmask','opt','loc');
addpath(genpath('/home/grad/rachel/matlab/FFST'));
path = [result_path 'result/job' job_id];
precompute = [result_path 'result/'];
% load([result_path 'result/precompute.mat']);
% FMstatus = FMstatus_vertical((3*str2num(job_id)-2):3*str2num(job_id));

colind = 9;% change the column index here!
[Y, Ymask, opt, Yind] = spfactorInputData(precompute,str2num(job_id),1,[3,45]);
    path = ['/gtmp/rachel/SparseFactorModel/19a_horizontal/job' job_id];
    if ~exist(path,'dir')
        mkdir(path);
    end
if ~isempty(Y)
    opt.method = 'supervised';
    opt.matname = [path '/output'];
    spfactcovest_mgploadings(Y,Ymask,opt);
    
    ind = sub2ind([15,9],str2num(job_id),colind); 
    [subind,jobind] = ind2sub([3,45],ind);
    load(['/home/grad/rachel/Documents/MATLAB/Preparation2/result/job' num2str(jobind) '/output.mat'],'Y','opt');
    opt = opt{subind};
    z = Y{subind};
    zout = spfactor_post_inference([path '/output_1.mat'],z,opt);
    finaloutput = struct();
    if isfield(zout,'cradleImg')
        finaloutput.cradleImg = zout.cradleImg;
        if strcmp(opt.direction,'vertical')
            finaloutput.cradleNresidualImg = (opt.originImg - zout.noncradleImg).*opt.mask_v;
        else
            finaloutput.cradleNresidualImg = (opt.originImg - zout.noncradleImg).*opt.mask_h;
        end
    else
        finaloutput = zout;
        finaloutput.originImg = opt.originImg;
    end
else
    ind = sub2ind([15,9],str2num(job_id),colind);
    [subind,jobind] = ind2sub([3,45],ind);
    load(['/home/grad/rachel/Documents/MATLAB/Preparation2/result/job' num2str(jobind) '/output.mat'],'Y','opt');
    opt = opt{subind};
    finaloutput = struct();
    finaloutput.cradleImg = zeros(opt.imgsize);
    finaloutput.cradleNresidualImg = zeros(opt.imgsize);
end
path = [path '/FinalOutput.mat'];
save(path,'-struct','finaloutput');




% savesample = savesample((3*str2num(job_id)-2):3*str2num(job_id));
% load([path '/output']);
% Opt = opt;
% for i = 1:3
%     opt = Opt{i};
%     y = Y{i};
%     pathtmp = ['/gtmp/rachel/SparseFactorModel/19a_vertical/job' job_id '/' num2str(i)];%[path '/' num2str(i)];
% if ~exist(pathtmp,'dir')
%     mkdir(pathtmp);
% end
% if strcmp(FMstatus{i},'apply')
%     mask = [ones(size(featureST{i}.cradle,1),1);zeros(size(featureST{i}.free,1),1)];
%     loc = [opt.vfeatureind.cradle;opt.vfeatureind.free];
% opt.matname = [pathtmp '/FMoutput'];
% opt.savesample = savesample(i);
%
% % ---  fixed rho value --- %
% opt.updaterho = 0;
% % opt.mrho = 0.5 + 0.1 * str2num(job_id);
% opt.updateLambda = 'single';
% opt.method = 'supervised'
%
% opt
%
% % opt.newchain = 0;
% output = spfactcovest_mgploadings(y,mask,opt,[],loc);
%
% % % analyze output
% % i = 1;
% % while exist([opt.matname '_' num2str(i) '.mat'],'file')
% %     path = [opt.matname '_' num2str(i) '.mat']
% %     factorloading2dictionary(path,opt);
% %     load(path,'*CW');
% %     tmp = struct();
% %     try
% %         tmp.noncradleImg = invfeatureCW(noncradleCW,opt);
% %         tmp.cradleImg = invfeatureCW(cradleCW,opt);
% %         tmp.residualImg = invfeatureCW(residualCW,opt);
% %         disp('scussessful!')
% %     catch
% %         disp('unscussessful!')
% %     end
% %     save(path,'-struct','tmp','-append');
% %     i = i+1;
% % end
%
% % path = [result_path 'result/job' job_id];
% % if ~exist(path,'dir')
% %     mkdir(path);
% % end
% % if ~exist([pathtmp '/wkspOutput.mat'], 'file')
% %     save([pathtmp '/wkspOutput.mat'],'-struct','output');
% % end
% end
% end
% save([result_path 'result/job' job_id '/output.mat'],'FMstatus','-append');