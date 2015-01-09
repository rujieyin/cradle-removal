function zout = spfactor_post_inference(output, z, opt, NCp)
% this function does posterior inference of z using sparse factor model
% MCMC output

% for data z from different subpatch from that of output, opt should be the
% one for z's subpatch

nz = size(z,1);
p = size(z,2);
if nargin < 4
    if ~isfield(opt,'method') || strcmp(opt.method, 'supervised')
        switch opt.direction
            case 'vertical'
        NCp = [zeros(size(opt.vfeatureind.cradle,1),1);ones(size(opt.vfeatureind.free,1),1)];
            case 'horizontal'
                NCp = [zeros(size(opt.hfeatureind.cradle,1),1);ones(size(opt.hfeatureind.free,1),1)];
        end
    else
    NCp = ones(size(z,1),1)*.5;
    end
end
load(output,'Lambda','Gamma','ps','rho','kappa','M','VY');
% -- renormalize z -- %
z = bsxfun(@minus,z,M);
z = bsxfun(@times,z,1./sqrt(VY));
zout = struct();
% -- compute conditional prob of x -- %
px1 = cellfun(@noncradle_cond,Lambda,ps,'UniformOutput',0);
px1 = cell2mat(px1')';
px2 = cellfun(@cradle_cond,Lambda,num2cell(rho),Gamma,kappa,ps,'UniformOutput',0);
px2 = cell2mat(px2')';
px1 = bsxfun(@times,px1,NCp(:)');
px2 = bsxfun(@times,px2,(1-NCp(:)'));
px = px1./(px1+px2);
% -- compute conditional prob of z_nc, z_c, z_e -- %
nsample = size(px,1);
k1 = size(Lambda{1},2);
k2 = size(Gamma{1},2);
zout.x = logical(binornd(1,px,nsample,nz));        % nsample-by-nz
zout.eta = cell(nsample,1);              % nz - by - k1
zout.xi = cell(nsample,1);               % nz - by - k2
zout.noncradle = zeros(nsample,nz*p);
zout.cradle = zeros(nsample,nz*p);
for i = 1:nsample
    x = zout.x(i,:);
    Nnoncradle = sum(x);
    Ncradle = sum(~x);
    % -- sample eta -- %
    eta = zeros(nz,k1);
    % non-cradle signal
    Lmsg = bsxfun(@times,Lambda{i},ps{i});
    Veta1 = eye(k1) + Lmsg'*Lambda{i};
    T = cholcov(Veta1); [~,sR] = qr(T);
    S = inv(sR); Veta = S*S';                   % Veta = inv(Veta1)
    Meta = z(x,:)*Lmsg*Veta;
    eta(x,:) = Meta + normrnd(0,1,[Nnoncradle,k1])*S';
    % cradle signal
    Lmsg = inv(Gamma{i}*Gamma{i}'+diag(1./ps{i}))*Lambda{i};
    Veta1 = eye(k1) + rho(i)^2*Lmsg'*Lambda{i};
    T = cholcov(Veta1); [~,sR] = qr(T);
    S = inv(sR); Veta = S*S';                   % Veta = inv(Veta1)
    Meta = z(~x,:)*Lmsg*Veta;
    eta(~x,:) = Meta + normrnd(0,1,[Ncradle,k1])*S';
    zout.eta{i} = eta;
    % -- compute z_nc -- %
    z_nc = zeros(nz,p);
    z_nc(x,:) = eta(x,:)*Lambda{i}';
    z_nc(~x,:) = eta(~x,:)*Lambda{i}'*rho(i);
    zout.noncradle(i,:) = z_nc(:)';
    % -- sample xi -- %
    xi = zeros(nz,k2);
    Gmsg = bsxfun(@times,Gamma{i},ps{i});
    Vxi1 = eye(k2) + Gmsg'*Gamma{i};
    T = cholcov(Vxi1); [~,sR] = qr(T);
    S = inv(sR); Vxi = S*S';
    Mxi = bsxfun(@plus,(z(~x,:) - z_nc(~x,:))*Gmsg*Vxi,kappa{i}*Vxi);
    xi(~x,:) = Mxi + normrnd(0,1,[Ncradle,k2])*S';
    zout.xi{i} = xi;
    % -- compute z_c -- %
    z_c = zeros(nz,p);
    z_c(~x,:) = xi(~x,:)*Gamma{i}';
    zout.cradle(i,:) = z_c(:)';
end

zout.noncradle = mean(zout.noncradle,1);
zout.noncradle = reshape(zout.noncradle,[nz,p]);
zout.cradle = mean(zout.cradle,1);
zout.cradle = reshape(zout.cradle,[nz,p]);
zout.residual = z - zout.cradle - zout.noncradle;
noncradle = bsxfun(@times,zout.noncradle,sqrt(VY));
noncradle = bsxfun(@plus,noncradle,M);
cradle = bsxfun(@times,zout.cradle,sqrt(VY));
residual = bsxfun(@times,zout.residual,sqrt(VY));
switch opt.transform
    case 'curvelet'
cwsiz = compute_original_cw_size;
zout.noncradleCW = collapse_cwfeature(noncradle);
zout.cradleCW = collapse_cwfeature(cradle);
zout.residualCW = collapse_cwfeature(residual);
zout.noncradleImg = invfeatureCW(zout.noncradleCW,opt);
zout.cradleImg = invfeatureCW(zout.cradleCW,opt);
zout.residualImg = invfeatureCW(zout.residualCW,opt);
    case 'shearlet'
        zout.noncradleST = noncradle;
        zout.cradleST = cradle;
        zout.noncradleImg = invfeatureST(noncradle,opt,1);
        zout.cradleImg = invfeatureST(cradle,opt,1);                   
        try
        zout.residualImg = opt.originImg - zout.noncradleImg - zout.cradleImg;
        end
end

px = mean(px);
zout.px = px;
zout.px1 = px1;
zout.px2 = px2;


    function prob = noncradle_cond(Lambda,ps)
        prob = mvnpdf(z,zeros(1,p),Lambda*Lambda'+diag(1./ps));
    end

    function prob = cradle_cond(Lambda,rho,Gamma,kappa,ps);
        prob = mvnpdf(z,kappa*Gamma',Gamma*Gamma' + rho^2*Lambda*Lambda' + diag(1./ps));
    end

    function CW = collapse_cwfeature(Z,locind)
        if nargin < 2
            locind = opt.unsupervised_setup.IndZ;
        end
        CW = cell(1,p);
        for ii = 1:p
            cwtmp = zeros(opt.vdownsamplesize);
            cwtmp(locind) = Z(:,ii);
            CW{ii} = imresize(cwtmp,cwsiz{ii},'bilinear');
        end
    end

    function siz = compute_original_cw_size
        angleind = opt.featureangleind;
        L = length(angleind);
        cw = fdct_wrapping(zeros(opt.imgsize),opt.curveletisreal,1,L);
        siz = cell(p,1);
        k = 1;
        for ii = 1:L
            for jj = 1:length(angleind{ii})
                siz{k} = size(cw{ii}{angleind{ii}(jj)});
                k = k + 1;
            end
        end
    end


end