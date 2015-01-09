% -- Rachel Yin -- %

% Gibbs sampler for covariance estimation
% using mgps prior on factor loadings


function output = spfactcovest_mgploadings(Y,mask,opt,Z,loc)
% clear;clc;close all;
tic;

% --- choose type of methods --- %
if isempty(opt)
    opt = struct();
end

if isfield(opt,'method')
    method = opt.method;
else
    if nargin >= 4 && ~isempty(Z)
        method = 'unsupervised';
    else
        if sum(mask) > 0
            method = 'supervised';
        else
            method = 'vanilla';
        end
    end
end


% --- start a new chain or not --- %
if ~isfield(opt,'newchain')
    opt.newchain = 1;
end


% --- check input data --- %
if ~isempty(mask)
    mask = mask>0;
    Ncradle = sum(mask);
    Nnoncradle = length(mask) - Ncradle;
    if sum(~mask) == 0
        disp('Should include at least one non-cradling data point.')
        return
    end
end

% --- renormalize data to get zero mean of non-cradle component --- %

switch method
    case 'vanilla'
        M = mean(Y);
        VY = var(Y);
        Y = bsxfun(@minus,Y,M);                     % center the training data
        Y = bsxfun(@times,Y,1./sqrt(VY));           % scale the training data
    otherwise
        M = mean(Y(~mask,:));                       % make non-cradling part mean-zero
        VY= var(Y(~mask,:));
        Y = bsxfun(@minus,Y,M);                     % center the training data
        Y = bsxfun(@times,Y,1./sqrt(VY));           % scale the training data
        if nargin > 3 && ~isempty(Z)
            Z = bsxfun(@minus,Z,M);
            Z = bsxfun(@times,Z,1./sqrt(VY));
        end
end

if strcmp(opt.transform, 'curvelet')
    % --- compute original size of curvelet matrices --- %
    cwsiz = compute_original_cw_size;
end

% --- get location index of data in the image domain --- %
if nargin == 5
    if strcmp(method,'supervised')
        Loc = loc.Y;
        Z = [];
    else
        Loc = [loc.Y(:);loc.Z(:)];
    end
%     Loc = loc;
else
    Loc = [];
end

% --- compute origin image if not provided -- %
% if isfield(opt,'originImg')
%     originImg = opt.originImg;
% else
%     switch opt.transform
%         case 'curvelet'
%             originImg = collapse_CWfeature_vector([Y;Z],mask,opt);
%             originImg = ifdct_wrapping(originImg,opt.curveletisreal,opt.imgsize(1),opt.imgsize(2));
%         case 'shearlet'
%             originImg = invfeatureST([Y;Z],opt,1,Loc)
%     end
%     opt.originImg = originImg;
% end

% --- define global constants --- %
if opt.newchain
    nrun = 2400; burn = 2000; thin = 1; sp = (nrun - burn)/thin; % number of posterior samples
else
    nrun = 5000; burn = 0; thin = 1; sp = (nrun-burn)/thin;
end
nchunk = 400; % size of chunk saved to local file
if nrun - burn > nchunk
    nsave = nchunk;
else
    nsave = nrun - burn;
end

p = size(Y,2); rep = 1; n = size(Y,1);
if nargin > 3
    nz = size(Z,1);
end

kinit = repmat(floor(log(p)*3),rep,1);                 % number of factors to start with
b0 = 1; b1 = 0.0005;
epsilon = 1e-2;%  1e-1;                                      % threshold limit
prop = .95;                                           % proportion of redundant elements within columns


%---- define output files across replicates-----%

% --- .mat file for mcmc draws --- %
if isfield(opt,'matname')
    path = opt.matname;
else
    path = 'output';
end
name = 'spfactcovest_mgploadings';

% ---- save samples or output Img's only --- %
if ~isfield(opt,'savesample')
    opt.savesample = 1;
end


output = struct();
output.method = method;
output.nrun = nrun;
output.burn = burn;
output.nchunk = nchunk;
output.newchain = 0;
output.transform = opt.transform;
output.matname = opt.matname;

% %% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% save(path,'-v7.3','name');
% x = cell(nrun,1);
% save(path,'x','-append');
% matobj = matfile(path,'Writable',logical(1));
%%
% if nsave < nrun - burn
%     Lambda = cell(nrun-burn,1);save(path,'Lambda','-append');clear Lambda;
%     eta = cell(nrun-burn,1);save(path,'eta','-append');clear eta;
%     ps = cell(nrun-burn,1);save(path,'ps','-append');clear ps;
%     if strcmp(method,'vanilla')
%         betac = cell(nrun-burn,1);save(path,'betac','-append');clear betac;
%         betan = cell(nrun-burn,1);save(path,'betan','-append');clear betan;
%         factorCW = cell(nrun-burn,p);save(path,'factorCW','-append');clear factorCW;
%     else
%         if strcmp(method,'unsupervised')
%             beta = cell(nrun-burn,1);save(path,'beta','-append');clear beta;
%         end
%         xi = cell(nrun-burn,1);save(path,'xi','-append');clear xi;
%         kappa = cell(nrun-burn,1);save(path,'kappa','-append');clear kappa;
%         Gamma = cell(nrun-burn,1);save(path,'Gamma','-append');clear Gamma;
%         rho = zeros(nrun-burn,1);save(path,'rho','-append');clear rho;
%         noncradleCW = cell(nrun-burn,p);save(path,'noncradleCW','-append');clear noncradleCW;  % noncradle component curvelet coefficients
%         cradleCW = cell(nrun-burn,p);save(path,'cradleCW','-append');clear cradleCW;           % cradle component curvelet coefficients
%     end
%     residualCW = cell(nrun-burn,p);save(path,'residualCW','-append');clear residualCW
% end

output.Lambda = cell(nsave,1);
output.ps = cell(nsave,1);
output.eta = cell(nsave,1);
if strcmp(method,'vanilla')
    if strcmp(opt.transform,'curvelet')
        output.factorCW = cell(nsave,p);
    elseif strcmp(opt.transform,'shearlet')
        output.factorST = cell(nsave,1);
    end
    output.betac = cell(nsave,1);
    output.betan = cell(nsave,1);
else
    if strcmp(method,'unsupervised')
        output.etaz = cell(nsave,1);                 % eta value of data in Z
        output.xiz = cell(nsave,1);                  % xi value of data in Z
        output.x = cell(nsave,1);                    % indicator of noncradle/cradle
        %         output.beta = cell(nsave,1);
    end
    output.xi = cell(nsave,1);
    output.kappa = cell(nsave,1);
    output.Gamma = cell(nsave,1);
    %         output.Vxi = cell(nrun-burn,1);
    output.rho = zeros(nsave,1);
    if strcmp(opt.transform,'curvelet')
        output.noncradleCW = cell(nsave,p);
        output.cradleCW = cell(nsave,p);
    elseif strcmp(opt.transform,'shearlet')
        output.noncradleST = cell(nsave,1);
        output.cradleST = cell(nsave,1);
    end
end
%     output.Veta = cell(nrun-burn,1);




for g = 1:rep
    
    disp(['start replicate','',num2str(g)]);
    disp('--------------------');
    
    % ------read data--------%
    
    num = 0;                                    % no. of redundant columns
    %     k=kinit(g);                                 % no. of factors to start with
    k1 = size(Y,2);
    k2 = floor(log(p)*3);
    
    % --- Define hyperparameter values --- %
    
    as = 1;bs = 0.3;                           % gamma hyperparameters for residual precision
    df = 3;                                    % gamma hyperparameters for t_{ij}
    ad1 = 2.1;bd1 = 1;                         % gamma hyperparameters for delta_1
    ad2 = 3.1;bd2 = 1;                         % gamma hyperparameters delta_h, h >= 2
    Ad1 = ad1;
    Ad2 = ad2;
    adf = 1; bdf = 1;                          % gamma hyperparameters for ad1 and ad2 or df
    if ~strcmp(method,'vanilla')
        mrho = .8; vrho = .3^2; rho_lowerbound = 0;                                  % mean and variance of truncated normal of rho
        output.mrho = mrho; output.vrho = vrho; output.rho_lowerbound = rho_lowerbound;
    end
    
    
    % --- Initial values --- %
    if opt.newchain
        ps = 10*ones(p,1);%%%*****ps=gamrnd(as,1/bs,p,1); 
        Sigma=diag(1./ps);                  % Sigma = diagonal residual covariance
        ps2 = ps;%%%******ps2 = gamrnd(as,1/bs,p,1);
        Sigma2 = diag(1./ps2);
        switch method
            case 'vanilla'
                if sum(mask) > 0
                    betac = normrnd(0,1,[1,k1]);                             % mean of eta for cradle
                else
                    betac = zeros(1,k1);
                end
                betan = normrnd(0,1,[1,k1]);                             % mean of eta for noncradle
            case 'unsupervised'
                NCp = .5;                                                % bernoulli parameter for x's
                ap1 = round(nz/2); ap2 = round(nz/2);                    % beta hyperparameters for NCp
                output.ap1 = ap1; output.ap2 = ap2;
                x = logical(binornd(1,NCp,nz,1)); x = x(:);
                etaz = normrnd(0,1,[nz,k1]);
                %                 beta = normrnd(0,1,[1,k1]);                              % mean of eta
        end
        Lambda = zeros(p,k1); % eta =  normrnd(0,1,[n,k1]);             % factor loadings & latent factors
        %     meta = zeros(n,k1); veta = eye(k1);                           % latent factor distribution = standrad normal
        
        psijh1 = gamrnd(df/2,2/df,[p,k1]);                            % local shrinkage coefficients
        delta1 = ...
            [gamrnd(ad1,bd1);gamrnd(ad2,bd2,[k1-1,1])];              % gobal shrinkage coefficients multilpliers
        tauh1 = cumprod(delta1);                                      % global shrinkage coefficients
        Plam = bsxfun(@times,psijh1,tauh1');                          % precision of loadings rows
        if ~strcmp(method,'vanilla')
            if isfield(opt,'mrho')
                rho = opt.mrho;
            else
                rho = mrho;
            end
            Gamma = zeros(p,k2);                                     % cradle factor loadings
            psijh2 = gamrnd(df/2,2/df,[p,k2]);
            delta2 = ...
                [gamrnd(Ad1,bd1);gamrnd(Ad2,bd2,[k2-1,1])];
            tauh2 = cumprod(delta2);
            Pgam = bsxfun(@times,psijh2,tauh2');
            kappa = normrnd(0,1,[1,k2]);                             % mean of kappa
            xi = bsxfun(@plus, normrnd(0,1,[Ncradle,k2]),kappa);
            if strcmp(method,'unsupervised')
                xiz = zeros(nz,k2);
                xiz(~x,:) = bsxfun(@plus,normrnd(0,1,[nz-sum(x),k2]),kappa);
            end
        end
    else
        k1 = size(opt.Lambda{end},2);                            % number of factors
        ps = opt.ps{end};
        Sigma = diag(1./opt.ps{end});
        if strcmp(method,'vanilla')
            betac = opt.betac{end};
            betan = opt.betan{end};
        else
            if strcmp(method,'unsupervised')
                %                 beta = opt.beta{end};
                x = opt.x{end};
                etaz = opt.etaz{end};
                xiz = opt.xiz{end};
                ap1 = opt.ap1;
                ap2 = opt.ap2;
            end
            k2 = size(opt.Gamma{end},2);
            kappa = opt.kappa{end};
            Gamma = opt.Gamma{end};
            psijh2 = opt.psijh2;
            delta2 = opt.delta2;
            tauh2 = opt.tauh2;
            Pgam = opt.Pgam;
            xi = opt.xi{end};
            rho = opt.rho(end);
            mrho = opt.mrho;
            vrho = opt.vrho;
            rho_lowerbound = opt.rho_lowerbound;
        end
        Lambda = opt.Lambda{end};
        psijh1 = opt.psijh1;
        delta1 = opt.delta1;
        tauh1 = opt.tauh1;
        Plam = opt.Plam;
    end
    
    
    % --- Define output files specific to replicate --- %
    nofout1 = zeros(nrun+1,1);                  % number of factors across iteartions
    nofout1(1) = k1;
    nofout2 = zeros(nrun+1,1);
    nofout2(1) = k2;
    nof1out = zeros(sp,1);
    Omegaout = zeros(p^2,1);
    Omega1out = zeros(p^2,1);
    
    
    
    
    %%
    %------start gibbs sampling-----%
    
    for i = 1:nrun
        
        % -- Update eta, etaz -- %
        Lmsg = bsxfun(@times,Lambda,ps);
        Veta1 = eye(k1) + Lmsg'*Lambda;
        T = cholcov(Veta1); [~,sR] = qr(T);
        S = inv(sR); Veta = S*S';                   % Veta = inv(Veta1)
        switch method
            case {'supervised','unsupervised'}
                Meta = Y(~mask,:)*Lmsg*Veta;
                eta(~mask,:) = Meta + normrnd(0,1,[Nnoncradle,k1])*S';
                if strcmp(method,'unsupervised')
                    Meta = Z(x,:)*Lmsg*Veta;
                    etaz(x,:) = Meta + normrnd(0,1,[sum(x),k1])*S';
                end
                Veta1 = eye(k1) + rho^2*Lmsg'*Lambda;
                T = cholcov(Veta1); [~,sR] = qr(T);
                S = inv(sR); Veta = S*S';
                Meta = (Y(mask,:) - xi*Gamma')*rho*Lmsg*Veta;
                eta(mask,:) = Meta + normrnd(0,1,[Ncradle,k1])*S';
                if strcmp(method,'unsupervised')
                    Meta = (Z(~x,:) - xiz(~x,:)*Gamma')*rho*Lmsg*Veta;
                    etaz(~x,:) = Meta + normrnd(0,1,[nz-sum(x),k1])*S';
                end
            case 'vanilla'
                Meta = (Y*Lmsg + mask*betac + ~mask*betan)*Veta;   % n x k1
                eta = Meta + normrnd(0,1,[n,k1])*S';        % update eta in a block
        end
        
        
        
        % -- Update beta's -- %
        switch method
            case 'vanilla'
                if sum(mask) > 0
                    betac = arrayfun(@(x)normrnd(x,1/Ncradle,[1,1]),sum(eta(mask,:),1)/Ncradle);
                end
                betan = arrayfun(@(x)normrnd(x,1/Nnoncradle,[1,1]),sum(eta(~mask,:),1)/Nnoncradle);
                %             case 'unsupervised'
                %                 beta = arrayfun(@(x)normrnd(x,1/n,[1,1]),sum(eta,1)/n);
        end
        
        
        % -- update Lambda (rue & held) -- %
        switch method
            case 'vanilla'
                eta2 = eta'*eta;
                alam = eta'*Y;
            case 'supervised'
                eta2 = eta(~mask,:)'*eta(~mask,:);
                alam = eta(~mask,:)'*Y(~mask,:);
                if isfield(opt,'updateLambda') && strcmp(opt.updateLambda,'both')
                    eta2 = eta2 + eta(mask,:)'*eta(mask,:)*rho^2;
                    alam = alam + rho*eta(mask,:)'*(Y(mask,:) - xi*Gamma');
                end
            case 'unsupervised'
                eta2 = eta(~mask,:)'*eta(~mask,:) + etaz(x,:)'*etaz(x,:);
                alam = eta(~mask,:)'*Y(~mask,:) + etaz(x,:)'*Z(x,:);
        end
        for j = 1:p
            Qlam = diag(Plam(j,:)) + ps(j)*eta2; blam = ps(j)*alam(:,j);
            Llam = chol(Qlam,'lower'); zlam = normrnd(0,1,k1,1);
            vlam = Llam\blam; mlam = Llam'\vlam; ylam = Llam'\zlam;
            Lambda(j,:) = (ylam + mlam)';
        end
        
        
        %------Update psi_{jh}'s------%
        psijh1 = gamrnd(df/2 + 0.5,1./(df/2 + bsxfun(@times,Lambda.^2,tauh1')));
        
        %------Update delta & tauh------%
        mat = bsxfun(@times,psijh1,Lambda.^2);
        ad = ad1 + 0.5*p*k1; bd = bd1 + 0.5*(1/delta1(1))*sum(tauh1.*sum(mat)');
        delta1(1) = gamrnd(ad,1/bd);
        tauh1 = cumprod(delta1);
        
        for h = 2:k1
            ad = ad2 + 0.5*p*(k1-h+1); bd = bd2 + 0.5*(1/delta1(h))*sum(tauh1(h:end).*sum(mat(:,h:end))');
            delta1(h) = gamrnd(ad,1/bd); tauh1 = cumprod(delta1);
        end
        
        
        if ~strcmp(method,'vanilla')
            
            % -- Update xi, xiz -- %
            Gmsg = bsxfun(@times,Gamma,ps);
            Vxi1 = eye(k2) + Gmsg'*Gamma;
            T = cholcov(Vxi1); [~,sR] = qr(T);
            S = inv(sR); Vxi = S*S';
            Mxi = bsxfun(@plus , (Y(mask,:) - rho*eta(mask,:)*Lambda')*Gmsg*Vxi , kappa*Vxi);
            xi = Mxi + normrnd(0,1,[Ncradle,k2])*S';
            if strcmp(method,'unsupervised')
                Mxi = bsxfun(@plus,(Z(~x,:) - rho*etaz(~x,:)*Lambda')*Gmsg*Vxi,kappa*Vxi);
                xiz(~x,:) = Mxi + normrnd(0,1,[nz - sum(x),k2])*S';
            end
            
            % -- Update kappa --%
            switch method
                case 'supervised'
                    kappa = arrayfun(@(x)normrnd(x,1/Ncradle,[1,1]),sum(xi,1)/Ncradle);
                case 'unsupervised'
                    kappa = arrayfun(@(x)normrnd(x,1/(Ncradle+nz-sum(x)),[1,1]),sum([xi;xiz(~x,:)],1)/(Ncradle+nz-sum(x)));
            end
            if sum(isnan(kappa))
                dbstop
            end
            
            % -- update Gamma -- %
            xi2 = xi'*xi;
            agam = xi'*(Y(mask,:)-rho*eta(mask,:)*Lambda');
            if strcmp(method,'unsupervised')
                xi2 = xi2 + xiz(~x,:)'*xiz(~x,:);
                agam = agam + xiz(~x,:)'*(Z(~x,:)-rho*etaz(~x,:)*Lambda');
            end
            for j = 1:p
                Qgam = diag(Pgam(j,:)) + ps2(j)*xi2;
                bgam = ps2(j)*agam(:,j);% + (Pgam(j,:).*kappa)';
                Lgam = chol(Qgam,'lower'); zgam = normrnd(0,1,k2,1);
                vgam = Lgam\bgam; mgam = Lgam'\vgam; ygam = Lgam'\zgam;
                Gamma(j,:) = (ygam + mgam)';
            end
            
            %-----Update psi2_{jh}'s ------%
            psijh2 = gamrnd(df/2 + .5, 1./(df/2 + bsxfun(@times,Gamma.^2,tauh2')));
            
            %-----Update delta2 & tauh2 ----%
            mat = bsxfun(@times,psijh2,Gamma.^2);
            ad = Ad1 + .5*p*k2; bd = bd1 + .5*(1/delta2(1))*sum(tauh2.*sum(mat)');
            delta2(1) = gamrnd(ad,1/bd);
            tauh2 = cumprod(delta2);
            for h = 2:k2
                ad = Ad2 + 0.5*p*(k2-h+1); bd = bd2 + 0.5*(1/delta2(h))*sum(tauh2(h:end).*sum(mat(:,h:end))');
                delta2(h) = gamrnd(ad,1/bd); tauh2 = cumprod(delta2);
            end
            
            Pgam = bsxfun(@times,psijh2,tauh2');
            
        
        if strcmp(method,'unsupervised')
            % -- update NCp -- %
            NCp = betarnd(ap1+sum(x),ap2+nz-sum(x));
            
            % -- update x -- %
            zp = mvnpdf(Z,zeros(1,p),Lambda*Lambda'+diag(1./ps));
            zp = NCp*zp./(NCp*zp + (1-NCp)*mvnpdf(Z,kappa*Gamma',rho^2*Lambda*Lambda'+Gamma*Gamma'+diag(1./ps)));
            xnew = binornd(1,zp); xnew = logical(xnew(:));
            xiz(xnew & ~x,:) = 0;
            xiz(~xnew & x,:) = normrnd(0,1,[sum(~xnew & x),k2]);
            x = xnew;
            
        end
        % -- Update Sigma -- %
        switch method
            case 'vanilla'
                Y_f = eta*Lambda';
                Ytil = Y - Y_f;
            case 'supervised'
                %=== update using non-cradle data only ===%
                Ytil = zeros(sum(~mask),size(Y,2));
                Ytil = Y(~mask,:) - eta(~mask,:)*Lambda';
                Y_nc = zeros(size(Y));
                Y_nc(mask,:) = rho*eta(mask,:)*Lambda';
                Y_nc(~mask,:) = eta(~mask,:)*Lambda';
                Y_c = xi*Gamma';
                %=== set different residual of cradle data ===%
                Ytil2 = zeros(sum(mask),size(Y,2));
                Ytil2 = Y(mask,:) - Y_nc(mask,:) - Y_c;
%                 % === update using full data ===%
%                 Ytil = zeros(size(Y));
%                 Ytil(mask,:) = rho*eta(mask,:)*Lambda';
%                 Ytil(~mask,:) = eta(~mask,:)*Lambda';
%                 Y_nc = Ytil;
%                 Y_c = xi*Gamma';
%                 Ytil = Y - Y_nc;
%                 Ytil(mask,:) = Ytil(mask,:) - Y_c;
            case 'unsupervised'
                Ytil = zeros(n+nz,p);
                Ytil([mask;~x],:) = rho*[eta(mask,:);etaz(~x,:)]*Lambda';
                Ytil([~mask;x],:) = [eta(~mask,:);etaz(x,:)]*Lambda';
                Y_nc = Ytil;
                Y_c = [xi;xiz(~x,:)]*Gamma';
                Ytil = [Y;Z] - Y_nc;
                Ytil([mask;~x],:) = Ytil([mask;~x],:) - Y_c;
                %                 Ytil(mask,:) = Y(mask,:) - rho*eta(mask,:)*Lambda' - xi*Gamma';
                %                 Ytil(~mask,:) = Y(~mask,:) - eta(~mask,:)*Lambda';
        end
        %         Yres = sum(Ytil.^2,2);
        
        % ===DONOT update the hyper-parameter for residual (ps, Sigma) === %
% %         ps=gamrnd(as + 0.5*size(Ytil,1),1./(bs+0.5*sum(Ytil.^2)))';
% %         try
% %             ps2 = gamrnd(as + .5*size(Ytil2,1),1./(bs+.5*sum(Ytil2.^2)))';
% %         catch
% %              ps2 = ps;
% %         end
% %         Sigma=diag(1./ps);
        
        
                    %-----Update rho ----%
            if isfield(opt, 'updaterho') && opt.updaterho
                if isfield(opt,'rhoUpdateIndex')
                    UpdateIndex = opt.rhoUpdateIndex;
                else
                    UpdateIndex = 2/3;
                end
                if i >  burn*UpdateIndex
                    switch method
                        case 'supervised'
                            inVrho = inv(Gamma*Gamma'+diag(1./ps2));
                            vrho1 = eta(mask,:)*Lambda'*inVrho*Lambda.*eta(mask,:);
                            mrho1 = bsxfun(@minus,Y(mask,:), kappa*Gamma')*inVrho*Lambda.*eta(mask,:);
%                             vrho1 = eta(mask,:)*Lambda'*diag(ps)*Lambda.*eta(mask,:);
%                             mrho1 = (Y(mask,:) - xi*Gamma')*diag(ps)*Lambda.*eta(mask,:);
                        case 'unsupervised'
                            vrho1 = [eta(mask,:);etaz(~x,:)]*Lambda'*diag(ps)*Lambda.*[eta(mask,:);etaz(~x,:)];
                            mrho1 = ([Y(mask,:);Z(~x,:)] - [xi;xiz(~x,:)]*Gamma')*diag(ps)*Lambda.*[eta(mask,:);etaz(~x,:)];
                    end
                    vrho1 = sum(vrho1(:)) + 1/vrho;%vrho1 = sum(vrho1(:).^2) + 1/vrho;
                    mrho1 = sum(mrho1(:)) + mrho/vrho;%mrho1 = sum(mrho1(:).^2) + mrho/vrho;
                    mrhonew = mrho1/vrho1;
                    vrhonew = 1/vrho1;
                    lcdf = normcdf(rho_lowerbound,mrhonew,sqrt(vrhonew));
                    rcdf = normcdf(1,mrhonew,sqrt(vrhonew));
                    rho1 = rand(1);
                    rho = norminv((rcdf-lcdf)*rho1+lcdf,mrhonew,sqrt(vrhonew));              % truncated normal distribution
                    if isinf(rho)
                        disp(mrhonew);
                        if mrhonew < rho_lowerbound
                            rho = rho_lowerbound
                        else
                            rho = 1;
                        end
                        %                         disp(['Itr ',num2str(i),'rho is Inf, reset to 1.']);
                    end
                    disp(['Itr ', num2str(i), 'mrho :',num2str(mrhonew),'vrho :',num2str(vrhonew),' rho is updated to:', num2str(rho)]);
                    %         rho = mrho + normrnd(0,1,[1,1])*sqrt(vrho);
                end
            end
            
        end

        %---update precision parameters----%
        Plam = bsxfun(@times,psijh1,tauh1');
        
        
        if i < burn%opt.newchain
            % ----- make adaptations for non-cradle parameters ----%
            prob = 1/exp(b0 + b1*i);                % probability of adapting
            uu = rand;
            lind = sum(abs(Lambda) < epsilon)/p;    % proportion of elements in each column less than eps in magnitude
            vec = lind >=prop;num = sum(vec);       % number of redundant columns
            
            if uu < prob
                if  i > 20 && num == 0 && all(lind < 0.995)
                    k1 = k1 + 1;
                    Lambda(:,k1) = zeros(p,1);
                    
                    if strcmp(method,'vanilla')
                        if sum(mask) > 0
                            betac(k1) = normrnd(0,1,[1,1]);
                        else
                            betac(k1) = 0;
                        end
                        betan(k1) = normrnd(0,1,[1,1]);
                        eta(mask,k1) = normrnd(betac(k1),1,[Ncradle,1]);
                        eta(~mask,k1) = normrnd(betan(k1),1,[Nnoncradle,1]);
                    else
                        eta(:,k1) = normrnd(0,1,[n,1]);
                        if strcmp(method,'unsupervised')
                            etaz(:,k1) = normrnd(0,1,[nz,1]);
                        end
                    end
                    psijh1(:,k1) = gamrnd(df/2,2/df,[p,1]);
                    delta1 = [delta1;gamrnd(ad2,1/bd2)];
%                     delta1(k1) = gamrnd(ad2,1/bd2);
                    tauh1 = cumprod(delta1);
                    Plam = bsxfun(@times,psijh1,tauh1');
                elseif num > 0 && num < k1
                    nonred = setdiff(1:k1,find(vec)); % non-redundant loadings columns
                    k1 = max(k1 - num,1);
                    Lambda = Lambda(:,nonred);
                    psijh1 = psijh1(:,nonred);
                    eta = eta(:,nonred);
                    switch method
                        case 'vanilla'
                            betac = betac(nonred);
                            betan = betan(nonred);
                        case 'unsupervised'
                            etaz = etaz(:,nonred);
                    end
                    delta1 = delta1(nonred);
                    tauh1 = cumprod(delta1);
                    Plam = bsxfun(@times,psijh1,tauh1');
                end
            end
            nofout1(i+1) = k1;
            
            if ~strcmp(method,'vanilla')
                % ----- make adaptations for cradle parameters----%
                prob = 1/exp(b0 + b1*i);                % probability of adapting
                uu = rand;
                lind = sum(abs(Gamma) < epsilon)/p;    % proportion of elements in each column less than eps in magnitude
                vec = lind >=prop;num = sum(vec);       % number of redundant columns
                
                if uu < prob
                    if  i > 20 && num == 0 && all(lind < .995)
                      %  disp(['number of columns of Gamma is: ', num2str(k2+1)]);
                        k2 = k2 + 1;
                        Gamma(:,k2) = zeros(p,1);
                        kappa(k2) = normrnd(0,1,[1,1]);
                        xi(:,k2) = normrnd(kappa(k2),1,[Ncradle,1]);
                        psijh2(:,k2) = gamrnd(df/2,2/df,[p,1]);
                        delta2 = [delta2;gamrnd(Ad2,1/bd2)];
%                         delta2(k2) = gamrnd(Ad2,1/bd2);
                        tauh2 = cumprod(delta2);
                        Pgam = bsxfun(@times,psijh2,tauh2');
                        if strcmp(method,'unsupervised')
                            xiz(:,k2) = zeros(nz,1);
                            xiz(~x,k2) = normrnd(kappa(k2),1,[sum(~x),1]);
                        end
                    elseif num > 0 && num < k2
                        nonred = setdiff(1:k2,find(vec)); % non-redundant loadings columns
                        k2 = max(k2 - num,1);
                        Gamma = Gamma(:,nonred);
                        psijh2 = psijh2(:,nonred);
                        xi = xi(:,nonred);
                        kappa = kappa(nonred);
                        delta2 = delta2(nonred);
                        tauh2 = cumprod(delta2);
                        Pgam = bsxfun(@times,psijh2,tauh2');
                        if strcmp(method,'unsupervised')
                            xiz = xiz(:,nonred);
                        end
                    end
                end
                nofout2(i+1) = k2;
            end
        end
        
        if mod(i,1000) == 0
            disp(i);
        end
        
        % -- save sampled values (without thinning) -- %
        %         %% ++++++++++++++++++++++++++++++++++++++++++++
        %         ind = mod(i,nsave);
        %         if ind == 0, ind = nsave;end
        %         output.x{ind} = x;
        %         if mod(i,nsave) == 0
        %             matobj.x(i-nsave+1:i,1) = output.x;
        %         end
        %%
        if i > burn
            ind = mod(i-burn,nsave);
            if ind == 0, ind = nsave; end
            output.Lambda{ind} = Lambda;
            output.ps{ind} = ps;
            output.eta{ind} = eta;
            switch method
                case 'vanilla'
                    output.betac{ind} = betac;
                    output.betan{ind} = betan;
                    if ~isempty(Loc)
                        Y_f = bsxfun(@times,Y_f,sqrt(VY));
                        Y_f = bsxfun(@plus,Y_f,M);
                        if strcmp(opt.transform,'curvelet')
                            output.factorCW(ind,:) = collapse_cwfeature(Y_f,Loc);
                        elseif strcmp(opt.transform,'shearlet')
                            output.factorST{ind} = Y_f;
                        end
                    end
                case {'unsupervised','supervised'}
                    if strcmp(method,'unsupervised')
                        output.x{ind} = x;
                        output.etaz{ind} = etaz;
                        output.xiz{ind} = xiz;
                        %                         output.beta{ind} = beta;
                    else
                        x = [];
                    end
                    output.kappa{ind} = kappa;
                    output.Gamma{ind} = Gamma;
                    %                     output.Vxi{i-burn} = Vxi;
                    output.rho(ind) = rho;
                    output.xi{ind} = xi;
                    if ~isempty(Loc)
                        Y_c = bsxfun(@times,Y_c,sqrt(VY));
                        Y_nc = bsxfun(@times,Y_nc,sqrt(VY));
                        Y_nc = bsxfun(@plus,Y_nc,M);
                        if strcmp(opt.transform,'curvelet')
                            output.cradleCW(ind,:) = collapse_cwfeature(Y_c,Loc([mask;~x]));
                            output.noncradleCW(ind,:) = collapse_cwfeature(Y_nc,Loc);
                        elseif strcmp(opt.transform,'shearlet')
                            output.cradleST{ind} = Y_c;
                            output.noncradleST{ind} = Y_nc;
                        end
                    end
            end
            %             output.Veta{i-burn} = Veta;
            
            if mod(i-burn,nsave) == 0
                output.psijh1 = psijh1;
                output.delta1 = delta1;
                output.tauh1 = tauh1;
                output.Plam = Plam;
                if ~strcmp(method,'vanilla')
                    output.psijh2 = psijh2;
                    output.delta2 = delta2;
                    output.tauh2 = tauh2;
                    output.Pgam = Pgam;
                end

                tic
                index = round((i-burn)/nsave);
%                 try
                    if opt.savesample
                    save([path,'_',num2str(index),'.mat'],'-struct','output');
                    save([path,'_',num2str(index),'.mat'],'M','VY','-append');
                    else
                        save([path,'_',num2str(index),'.mat'],'method','mrho');
                    end
                    try
                        switch opt.transform
                            case 'curvelet'
                        noncradleImg = invfeatureCW(output.noncradleCW,opt);
                        cradleImg = invfeatureCW(output.cradleCW,opt);
                            case 'shearlet'
                                noncradleImg = invfeatureST(output.noncradleST,opt,1,Loc);
                        cradleImg = invfeatureST(output.cradleST,opt,1,Loc([mask;~x]));
                        end
                    save([path,'_',num2str(index),'.mat'],'noncradleImg','cradleImg','-append');
                    
                    if isfield(opt,'originImg')
                        residualImg = opt.originImg - noncradleImg - cradleImg;
                     save([path,'_',num2str(index),'.mat'],'residualImg','-append');
                    end
                    end
%                 end
                %                 ind1 = i - burn - nsave + 1;
                %                 ind2 = i - burn;
                %                 matobj.Lambda(ind1:ind2,1) = output.Lambda;
                %                 matobj.ps(ind1:ind2,1) = output.ps;
                %                 matobj.eta(ind1:ind2,1) = output.eta;
                %                 switch method
                %                     case 'vanilla'
                %                         matobj.betac(ind1:ind2,1) = output.betac;
                %                         matobj.betan(ind1:ind2,1) = output.betan;
                %                         matobj.factorCW(ind1:ind2,1:p) = output.factorCW;
                %                     case {'unsupervised','supervised'}
                %                         if strcmp(method,'unsupervised')
                %                             matobj.beta(ind1:ind2,1) = output.beta;
                %                         end
                %                         matobj.kappa(ind1:ind2,1) = output.kappa;
                %                         matobj.Gamma(ind1:ind2,1) = output.Gamma;
                %                         matobj.rho(ind1:ind2,1) = output.rho;
                %                         matobj.xi(ind1:ind2,1) = output.xi;
                %                         matobj.cradleCW(ind1:ind2,1:p) = output.cradleCW;
                %                         matobj.noncradleCW(ind1:ind2,1:p) = output.noncradleCW;
                %                 end
                %                 matobj.residualCW(ind1:ind2,1:p) = output.residualCW;
                toc;
            end
        end
        
    end
    
    
    %w. evolution of factors
    nofrep(g,:) = nof1out';
    
    disp(['end replicate','',num2str(g)]);
    disp('--------------------');
    
    
    
end

%---- save parameters for restart the MCMC chain ----%
output.mean = M;
output.VY = VY;
output.psijh1 = psijh1;
output.delta1 = delta1;
output.tauh1 = tauh1;
output.Plam = Plam;
if ~strcmp(method,'vanilla')
    output.psijh2 = psijh2;
    output.delta2 = delta2;
    output.tauh2 = tauh2;
    output.Pgam = Pgam;
    if strcmp(method,'unsupervised')
        output.NCp = NCp;
    end
end

% clearvars -except 'output' 'data' 'Output' 'featureCW' 'Y' 'mask' 'test'


toc;
%%

    function c = Mat2Cell(m)
        c = mat2cell(m,size(m,1),size(m,2));
    end


    function siz = compute_original_cw_size
        angleind = opt.featureangleind;
        L = length(angleind);
        cw = fdct_wrapping(zeros(opt.imgsize),opt.curveletisreal,1,L);
        siz = cell(size(Y,2),1);
        k = 1;
        for ii = 1:L
            for jj = 1:length(angleind{ii})
                siz{k} = size(cw{ii}{angleind{ii}(jj)});
                k = k + 1;
            end
        end
    end

    function CW = collapse_cwfeature(YY,locind)
        CW = cell(1,p);
        for ii = 1:p
            cwtmp = zeros(opt.vdownsamplesize);
            cwtmp(locind) = YY(:,ii);
            %             if cradleonly
            %                 cwtmp(opt.vfeatureind.cradle) = YY(:,ii);
            %             else
            %                 cwtmp(opt.vfeatureind.free) = YY(~mask,ii);
            %                 cwtmp(opt.vfeatureind.cradle) = YY(mask,ii);
            %             end
            CW{ii} = imresize(cwtmp,cwsiz{ii},'nearest');
        end
    end

end

