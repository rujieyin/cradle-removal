function [p,pthreshold] = labelpostprob(output,pi,Ynew)

% this function compute the conditional posterior probability of response
% of new data points Ynew

        M=output.mean;VY=output.VY;
        Ynew = bsxfun(@minus,Ynew,M);                     % center the training data
        Ynew = bsxfun(@times,Ynew,1./sqrt(VY));           % scale the training data
        sigma = cellfun(@(x,y)(x*x'+diag(1./y)),output.Lambda,output.ps,'UniformOutput',0);
if length(pi) > 1 % vanilla mixed non-zero mean factor model
        mu1 = cellfun(@(x,y)y*x',output.Lambda,output.betac,'UniformOutput',0);
        mu2 = cellfun(@(x,y)y*x',output.Lambda,output.betan,'UniformOutput',0);
        p = cellfun(@(mu1,mu2,Sigma,pi)multinormprob(Ynew,mu1,mu2,Sigma,pi),mu1,mu2,sigma,num2cell(pi'),'UniformOutput',0);
        p = cell2mat(p(:)');
        p = sum(p,2)/size(p,2);
        p = reshape(p,[],2);
elseif pi == 0 % vanilla model with only non-cradling input
    p = cellfun(@(x)logmultinormprob(Ynew,x),sigma,'UniformOutput',0); % average over the confidence region
    p = cell2mat(p(:)');
    p = sum(p,2)/size(p,2);
    alpha = .05; % .95 confidence interval
    pthreshold = chi2inv(1-alpha,size(output.Lambda{1},1));
end
    
        
    function prob = multinormprob(Y,mu1,mu2,Sigma,pi)
        [U,S] = svd(Sigma);
        S = sqrt(diag(S));
        prob1 = bsxfun(@minus,Y,mu1)*U*diag(1./S);
        prob1 = sum(prob1.^2,2);                          % n x 1
        prob2 = bsxfun(@minus,Y,mu2)*U*diag(1./S);
        prob2 = sum(prob2.^2,2);
        prob = [prob1,prob2];
%         prob1 = arrayfun(@(x)sum(((Y(x,:)-mu)*U*diag(1./S)).^2),1:size(Y,1));
%         prob2 = arrayfun(@(x)sum(Y(x,:)*U*diag(1./S).^2),1:size(Y,1));
        prob = [prob2-prob1+log(1-pi)-log(pi),prob1-prob2+log(pi)-log(1-pi)];
        prob = exp(-.5*prob);
        prob = 1./(1+prob);
        prob = prob(:);
    end
    
    function prob = logmultinormprob(Y,Sigma)
    [U,S] = svd(Sigma);
    S = sqrt(diag(S));
    prob = Y*U*diag(1./S);
    prob = sum(prob.^2,2);
    end
    
end
        