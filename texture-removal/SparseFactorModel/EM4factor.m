function [factor,pi] = EM4factor (Ynew,method, output,mask)

% this function estimate the latent factors of new data points Ynew using
% maximum posterior likelihood


        M=output.mean;VY=output.VY;
if strcmp(method,'binary')
    Ynew1 = [Ynew,ones(size(Ynew,1),1)];
    Ynew1 = rescale(Ynew1,M,VY);
    Ynew2 = [Ynew,zeros(size(Ynew,1),1)];
    Ynew2 = rescale(Ynew2,M,VY);
else
    Ynew = rescale(Ynew,M,VY);
end

ndraw = length(output.Lambda);
a = 1; b = 1;
pi = betarnd(a+sum(mask),b+sum(~mask),[1,ndraw]);



if strcmp(method,'nonzero-mean')

    factor = cellfun(@(Veta,Lambda,ps,betac,betan)getmean(Veta,Lambda,ps,method,Ynew,betac,betan),...
        output.Veta,output.Lambda,output.ps,output.betac,output.betan,'UniformOutput',0);
elseif strcmp(method,'binary')
        factor = cellfun(@(Veta,Lambda,ps)getmean(Veta,Lambda,ps,method,[Ynew1;Ynew2],[]),...
        output.Veta,output.Lambda,output.ps,'UniformOutput',0);
end
    factor = cell2mat(factor(:)');                    % k*2n x ndraw
    factor = factor(1:size(factor,1)/2,:)*pi(:)+ factor(size(factor,1)/2+1:end,:)*(1-pi(:));
    factor = sum(factor,2)/size(factor,2);
    factor = reshape(factor,size(Ynew,1),[]);

    

    function Ynew = rescale(Y,M,VY)
        Ynew = bsxfun(@minus,Y,M);                     % center the training data
        Ynew = bsxfun(@times,Ynew,1./sqrt(VY));           % scale the training data
    end

    function Meta = getmean(Veta, Lambda, ps,method,Y,betac,betan)
        Lmsg = bsxfun(@times,Lambda,ps);
        if strcmp(method,'binary')
            Meta = Y*Lmsg*Veta;                        % 2n x k
        elseif strcmp(method,'nonzero-mean')
            Meta = Y*Lmsg*Veta;
            Meta = [Meta+repmat(betac*Veta,size(Meta,1),1); Meta + repmat(betan*Veta,size(Meta,1),1)];
        end
        Meta = Meta';                                  % k x 2n
        Meta = Meta(:);
    end

end