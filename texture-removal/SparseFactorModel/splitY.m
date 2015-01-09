function output = splitY(output,Y,mask,struct_output)

% this function split input Y by first normalizing it and then split
% according to the posterior draws of MCMC

% Input :  output:   output of spfactcovest_mgploading after truncated ,
%                  either a struct OR a string of the path2saved .mat file
%          struct_output: (required if output is an .io.MatFile) a struct
%                         that contains mean and VY

if isstruct(output)
    M = output.mean;
    VY = output.VY;
    method = output.method;
else
    if isstr(output)
        path = output;
        output = matfile(output,'Writable',logical(1));
    end
    M = struct_output.mean;
    VY = struct_output.VY;
    method = struct_output.method;
end

% ----- renormalize the data ---- %
Y = bsxfun(@minus,Y,M);                     % center the training data
Y = bsxfun(@times,Y,1./sqrt(VY));           % scale the training data

nsample = 100;%length(output.Lambda);
Ys = struct();
if isstruct(output)
    output.Ys_noncradle = cell(nsample,1);
    output.Ys_residual = cell(nsample,1);
    if isfield(output,'Gamma')
        output.Ys_cradle = cell(nsample,1);
    end
else
    Ys_residual = cell(nsample,1);
    save(path,'Ys_residual','-append');
    clear Ys_residual
    average_residual = zeros(size(Y));
    save(path,'average_residual','-append');
    clear average_residual
    if ~isempty(output.Gamma)
        Ys_noncradle = cell(nsample,1);
        save(path,'Ys_noncradle','-append');
        clear Ys_noncradle
        Ys_cradle = cell(nsample,1);
        save(path,'Ys_cradle','-append');
        clear Ys_cradle
        average_cradle = zeros(size(Y));
        save(path,'average_cradle','-append');
        clear average_cradle
        average_noncradle = zeros(size(Y));
        save(path,'average_noncradle','-append');
        clear average_noncradle
    else
        Ys_factor = cell(nsample,1);
        save(path,'Ys_factor','-append');
        clear Ys_factor
        average_factor = zeros(size(Y));
        save(path,'average_factor','-append');
        clear average_factor
    end
end

nsets = ceil(nsample/1000);

for i = 1:nsets
    ind1 = 1000*(i-1) + 1;
    if i < nsets
        ind2 = 1000*i;
    else
        ind2 = nsample;
    end
    switch method
        case 'vanilla'
            output.Ys_factor(ind1:ind2,1) = cellfun(@(x,y)x*y',output.eta(ind1:ind2,1),output.Lambda(ind1:ind2,1),'UniformOutput',0);
            output.Ys_residual(ind1:ind2,1) = cellfun(@(x)(Y - x), output.Ys_noncradle(ind1:ind2,1),'UniformOutput',0);
        otherwise
            mask = mask(:)';
            rho = output.rho(ind1:ind2,1)*mask + ones(length(output.rho(ind1:ind2,1)),1)*(~mask);
            rho = reshape(rho',[],1);
            eta = cell2mat(output.eta(ind1:ind2,1));
            eta = bsxfun(@times,eta,rho);
            eta = mat2cell(eta,size(eta,1)/(ind2-ind1+1)*ones(ind2-ind1+1,1),size(eta,2));
            output.Ys_noncradle(ind1:ind2,1) = cellfun(@(x,y)x*y',eta,output.Lambda(ind1:ind2,1),'UniformOutput',0);
            output.Ys_cradle(ind1:ind2,1) = cellfun(@(x,y)x*y',output.xi(ind1:ind2,1),output.Gamma(ind1:ind2,1),'UniformOutput',0);
            output.Ys_residual(ind1:ind2,1) = cellfun(@(x)(Y-x),output.Ys_noncradle(ind1:ind2,1),'UniformOutput',0);
            matresidual = cell2mat(output.Ys_residual(ind1:ind2,1)');
            matresidual(mask,:) = matresidual(mask,:) - cell2mat(output.Ys_cradle(ind1:ind2,1)');
            output.Ys_residual(ind1:ind2,1) = mat2cell(matresidual,size(matresidual,1),...
                size(matresidual,2)/(ind2-ind1+1)*ones(ind2-ind1+1,1))';
    end
end

% ---- compute the average and scale back the data ----- %
if ~strcmp(method,'vanilla')
    average_noncradle = reshape(cell2mat(output.Ys_noncradle'),[],nsample);
    average_noncradle = reshape(mean(average_noncradle,2),size(Y,1),size(Y,2));
    average_noncradle = bsxfun(@times,average_noncradle,sqrt(VY));
    average_noncradle = bsxfun(@plus,average_noncradle,M);
    average_cradle = reshape(cell2mat(output.Ys_cradle'),[],nsample);
    average_cradle = reshape(mean(average_cradle,2),[],size(Y,2));
    average_cradle = bsxfun(@times,average_cradle,sqrt(VY));
    output.average_noncradle = average_noncradle;
    output.average_cradle = average_cradle;
else
    average_factor = reshape(cell2mat(output.Ys_factor'),[],nsample);
    average_factor = reshape(mean(average_factor,2),size(Y,1),size(Y,2));
    average_factor = bsxfun(@times,average_factor,sqrt(VY));
    average_factor = bsxfun(@plus,average_factor,M);
    output.average_factor = average_factor;
end
average_residual = reshape(cell2mat(output.Ys_residual'),[],nsample);
average_residual = reshape(mean(average_residual,2),size(Y,1),size(Y,2));
average_residual = bsxfun(@times,average_residual,sqrt(VY));

output.average_residual = average_residual;

