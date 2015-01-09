function output = factorloading2dictionary(output,opt)

% this function compute the mean of factor loading matrix across posterior
% draws as dictionaries
% Input:  output :  a struct OR path2saved .mat file where posterior draws stored

if isstruct(output)
    p = size(output.Lambda{end},1);
    k1 = size(output.Lambda{end},2);
    Lambda = output.Lambda;
    if isfield(output,'Gamma')
        Gamma = output.Gamma;
    end
else
    if isstr(output)
        path = output;
        output = struct();
        output.path = path;
    end
    load(path,'Lambda');
    try
    load(path,'Gamma');
    end
    tmp = Lambda(1);
    p = size(tmp{1},1);
    k1 = size(tmp{1},2);
    Output = struct();
end
Nsample = length(Lambda);

% ----- plot mean and std of each column ------ %
for i = 1:k1
    tmp = cellfun(@(x)x(:,i),Lambda,'UniformOutput',0);
    tmp = cell2mat(tmp');
    if isfield(opt,'display') && opt.display
        figure;errorbar(1:p,mean(tmp'),std(tmp'));title(['Lambda Dictionary #',num2str(i)]);
    end
end


% ------ compute DicLambda, DicGamma ---------- %

tmp = reshape(cell2mat(Lambda'),[p,k1,Nsample]);
tmp = sum(tmp,3)/Nsample;
if isstruct(output)
    output.DicLambda = tmp;
    output.noncradleAtom = arrayfun(@(x)invCWfeature_vector(opt,tmp(:,x)),1:k1,'UniformOutput',0);
    % figure;imshow(cell2mat(output.noncradleAtom),[]);
else
    Output.DicLambda = tmp;
    Output.noncradleAtom = arrayfun(@(x)invCWfeature_vector(opt,tmp(:,x)),1:k1,'UniformOutput',0);
    %     figure;imshow(cell2mat(Output.noncradleAtom),[]);
end

if length(who('Gamma')) > 0
        k2 = size(Gamma{end},2);
    for i = 1:k2
        tmp = cellfun(@(x)x(:,i),Gamma,'UniformOutput',0);
        tmp = cell2mat(tmp');
        %         figure;errorbar(1:p,mean(tmp'),std(tmp'));title(['Gamma Dictionary #',num2str(i)]);
    end
    tmp = reshape(cell2mat(Gamma'),[p,k2,Nsample]);
    tmp = sum(tmp,3)/Nsample;
    if isstruct(output)
        output.DicGamma = tmp;
        output.cradleAtom = arrayfun(@(x)invCWfeature_vector(opt,tmp(:,x)),1:k2,'UniformOutput',0);
        if isfield(opt,'display') && opt.display
            figure;imshow(cell2mat(output.cradleAtom),[])
        end
    else
        Output.DicGamma = tmp;
        Output.cradleAtom = arrayfun(@(x)invCWfeature_vector(opt,tmp(:,x)),1:k2,'UniformOutput',0);
        Output.source = output;
        output = Output;
    end
end

if isfield(output,'path')
    save(path,'-struct','output','-append');
end


