function [Y, Ymask, opt, Yind] = spfactorInputData(path,i,j,jobsiz)

% this function compute the input data of sparse factore model

% Inpute:     precompute:  path to the precompute.mat file

precompute = [path 'precompute.mat'];
load(precompute,'Nnoncradle','Ncradle','sample_noncradle_region','direction');
% Nnoncradle, Ncradle, sample_noncradle_region, FMstatus_vertical
switch direction
    case 'vertical'
        load(precompute,'FMstatus_vertical');
        status = FMstatus_vertical;
        orth_dir = 'horizontal';
    case 'horizontal'
        load(precompute,'FMstatus_horizontal');
        status = FMstatus_horizontal;
        orth_dir = 'vertical';
end
siz = size(status);

if strcmp(status{i,j},'unnecessary')
    Y = [];
    Ymask = [];
    opt = [];
    Yind = [];
    inferind = [];
    return;
end

r = 2; % radius of neighborhood
num_NC_region = 5; % total number of non-cradle regions to sample from
num_NCsample = 500; % number of non-cradle samples draw per region
num_Csample = 1500; % number of cradle samples draw per region
% --- get non-cradle samples -- %
NCind = get2dneighbor(i,j,siz(1),siz(2),r);
NCind = NCind(sample_noncradle_region(NCind));
disp([num2str(length(NCind)) 'non cradle regions found.'])
if length(NCind) > num_NC_region
    NCind = randsample(NCind,num_NC_region);
end
if length(NCind) < num_NC_region
    tmp = setdiff(find(sample_noncradle_region),NCind);
    tmp = randsample(tmp,num_NC_region-length(NCind));
    NCind = [NCind(:);tmp(:)];
end
[row,col] = ind2sub(jobsiz,NCind);
Ync = cell(num_NC_region,1);

for k = 1:num_NC_region
    load([path 'job' num2str(col(k)) '/output.mat'],'featureST');
    Ync{k} = featureST{row(k)}.free(randsample(Nnoncradle(NCind(k)),num_NCsample),:);
    if Nnoncradle(NCind(k)) ~= size(featureST{row(k)}.free,1)
        disp('Noncradle wrong!');
    end
end
Ync = cell2mat(Ync);
NCind = [row(:) col(:)];

% --- get cradle samples -- %
Cind = get1dneighbor(i,j,siz(1),siz(2),direction);

[row,col] = ind2sub(jobsiz,Cind);
Yc = cell(3,1);
for k = 1:3
    load([path 'job' num2str(col(k)) '/output.mat'],'featureST');
    num_sample = min(num_Csample,Ncradle(Cind(k)));
    Yc{k} = featureST{row(k)}.cradle(randsample(Ncradle(Cind(k)),num_sample),:);
    if Ncradle(Cind(k)) ~= size(featureST{row(k)}.cradle,1)
        disp('Cradle wrong!');
    end
    
end
Yc = cell2mat(Yc);
Cind = [row(:) col(:)];
Y = [Yc;Ync];
Ymask = [ones(size(Yc,1),1);zeros(num_NCsample*num_NC_region,1)];
Yind = struct();
Yind.NCind = NCind;
Yind.Cind = Cind;

% ---- define opt ---- %
opt = struct();
opt.direction = direction;
opt.method = 'supervised';
opt.transform = 'shearlet';

% inferind = get1dneighbor(i,j,siz(1),siz(2),orth_dir);
% inferind = inferind(~strcmp(status,'apply') && ~strcmp(status,'unnecessary'));
% inferind = [inferind(:);sub2ind(siz,i,j)];


    function neighind = get1dneighbor(i,j,m,n,dir)
        
        switch dir
            case 'vertical'
                if i > 1 && i < m
                    neighind = [i-1, i, i+1;j,j,j];
                else
                    if i == 1
                        neighind = [1,2,3;j,j,j];
                    else
                        neighind = [(m-2):m;j,j,j];
                    end
                end
            case 'horizontal'
                if j > 1 && j < n
                    neighind = [i,i,i;j-1,j,j+1];
                else
                    if j == 1
                        neighind = [i,i,i;1,2,3];
                    else
                        neighind = [i,i,i;n-2:n];
                    end
                end
        end
        neighind = sub2ind([m,n],neighind(1,:),neighind(2,:));
        
    end

    function neighind = get2dneighbor(i,j,m,n,r)
        mask = zeros(m,n);
        a1 = max(i-r,1);
        a2 = min(i+r,m);
        b1 = max(j-r,1);
        b2 = min(j+r,n);
        mask(a1:a2,b1:b2) = 1;
        neighind = find(mask);
    end

end