cd /gtmp/rachel/SparseFactorModel/mcmc_spfact/result/

i = 1;
while exist(['job' num2str(i)],'dir')
    cd(['job' num2str(i)])
    j = 1;
    separation_result = cell(1);
    Atoms = cell(1);
    while exist(['output_' num2str(j) '.mat'],'file')
        separation_result{j} = load(['output_' num2str(j) '.mat'],'*Img');
        Atoms{j} = load(['output_' num2str(j) '.mat'],'*Atom');
        Atoms{j}.cradleAtom = cell2mat(Atoms{j}.cradleAtom);
        Atoms{j}.noncradleAtom = cell2mat(Atoms{j}.noncradleAtom);
        j = j + 1;
    end

n = length(separation_result);
noncradleImg = zeros(512,512,n);
cradleImg = zeros(512,512,n);
residualImg = zeros(512,512,n);
cradleAtom = zeros(size(Atoms{1}.cradleAtom,1),size(Atoms{1}.cradleAtom,2),n);
noncradleAtom = zeros(size(Atoms{1}.noncradleAtom,1),size(Atoms{1}.noncradleAtom,2),n);

for j = 1:n
    noncradleImg(:,:,j) = separation_result{j}.noncradleImg;
    cradleImg(:,:,j) = separation_result{j}.cradleImg;
    residualImg(:,:,j) = separation_result{j}.residualImg;
    cradleAtom(:,:,j) = Atoms{j}.cradleAtom;
    noncradleAtom(:,:,j) = Atoms{j}.noncradleAtom;    
end

noncradleImg = mean(noncradleImg,3);
cradleImg = mean(cradleImg,3);
residualImg = mean(residualImg,3);

cradleAtom = mean(cradleAtom,3);
noncradleAtom = mean(noncradleAtom,3);



save('output_final.mat','noncradleImg','cradleImg','residualImg','cradleAtom','noncradleAtom');

    cd ..
    i = i + 1;
end

% figure;imagesc(sum(noncradleImg,3)/5);axis off;axis image;colormap gray
% figure;imagesc(std(noncradleImg,[],3));axis off;axis image;colormap gray
