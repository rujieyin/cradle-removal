% function analyze_output

i = 2;
while (exist([opt.matname,'_',num2str(i),'.mat'],'file'))
    path = [opt.matname,'_',num2str(i),'.mat'];
    
    % --- look at average Dictionaries of chunks -- %
        Output = factorloading2dictionary(path,opt);
        nrow = max(factor(length(Output.cradleAtom)));
        figure;imagesc(cell2mat(reshape(Output.cradleAtom,nrow,[])));axis image;axis off;colormap gray;
        title(['Averaged Cradle Atom set#',num2str(i)]);
        nrow = max(factor(length(Output.noncradleAtom)));
        figure;imagesc(cell2mat(reshape(Output.noncradleAtom,nrow,[])));axis image;axis off;colormap gray;
        title(['Averaged Non Cradle Atom set#',num2str(i)]);
    
    % --- look at source separation of chunks -- %
    tmp = load(path,'*CW');
    try
        tmp.noncradleImg = invfeatureCW(tmp.noncradleCW,opt);
        rmfield(tmp,'noncradleCW');
        tmp.cradleImg = invfeatureCW(tmp.cradleCW,opt);
        rmfield(tmp,'cradleCW');
        tmp.residualImg = invfeatureCW(tmp.residualCW,opt);
        rmfield(tmp,'residualCW');
        
    end
    figure;imagesc([tmp.noncradleImg,tmp.cradleImg,tmp.residualImg]);
%     figure;imagesc(abs(cat(3,tmp.noncradleImg,tmp.cradleImg,tmp.residualImg))/2);title('superposition of 3 components');
    axis image;axis off;colormap gray;title(['noncradle -- cradle -- residual set #',num2str(i)])
    save(path,'-struct','tmp','-append');
    i = i+1;
end
