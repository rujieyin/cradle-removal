    function siz = compute_original_cw_size(opt)
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
