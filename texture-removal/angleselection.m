function cwnew = angleselection(cw,opt)

% Input: opt:  structure
%              .direction  'horizontal' or 'vertical' select the
%                          corresponding quads in frequency band

if nargin == 1
    opt = struct();
end

if ~isfield(opt,'direction')
    opt.direction = 'horizontal';
end

cwnew = cell(size(cw));

for i = 1:length(cwnew)
    cwnew{i} = cell(size(cw{i}));
    nangle = length(cwnew{i});
    if nangle == 1
        cwnew{i}{1} = zeros(size(cw{i}{1}));
    else
        %------ complex & real curvelet coefficients ------%
        % complex: cw{i}{j} = conj(cw{i}{j+nangle/2}
        % real: cw{i}{j} and cw{i}{j+nangle/2} are real and imaginary parts
        % of complex curvelet coefficients
        if strcmp(opt.direction,'horizontal')
%             if isfield(opt,'curveletisreal') & opt.curveletisreal
                cwnew{i}(1:nangle/4) = cw{i}(1:nangle/4);
                cwnew{i}(nangle/4+1:nangle/2) = cellfun(@(x)zeros(size(x)),cw{i}(nangle/4+1:nangle/2),'UniformOutput',0);
                cwnew{i}(nangle/2+1:nangle*3/4) = cw{i}(nangle/2+1:nangle*3/4);
                cwnew{i}(nangle*3/4+1:nangle) = cellfun(@(x)zeros(size(x)),cw{i}(nangle*3/4+1:nangle),'UniformOutput',0);
%             else
%                 cwnew{i}(1:nangle/2) = cw{i}(1:nangle/2);
%                 cwnew{i}(nangle/2+1:nangle) = cellfun(@(x)zeros(size(x)),cw{i}(nangle/2+1:nangle),'UniformOutput',0);
%             end
        else
%             if isfield(opt,'curveletisreal') & opt.curveletisreal
                cwnew{i}(1:nangle/4) = cellfun(@(x)zeros(size(x)),cw{i}(1:nangle/4),'UniformOutput',0);
                cwnew{i}(nangle/4+1:nangle/2) = cw{i}(nangle/4+1:nangle/2);
                cwnew{i}(nangle/2+1:nangle*3/4) = cellfun(@(x)zeros(size(x)),cw{i}(nangle/2+1:nangle*3/4),'UniformOutput',0);
                cwnew{i}(nangle*3/4+1:nangle) = cw{i}(nangle*3/4+1:nangle);
%             else
%                 cwnew{i}(1:nangle/2) = cellfun(@(x)zeros(size(x)),cw{i}(1:nangle/2),'UniformOutput',0);
%                 cwnew{i}(nangle/2+1:nangle) = cw{i}(nangle/2+1:nangle);
%             end
        end
    end
end

if isfield(opt,'imgsize')
    imgnew = ifdct_wrapping(cwnew,1,opt.imgsize(1),opt.imgsize(2));
    if isfield(opt,'display') && opt.display
        figure;imshow(imgnew,[]);
    end
end
