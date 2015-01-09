function cwnew = ApplyMask2CW(cwmask,cw)
% this function apply cwmask to cw

if length(cwmask) ~= length(cw)
    disp('mask size wrong')
    cwnew = [];
    return;
end

cwnew = cell(size(cw));
for i = 1:length(cw)
    if length(cwmask{i}) ~= length(cw{i})
        disp(sprintf('mask size wrong : %d', i ));
        return;
    end
    cwnew{i} = cellfun(@(x,y)x.*y, cw{i},cwmask{i},'UniformOutput',0);
end