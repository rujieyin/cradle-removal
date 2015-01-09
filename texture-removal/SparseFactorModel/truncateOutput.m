function output = truncateOutput(output,n)

% input:  output: structure that contains posterior draws
%                   OR  .io.MatFile stores massive posterior draws
%         n: integer (optional) number of draws to keep

if ~isstruct(output)
    output = matfile(output,'Writable',logical(1));
end
if nargin < 2
        k1 = cellfun(@(x)size(x,2),output.Lambda);
        k2 = cellfun(@(x)size(x,2),output.Gamma);
    n1 = length(k1) - find(k1 ~= k1(end),1,'last');
    n2 = length(k2) - find(k2 ~= k2(end),1,'last');
    n = min(n1,n2);
end

N = length(output.Lambda);

% output.Lambda = output.Lambda(N-n+1:N,1);
% % output.Veta = output.Veta(end-n+1:end);
% output.ps = output.ps(N-n+1:N,1);
% output.eta = output.eta(N-n+1:N,1);
% if isfield(output,'betac')
%     output.betac = output.betac(N-n+1:N,1);
%     output.betan = output.betan(N-n+1:N,1);
% else
%     if isfield(output,'beta')
%         output.beta = output.beta(N-n+1:N,1);
%     end
%     output.kappa = output.kappa(N-n+1:N,1);
%     output.Gamma = output.Gamma(N-n+1:N,1);
% %     output.Vxi = output.Vxi(end-n+1:end);
%     output.rho = output.rho(N-n+1:N,1);
%     output.xi = output.xi(N-n+1:N,1);
% end


output.Lambda(1:N-n,1) = []; 
output.ps(1:N-n,1) = [];
output.eta(1:N-n,1) = [];
if isfield(output,'betac')
    output.betac(1:N-n,1) = [];
    output.betan(1:N-n,1) = [];
else
    if isfield(output,'beta')
        output.beta(1:N-n,1) = [];
    end
    output.kappa(1:N-n,1) = [];
    output.Gamma(1:N-n,1) = [];
    output.rho(1:N-n,1) = [];
    output.xi(1:N-n,1) = [];
end

