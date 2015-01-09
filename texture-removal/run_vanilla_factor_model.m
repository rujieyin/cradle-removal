
% blocktexturefeature

Y = [cwfeature.vcradle';cwfeature.nonvcradle'];
cwmask = [ones(length(cwfeature.vcradle),1);zeros(length(cwfeature.nonvcradle),1)];
% case I : non-cradle training data alone
 output2 = spfactcovest_mgploadings(Y(cwmask==0,:),cwmask(cwmask==0),[]);
 output2 = truncateOutput(output2,6000); % remove burn-in
 % case II : supervised training
 output = spfactcovest_mgploadings(Y,cwmask,struct('method','vanilla'));
 output = truncateOutput(output,7000);
 
 pcw  = labelpostprob(output2,0,Y(1:7089,:));
 pcw2 = labelpostprob(output,pi,Y(1:7089,:));
 
 % generate random vector for case II
  pi = betarnd(1+sum(cwmask),1+sum(~cwmask),[1,7000]);
 
%compute the mean of the first column of loading matrix
tmp = cellfun(@(x)x(:,1),output2.Lambda,'UniformOutput',0);
tmp = cell2mat(tmp');
lambdamean = mean(tmp');
figure;imshow(invfeature(opt,lambdamean),[]);title('non-cradling dictionary')
% renormalize the input data
figure;imshow(invfeature(opt,(Y(1,:)-output2.mean)./output2.VY),[]);title('cradling')
% display the whole dictionaries ( 9 columns )
tmp = reshape(cell2mat(output2.Lambda'),[15,size(output2.Lambda{end},2),length(output2.Lambda)]);
tmp = sum(tmp,3)/length(output2.Lambda);
tmp = arrayfun(@(x)invfeature(opt,tmp(:,x)),1:size(output2.Lambda{end},2),'UniformOutput',0);
figure;imshow(cell2mat(tmp),[])
figure;imshow(cell2mat(reshape(tmp,[3,3])),[])


% case I : compute confidence interval
test = zeros(opt.vdownsamplesize);
flag = 1- chi2cdf(pcw,15);
test(opt.vfeatureind.cradle) = flag;
figure;imshow(test,[])
set(h,'AlphaData',mask.*flag.*(flag>0)/max(flag(:))/2);
flag = test;
flag = imresize(flag,[512,512]);
set(h,'AlphaData',mask.*flag.*(flag>0)/max(flag(:))/2);
figure;image(repmat((real(imgnew+20))/40,[1,1,3]));axis off;hold on
h = image(cat(3,mask,zeros(512),1-mask));
set(h,'AlphaData',mask.*flag.*(flag>0)/max(flag(:))/2);
% case II : compute conditional probabitlity
flag = pcw2(:,2)-pcw2(:,1);
flag = flag.*(flag>0);
test = zeros(opt.vdownsamplesize);
test(opt.vfeatureind.cradle) = flag;
figure;imshow(test,[])
flag = test;
flag = imresize(flag,[512,512]);
figure;image(repmat((real(imgnew+20))/40,[1,1,3]));axis off;hold on
h = image(cat(3,mask,zeros(512),1-mask));
set(h,'AlphaData',mask.*flag.*(flag>0)/max(flag(:))/2);

% plot the mean and std of the posterior draws of the first column of
% loading matrix
tmp = cellfun(@(x)x(:,3),output2.Lambda,'UniformOutput',0);
tmp = cell2mat(tmp');
size(std(tmp'))
figure;errorbar(1:15,mean(tmp'),std(tmp'))
