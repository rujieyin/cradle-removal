% clear vars
% img = subpatch{4,7};
% img = imread('~/Documents/simulations/cradle removal/cradleremoved.tif');
img = double(img(:,:,1));
addpath(genpath('/home/grad/rachel/matlab/MCALabWithUtilities'))
rmpath(genpath('/home/grad/rachel/matlab/MCALabWithUtilities/CurveletToolbox'))

% Dictionary stuff (here Curvelets + UDWT).
% dict1='CURVWRAP';pars11=3;pars12=0;pars13=0;
dict2 = 'CURVWRAP';pars21=3;pars22=0;pars23=0;
addpath(genpath('/home/grad/rachel/matlab/CurveLab-2.1.3/fdct_wrapping_matlab'))
% dict2='LDCT2';pars21=64;pars22=16/256;pars23=0; % Remove Low-frequencies < 16/256 from textured part.
dict1='CDWT2';pars11=6;pars12=0;pars13=0;
addpath(genpath('/home/grad/rachel/matlab/dtcwt'))
% dict2='LDCT2iv';pars21='Sine';pars22=32;pars23=128/512; % Remove Low-frequencies 128/512 from textured part.
dicts=MakeList(dict1,dict2);
pars1=MakeList(pars11,pars21);
pars2=MakeList(pars12,pars22);
pars3=MakeList(pars13,pars23);

% Call the MCA.
itermax 	= 100;
tvregparam 	= .2;
tvcomponent	= 0;
expdecrease	= 1;
lambdastop	= 1;
sigma		= 1E-6;
display		= 1;
imgsize = min(size(img),512);
img = img(1:imgsize(1),1:imgsize(2));
[parts,options]=MCA2_Bcr(img,dicts,pars1,pars2,pars3,itermax,tvregparam,tvcomponent,expdecrease,lambdastop,[],sigma,display);
options.inputdata = 'Input image: cradleremoved';
options
% [ST,I] = dbstack;
% name=eval(['which(''' ST(1).name ''')']);
% eval(sprintf('save %s options -V6',[name(1:end-2) 'metadata']));

% Display results.
% % figure;
% % set(gcf,'Name','MCA cradleremoved','NumberTitle','off');
% % subplot(221);
% % imagesc(img);axis image;rmaxis;
% % title('Original cradleremoved');
% % 
% % subplot(222);
% % imagesc(squeeze(sum(parts,3)));axis image;rmaxis;
% % title(sprintf('MCA cradleremoved PSNR=%.3g dB',psnr(img,squeeze(sum(parts,3)))));
% % 
% % subplot(223);
% % imagesc(squeeze(parts(:,:,1)));axis image;rmaxis;
% % title('cradleremoved Cartoon');
% % 
% % subplot(224);
% % imagesc(squeeze(parts(:,:,2)));axis image;rmaxis;
% % title('cradleremoved Texture');
% % 
% % colormap('gray');
