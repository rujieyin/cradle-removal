% notice that for smooth cradling and non-smooth cradling, different
% methods should be used!
tic
% step 1: estimate roughly the location of horizontal and vertical cradles
vn = 1; % number of vertical cradles
vrange = [1 , size(img,1)]; % the range of vertical cradles, get rid of zero background
hn = 1; % number of horizontal cradles
hrange = [ 1 , size(img,2)]; % the range of horizontal cradles, get rid of zero background
opt = struct();
opt.L = 20;
opt.s = 10;
opt.display = 0;
[verest,horest] = cradledetect(img,vn,vrange,hn,hrange,opt);
% step 2: estimate the accuracy of the estimation of verest & horest
opt.hs = 20;
opt.vs = 50;%20;%50;
opt.smoothbd = 0; % indicate if the boundary is smooth (noisy)
% opt.deltaHIntense = 25 ; % the intensity difference induced by horizontal cradle
% opt.deltaVIntense = 30 ; % the intensity difference induced by vertical cradle

% step 3: remove all horizontal cradles
% [imgnew,~,hinfo,hangle] = RmFullHorizontalCradles(img,horest,verest,opt);
% opt.HdI = ;
[imgnew,hinfo] = RmHorizontalCradle(img,horest,verest,opt);
% step 4: remove all vertical cradles
Iflag = [ 0, 0]; % indicate whether the very left or right cradle is adjacent to zero background
opt.edgecut = 1; % whether use horizontal cut to estimate the location of vertical boundary
[I ,vinfo] = RmRegularVerticalCradleSegmentation(imgnew,hinfo,verest,opt);
% % ============optional set edge location using vinfo =====================%
% % ****** using additive model ********%
% opt.edgeloc_v = 
% opt.VdI = 
% [I, vinfo] = RmRegularVerticalCradleSegmentation(imgnew,hinfo,verest,opt);
% % ******* using multiplicative attenuation model **********%
% i1 = 1; j1 = 3; % cradle to be reprocessed
% i2 = 2; j2 = 3; % cradle whose profile to be used
% [~,~,I] = PostRmCradleProfile(imgnew,vinfo,i1,j1,i2,j2,I);
% % ========================================================================%
vinfonew = connect_vinfo(vinfo);
% [I, info] = RmFullVerticalCradleSegmentation(imgnew, Iflag, verest,horest,hinfo,opt);
[Inew,crossinfo] = Rm_cross_section(I,hinfo,vinfonew);

% %==============================================================================================%
%     
% opt.Vangle = []; % set the angle of each whole vertical cradling according to vinfo{i,j}.angle
% opt.edgeloc = []; % set the edge location estimation of each whole vertical cradling  
% [I, info] = RmFullVerticalCradleSegmentation(imgnew, Iflag, verest,horest,hinfo,opt);
% % step 5: redo some vertical segmentation if necessary
% Ind = []; % index of the vertical segmentations need to be removed based on neighborhood info
% [imgnewnew,info] = RmSpecialVerticalCradleSegmentation(I,Ind,info,opt,imgnew);
% % step 6 : remove the cross section part ? or remove the texture first ?
toc