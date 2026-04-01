%close all
clc
clear
addpath(genpath(cd));

%% Video
load ('FramefileDO01_013.mat') %horse
f=cell(32,1);
for i =1:32
  f{i}=Frame{i};
end
Frame=f;
NumFrame = size(Frame,2);
 for iframe = 1:32
     Frame{iframe} = imresize(Frame{iframe},[64,64]);
 end
[hei,wid,~] = size(Frame{1});
fname = 'horse';
psnrVedio = [];
ssimVedio = [];

denseNOIS=0; % generate the random positions of S0, denoted by Omega1
if denseNOIS>0
    rng('shuffle');
    sizeNOIS=fix(denseNOIS*hei*wid);
    rindexNOIS=randi(hei,sizeNOIS,1);
    cindexNOIS=randi(wid,sizeNOIS,1);
    
end
% (1) smapling
denseNAN=0.8; % generate Omega of unoberserved part
if denseNAN>0
    rng('shuffle');
    sizeNAN=fix(denseNAN*wid*hei);
    rindexNAN=randi(hei,sizeNAN,1);
    cindexNAN=randi(wid,sizeNAN,1);
    
end

NumFrameT = 5;
A = cell(NumFrame,1);
for iframe =1:NumFrameT
    im = Frame{iframe};
    im = im(1:hei,1:wid,:);
    %% double   generate A
    
    % im=imresize(im,[100,100]);
    im=double(im)/255;
    
    [ma,na,~]=size(im);
    A0{iframe}=[zeros(size(im(:,:,1))),im(:,:,2),im(:,:,1),im(:,:,3)];
    
    %figure;imshowQ(A0{iframe}); title('Original')
    % figure; spy(A)
    %[A0,A1,A2,A3]=A2A0123(A0);
    [M,N]=size(A0{iframe});N=N/4;
    
    
    
    %% noising and sampling
    L0=A0{iframe}; %  original color image
    %(1)noising
    %denseNOIS=0.0; % generate the random positions of S0, denoted by Omega1
    if denseNOIS>0
        %         rng('shuffle');
        %         sizeNOIS=fix(denseNOIS*N*M);
        %         rindex=randi(M,sizeNOIS,1);
        %         cindex=randi(N,sizeNOIS,1);
        Omega1=zeros(M,N);
        for i=1:min(length(rindexNOIS),length(cindexNOIS))
            Omega1(rindexNOIS(i),cindexNOIS(i))=nan;
        end
        OMEGA1=[Omega1,Omega1,Omega1,Omega1];
        
    else
        OMEGA1=zeros(M,4*N);
    end
    POSITIONNOIS= isnan(OMEGA1);
    S0=full(sprandQ(M,N,[1.0,1.0,1.0,1.0])+0.01);% sparse quaternion matrix, the sparsity of each part of S0 is tcard!
    S0(~POSITIONNOIS)=0;
    % (1) smapling
    %denseNAN=0.8; % generate Omega of unoberserved part
    if denseNAN>0
        %         rng('shuffle');
        %         sizeNAN=fix(denseNAN*N*M);
        %         rindex=randi(M,sizeNAN,1);
        %         cindex=randi(N,sizeNAN,1);
        Omega=zeros(M,N);
        for i=1:min(length(rindexNAN),length(cindexNAN))
            Omega(rindexNAN(i),cindexNAN(i))=nan;
        end
        OMEGA=[Omega,Omega,Omega,Omega];
    else
        OMEGA=zeros(M,4*N);
    end
    unobserved= isnan(OMEGA);
    
    %(3) generate the inpute data
    A{iframe}=L0+S0+OMEGA;  %L0: original low rank quaternion matrix; S0: sparse quaterion matrix; OMEGA: missing poistions NAN
    
    %  fprintf(1, 'denseNOIS = %d\n', denseNOIS);
    %  fprintf(1, 'denseNAN = %d\n', denseNAN);
    
    %   figure;imshowQ(A{iframe}); title('Observed')
    
end

%% Block-patchedQMC-windowsize  2020
hatA=zeros(M,4*N);
hatA_patched=zeros(M,4*N);
% blksize1=27;horse
% blksize2=70;barbara
% blksize2=32;horse
blksize1=16;
blksize2=16;
blknum1=M/blksize1;
blknum2=N/blksize2;
%     blktime=[];
%     blktime_patched=[];
%     blkpsnr=[];
%     blkssim=[];
%     blkpsnr_patched=[];
%     blkssim_patched=[];
blktime=cell(NumFrameT,1);
blktime_patched=cell(NumFrameT,1);
blkpsnr=cell(NumFrameT,1);
blkssim=cell(NumFrameT,1);
blkpsnr_patched=cell(NumFrameT,1);
blkssim_patched=cell(NumFrameT,1);
blku_qmc=cell(NumFrameT,1);
blkA=cell(NumFrameT,1);
hatA=cell(NumFrameT,1);
hatA_patched=cell(NumFrameT,1);
blku_obs = cell(NumFrameT,1);
for i1=1:blknum1
    for i2=1:blknum2
        '(i1,i2)='
        blki1i2=[i1,i2]
        for iframe=1:NumFrameT
             'iframe='
             iframe
            blkr=((i1-1)*blksize1+1):i1*blksize1;
            blkc=[((i2-1)*blksize2+1):i2*blksize2,N+(((i2-1)*blksize2+1):i2*blksize2),...
                2*N+(((i2-1)*blksize2+1):i2*blksize2),3*N+(((i2-1)*blksize2+1):i2*blksize2)];
            blkA{iframe}=A{iframe}(blkr,blkc);
            blku_obs{iframe}=qm2im(blkA{iframe},blksize2);
            
            %(1)initial step
            'initial step by QMC'
            tol=1e-4;
            max_iter=100;
            lambda=1/sqrt((1-denseNAN)*max(blksize1,blksize2));
            mu=1*lambda;
            t0=cputime;
            [L1{iframe},S1,itsteps1] = rpcaQ_alm_lansvdQ(blkA{iframe}, lambda, mu, tol, max_iter,0);
            blkt1=cputime-t0;
            blku_qmc{iframe}=qm2im(L1{iframe},blksize2);
            
            blkA0{iframe}=A0{iframe}(blkr,blkc);
            blku_orig{iframe}=qm2im(blkA0{iframe});
            %         blkS_qmc=qm2im(S1,blksize2);
            %
            hatA{iframe}(blkr,blkc)=L1{iframe};
            blktime{iframe}=[blktime{iframe},blkt1];
            blkpsnr{iframe}=[blkpsnr{iframe},psnr(blku_orig{iframe},blku_qmc{iframe})];
            blkssim{iframe}=[blkssim{iframe},ssim(blku_orig{iframe},blku_qmc{iframe})];
        end
        
        % (2)main step
        for iframe=1:NumFrameT
            'Main step by patched_QMC'
              'iframe='
             iframe
            % blku_obs{iframe}=qm2im(blkA{iframe},blksize2);
            %input=blku_obs{iframe};
            input=blku_obs;
            Ini=blku_qmc;
            par.denseNAN=denseNAN;
            par.tol=1.0e-4;
            par.max_iter=100;%500;
            par.stepout=0;
            par.patsize1=6;%16;%7 size of patch(row)
            par.patsize2=6;%16;%7 size of patch(colum)
            par.patnum=60;%160;%55,60,65
            par.step1=floor((par.patsize1-1)); %step of choosing patch (row)
            par.step2=floor((par.patsize2-1)); %step of choosing patch (colum)
            par.SearchWin=20;% 20;%50;
            t0=cputime;
            blku_qmc_patched{iframe}= rpcaQ_alm_lansvdQ_patched_largescale_video(input,Ini,par,NumFrameT,iframe);
            blkt3=cputime-t0;
            L3{iframe}=im2qm(blku_qmc_patched{iframe});
            %
            hatA_patched{iframe}(blkr,blkc)=L3{iframe};
            blktime_patched{iframe}=[blktime_patched{iframe},blkt3];
            blkpsnr_patched{iframe}=[blkpsnr_patched{iframe},psnr(blku_orig{iframe},blku_qmc_patched{iframe})];
            blkssim_patched{iframe}=[blkssim_patched{iframe},ssim(blku_orig{iframe},blku_qmc_patched{iframe})];
        end
    end
end

%% Output:
%Data
for iframe =1:NumFrameT
u_orig=qm2im(A0{iframe});
u_obs=qm2im(A{iframe});
u_qmc=qm2im(hatA{iframe});
u_qmc_patched=qm2im(hatA_patched{iframe});
psnrValue=[psnr(u_orig,u_qmc),psnr(u_orig,u_qmc_patched)]
psnrVideo = [psnrVideo;psnrValue];
ssimValue=[ssim(u_orig,u_qmc),ssim(u_orig,u_qmc_patched)]
ssimVideo = [ssimVideo;ssimValue];
% Figure
figure
subplot(2,3,1);imshow(u_orig,[]),title('Original')
subplot(2,3,2);imshow(u_obs,[]),title('Observed')
subplot(2,3,4);imshow(u_qmc,[]),title('QMC')
subplot(2,3,6);imshow(u_qmc_patched,[]),title('QMC-patched')
end


