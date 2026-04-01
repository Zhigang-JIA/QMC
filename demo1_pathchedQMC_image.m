close all
clc
clear
addpath(genpath(cd));
%% standard image
 for iim =[1] %1:6 %[1 3 4 5][7 10 15]
    switch iim
        case 1  % Pepper
            im=imread('DATASET/image_Peppers512rgb.png');
            imname = '-Pepper512'; fname = 'Pepper512';                
    end 
   

%% double   generate A
im=imresize(im,[32,32]);
im=double(im)/255;

[ma,na,~]=size(im);
A0=[zeros(size(im(:,:,1))),im(:,:,2),im(:,:,1),im(:,:,3)];

figure;imshowQ(A0);
% figure; spy(A)
%[A0,A1,A2,A3]=A2A0123(A0);
[M,N]=size(A0);N=N/4;



%% main progra %%  quaternion rpca coded by Zhigang jia 

%Input: 
    %(1) generate low rank quaternion matrix / color image
        toy_rank=min(M,N);
        if toy_rank==min(M,N)
            L0=A0; %  original color image
        else
            L0=lowrankQ(M,4*N,toy_rank,A0);  % low rank quaternion matrix
        end

    %(2) add noise and missing poision
         denseNOIS=0.1; % generate the random positions of S0, denoted by Omega1
         if denseNOIS>0
            rng('shuffle');
            sizeNOIS=fix(denseNOIS*N*M);
            rindex=randi(M,sizeNOIS,1);
            cindex=randi(N,sizeNOIS,1);
            Omega1=zeros(M,N);
            for i=1:min(length(rindex),length(cindex))
                Omega1(rindex(i),cindex(i))=nan;
            end
            OMEGA1=[Omega1,Omega1,Omega1,Omega1];

         else
             OMEGA1=zeros(M,4*N);
         end
         POSITIONNOIS= isnan(OMEGA1);
         S0=full(sprandQ(M,N,[1.0,1.0,1.0,1.0])+0.01);% sparse quaternion matrix, the sparsity of each part of S0 is tcard!
         S0(~POSITIONNOIS)=0;
        
        denseNAN=0.1; % generate Omega of unoberserved part
            if denseNAN>0
                rng('shuffle');
                sizeNAN=fix(denseNAN*N*M);
                rindex=randi(M,sizeNAN,1);
                cindex=randi(N,sizeNAN,1);
                Omega=zeros(M,N);
                for i=1:min(length(rindex),length(cindex))
                    Omega(rindex(i),cindex(i))=nan;
                end
                OMEGA=[Omega,Omega,Omega,Omega];
            else
                OMEGA=zeros(M,4*N);
            end
            unobserved= isnan(OMEGA);
            
     % (3) generate the inpute data
        A=L0+S0+OMEGA;  %L0: original low rank quaternion matrix; S0: sparse quaterion matrix; OMEGA: missing poistions NAN
        
        fprintf(1, 'denseNOIS = %d\n', denseNOIS);
        fprintf(1, 'denseNAN = %d\n', denseNAN);

%load('50_50_demo_2.mat')
%% Solvers:
    %Method 1: QMC by Jia-Ng-Song 2018
    [fname,'-QMC']
        tol=1e-4;
        max_iter=500;
        lambda=1/sqrt((1-denseNAN)*max(M,N));
        mu=1*lambda;
        t0=cputime;
        [L1,S1,itsteps1] = rpcaQ_alm(A, lambda, mu, tol, max_iter,1);
        t1=cputime-t0;
        u_qmc=qm2im(L1,N);
        S_qmc=qm2im(S1,N);

% imshow(S_qmc)
 
        
     %% Method 2: patchedQMC  2018
     [fname,'-patchedQMC']
     input=qm2im(A,N);
     Ini=u_qmc;
     par.denseNAN=denseNAN;
     par.patsize1=16;%7 size of patch(row)
     par.patsize2=16;%  size of patch(colum)
     par.patnum=160;%55,60,65 number of each group 
     par.step1=floor((par.patsize1-1)); % step of choosing patch (row)
     par.step2=floor((par.patsize2-1)); % step of choosing patch (column)
     t0=cputime;
     u_qmc_patched= rpcaQ_alm_lansvdQ_patched(input,Ini,par);
     t3=cputime-t0;
     L3=im2qm(u_qmc_patched);

        
%% Output:
        %Data
        u_orig=qm2im(L0);
        psnrValue=[psnr(u_orig,u_qmc),psnr(u_orig,u_qmc_patched)]     
        ssimValue=[ssim(u_orig,u_qmc),ssim(u_orig,u_qmc_patched)]    
        %% Figure
        subplot(2,2,1);imshow(u_orig,[]),title('Original')
        subplot(2,2,2);imshow(u_obs,[]),title('Observed')
        subplot(2,2,3);imshow(u_qmc,[]),title('QMC')
        subplot(2,2,4);imshow(u_qmc_patched,[]),title('QMC-patched')
 end