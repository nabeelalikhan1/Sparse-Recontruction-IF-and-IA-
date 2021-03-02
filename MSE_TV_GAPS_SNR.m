clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;


Sig1 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(2*pi*(0*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(-1*pi*(30*t.^3))+1i*(1*pi*(100*t))); %300t»òÕß150t
%Sig4 = 1*exp(1i*(-1*pi*(40*t.^3))+1i*(1*pi*(115*t))); %300t»òÕß150t

Sig3 = exp(1i*(1*pi*(20*t +30*t.^3)));
Sig =1*Sig1 +0*Sig3+1*Sig2;
Sig=hamming(length(Sig)).'.*Sig;
SigO =Sig;
cccc=30*3;
IF_O(:,1)=cccc*t.^2/2;
IF_O(:,2)=-cccc*t.^2/2+100/2;
%IF_O(:,3)=cccc*t.^2/2+20/2;
%IF_O(:,4)=-cccc*t.^2/2+115/2;



%IF_O(:,3)=90*t.^2/2+15;
WN=64;
wind_step=32;

%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=100;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
iiii=0;
jjjj=0;
for snr=0:10:30
 jjjj=jjjj+1;
 iiii=0;
    for N_S=3:8 
         iiii=iiii+1;

    for k1=1:NS
        Sig=SigO;
        Sig=awgn(Sig,snr,'measured');
        p=[];
        for i=1:8
        pp = 16*(i-1)+ randperm(16-N_S-1,1);
        p1=pp:1:pp+N_S;
        p=[ p p1];
        end
        Sig(p)=0;
        [NA]=find(Sig~=0);
        for kkkkk=0:2
            
            % ORIGINAL
            delta=5;
            alpha = 5;
            if kkkkk==0   %ADTFD+VITERBI
                             [ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num,11,Sig,p,5,61);
               
            elseif kkkkk==1 %the new algorithm
                [ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,length(Sig)/2-1, num, 2,100,0,0,NA,5-2);
              %  [~,findex] = FAST_IF_ICCD_Sparse(ext_sig,length(Sig)/2-1, num, 2,100,0,0,NA,5);
             %   [ext_sig,~,~] = ICCD_sparse(Sig,1,findex,1,5,0,NA);
               % ext_sig=sum(ext_sig);
               % ext_sig(NA)=Sig(NA);
            else
               [ext_sig] = STFT_RECONSTRUCTION(Sig,WN,wind_step,NA);
            end
            
            if kkkkk==0
                                mse_FAST_IF_TF_FILTER(k1)=mean(abs(ext_sig-SigO));

            elseif kkkkk==1
                                mse_FAST_IF_ICDD(k1)=mean(abs(ext_sig-SigO));

            else
                mse_STFT(k1)=mean(abs(ext_sig-SigO));
            end
            
            
            
        end
        
        
    end
    mse_ICDD(jjjj,iiii)=mean(mse_FAST_IF_ICDD);
    mse_TF(jjjj,iiii)=mean(mse_FAST_IF_TF_FILTER);
    mse_ST(jjjj,iiii)=mean(mse_STFT);
    end
 mse_ICDD(jjjj,:)
  mse_ST(jjjj,:)

 mse_TF(jjjj,:)
end

figure;
plot(3:8,mse_ICDD(1,:),'k','linewidth',3);
hold on;
plot(3:8,mse_ST(1,:),'r','linewidth',3);
hold on;
plot(3:8,mse_TF(1,:),'b','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('a')
legend('The Proposed Method','Modified OMP','TF filtering');

figure;

plot(3:8,mse_ICDD(2,:),'k','linewidth',3);
hold on;
plot(3:8,mse_ST(2,:),'r','linewidth',3);
hold on;
plot(3:8,mse_TF(2,:),'b','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('b')
legend('The Proposed Method','Modified OMP','TF filtering');

figure;

plot(3:8,mse_ICDD(3,:),'k','linewidth',3);
hold on;
plot(3:8,mse_ST(3,:),'r','linewidth',3);
hold on;
plot(3:8,mse_TF(3,:),'b','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('c')
legend('The Proposed Method','Modified OMP','TF filtering');


figure;

plot(3:8,mse_ICDD(4,:),'k','linewidth',3);
hold on;
plot(3:8,mse_ST(4,:),'r','linewidth',3);

hold on;
plot(3:8,mse_TF(4,:),'b','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('d')
legend('The Proposed Method','Modified OMP','TF filtering');

