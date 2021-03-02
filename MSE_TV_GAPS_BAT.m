clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;

Sig=bat_signal.';

%load seizure;
Sig=hilbert(Sig);

%IF_O(:,3)=cccc*t.^2/2+20/2;
%IF_O(:,4)=-cccc*t.^2/2+115/2;

SigO=Sig;

%IF_O(:,3)=90*t.^2/2+15;
WN=64;
wind_step=32;

%Sig=Sig.*([1:128 128:-1:1]);
num=3;
NS=100;
% HADTFD BASED
iiii=0;
%for snr=-10:2:10
 for N_S=15:5:30
         iiii=iiii+1;
    for k1=1:NS
        Sig=SigO;
        p=[];
        for i=2:7
        pp = 50*(i-1)+ randperm(50-N_S-1,1);
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
                
                             [ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num,11,Sig,p,5,121);
                
            elseif kkkkk==1 %the new algorithm
                
                [ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,121, num, 2,100,0,0,NA,7);
              %  [ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,121, num, 2,50,0,0,NA,5);
ext_sig_iccd=ext_sig;
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
    mse_ICDD(iiii)=mean(mse_FAST_IF_ICDD)
    mse_TF(iiii)=mean(mse_FAST_IF_TF_FILTER)
    mse_ST(iiii)=mean(mse_STFT)
 end
 plot(real(Sig))
 hold on; plot(real(SigO),'r:');
 hold on; plot(real(ext_sig_iccd),'k-');
 
% mse_ICDD =
% 
%     0.0044    0.0065    0.0092    0.0163
% 
% 
% mse_TF =
% 
%     0.0046    0.0086    0.0156    0.0265
% 
% 
% mse_ST =
% 
%     0.0177    0.0286    0.0367    0.0442
