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
N_S=20;%:5:30
iiii=iiii+1;
Sig=SigO;
p=[];
for i=2:6
    pp = 50*(i-1)+ randperm(50-N_S-1,1);
    p1=pp:1:pp+N_S;
    p=[ p p1];
end
Sig(p)=0;
[NA]=find(Sig~=0);

% ORIGINAL
delta=5;
alpha = 5;

[ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num,11,Sig,p,5,121);




plot(real(Sig),'linewidth',1)
hold on; plot(real(SigO),'r:','linewidth',1);
hold on; plot(real(ext_sig),'k-','linewidth',1);
legend('Signal with gapped missing samples','Original Signal', 'Reconstructed Signal') 
xlabel('Samples')
ylabel('Amplitude')



Is=HTFD_new1(Sig,3,8,64);
Ir=HTFD_new1(ext_sig,3,8,64);
Io=HTFD_new1(SigO,3,8,64);
figure;imagesc(Io);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(Is);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('b','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(Ir);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(c)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);


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
