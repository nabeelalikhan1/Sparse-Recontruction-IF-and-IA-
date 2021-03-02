clc
clear all;
close all
SampFreq = 256/2;
addpath('E:\TFSA7\TFSA7');
t = 0:1/SampFreq:1-1/SampFreq;


Sig1 = 1*exp(1i*(1*pi*(25*t.^3))+1i*(1*pi*(10*t))); %300t»òÕß150t
%Sig4 = 1*exp(1i*(-1*pi*(40*t.^3))+1i*(1*pi*(115*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(1*pi*(25*t.^3))+1i*(1*pi*(40*t))); %300t»òÕß150t

Sig =1*Sig1 +1*Sig2;
%Sig=hamming(length(Sig)).'.*Sig;
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

p=[10:20   30:40   50:60   70:80  100:110];
Sig(p)=0;
[NA]=find(Sig~=0);
[ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,length(Sig)/2-1, num, 2,100,0,0,NA,8);
figure;
plot(real(SigO),'b--','linewidth',3);
hold on;plot(real(Sig),'r-.','linewidth',3);
hold on; plot(real(ext_sig),'k:','linewidth',3);
xlabel('Time(s)','FontSize',20,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',20,'FontName','Times New Roman');
legend('Original Signal','Sparse Signal','Reconstructed Signal');
axis([1  128  -2.5  2.5]);
set(gca,'FontSize',20);

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


% 
% [ext_sig] = STFT_RECONSTRUCTION(Sig,WN,wind_step,NA);
% %[ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num,11,Sig,p,5,61);
% 
% Ir=HTFD_new1(ext_sig,3,8,64);
% 
% figure;imagesc(Ir);
% xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
% set(gca,'YDir','normal');
% title('(c)','FontSize',24,'FontName','Times New Roman');
% set(gca,'FontSize',20);
% 
% 
% 
% figure;
% plot(real(SigO),'b--','linewidth',3);
% hold on;plot(real(Sig),'r-.','linewidth',3);
% hold on; plot(real(ext_sig),'k:','linewidth',3);
% xlabel('Time(s)','FontSize',20,'FontName','Times New Roman');
% ylabel('Amplitude','FontSize',20,'FontName','Times New Roman');
% legend('Original Signal','Sparse Signal','Reconstructed Signal');
% axis([1  128  -2.5  2.5]);
% set(gca,'FontSize',20);


%     0.0100    0.0276    0.0765    0.1415    0.2160    0.2866    0.3566    0.4306