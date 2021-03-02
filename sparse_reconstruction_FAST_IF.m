function [ext_sig, IF_est] = sparse_reconstruction_FAST_IF(Sig, num,bw,Seg,iii,iter,win_length)
WIN=length(Sig)/2;
if (isreal(Sig))
    Sig = hilbert(Sig);
end

ext_sig=Sig;
%mean(abs(real(Seg-Sig)).^2)
for ii=1:iter
    
    Sig=ext_sig;

        %[IF_est]=FAST_IF(Sig,win_length, num, 2,100,0,0,iii); 
        [IF_est]=FAST_IF(Sig,win_length, num, 2,100,0,0); 
        
%        [IF_est]=non_tfd_IF_new_display_sparse(Sig,win_length,2,2,100,0,0);

ext_sig(iii)=0;
    %Sig=ext_sig;
   % IF_est=IF_O.';
    for i=1:num
        
        IF=IF_est(i,:);
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);%/sqrt(length(Phase));
        L=bw;
        %TF filtering for each sensor
        s1 = Sig.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(Sig));
        s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L).*hanning(2*L+1)';%.*hanning(2*L+1)';%hamming(2*L+1)';
        s2(WIN-L:WIN+L)=0;
        extr_Sig1=ifft(ifftshift(s3)).*conj(s_dechirp);
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        ext_sig(iii)=ext_sig(iii)+extr_Sig1(iii);
              
    end

%     if ii==2
%         M1=mean(abs(extr_Sig1-Seg).^2)
%      %   [tfd,orienttfd]=HTFD_new1(Sigg,2,50,84);
%         
%      %   figure; SetFigDef(16,9,'Times',20); tfsapl(real(extr_Sig1),tfd,'YLabel','Sample Number','Title','(c)', 'TFfontSize' , 20,'grayscale','on');
%         %%%%%%%%%%%%%%%%
%         t=0:127;
%         figure
%         SetFigDef(16,9,'Times',20);
%         % set(gcf,'Color','w');
%         plot(IF_est(:,1:4:end),t(1:4:end),'o',IF_O',t,'-','linewidth',2);
%         ylabel('Sample Number');
%         xlabel('Instantaneous Frequency (Hz)');
%         title('(d)');
%         axis([0 0.5 0 128]);
%         
%     end
    
    
end
for ii=1:1*iter
    Sig=ext_sig;
    ext_sig(iii)=0;
    
    for i=1:num
        
        IF=IF_est(i,:);
        
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        L=bw;
        %TF filtering for each sensor
        s1 = Sig.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(Sig));
        s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L);%.*hamming(2*L+1)';
        s2(WIN-L:WIN+L)=0;
        extr_Sig1=ifft(ifftshift(s3)).*conj(s_dechirp);
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        ext_sig(iii)=ext_sig(iii)+extr_Sig1(iii);
       
       
    end
       
  
end



