function [Sig_r]=STFT_RECONSTRUCTION(Sig,WN,wind_step,NA)
%WN=64;
%wind_step=32;
%Sig;
NA=NA+WN;
N=length(Sig);
x_p=[zeros(1,WN) Sig zeros(1,WN)];
it=0;
STFT=[];
STFT_R=[];
Sig_r=zeros(1,length(x_p));

for tt=0:WN/wind_step:length(x_p)-WN-1    % Signal length (1024+256) with a step of 32
    it=it+1;
    K=24;                % Sparsity
    x_w=x_p(tt+1:tt+WN).'.*hanning(WN,'periodic');
    X=(fft(x_w));                    % FT of the original (nonsparse) signal
    
    STFT=[STFT, X];              % STFT of the original (nonsparse) signal
      NA_v=NA(and(NA<tt+WN,NA>tt))-tt; % Random positions of the available samples
  
    % Generate the measurement matrix:
    IDFT_mtx=conj(dftmtx(WN))/WN; 
    A=IDFT_mtx(NA_v,:);             % Partial IDFT matrix
    y=A*X;                          % Available samples
  
    % Reconstruction with iterative OMP:
    KB=[];  y0=y; Kstep=1;
    for iter=1:round(K/Kstep)
        
        if iter>1; y=y0-x_R; end                 % Remove the reconstructed component from the array
        
        X0=(A*WN)'*y;                            % Initial estimate
        [Xv,Kv]=sort(abs(X0));                   % Sort the components
        KB=union(KB,Kv(end-Kstep+1:end));        % Take the largest component
        A_K=A(:,KB);                             % Measurement matrix for the largest component
        X_R=zeros(WN,1);
        X_R(KB)=pinv(A_K)*y0;                    % FT reconstruction of the components
        x_R=A*X_R;                               % Inverse of X_R
    end
    
    STFT_R=[STFT_R, X_R]; % Reconstructed STFT
    Sig_r(tt+1:tt+WN)=Sig_r(tt+1:tt+WN)+ifft(X_R).'/(wind_step/2);
    
    % Statistical error per one window:        
    E_stat(it)=10*log10(sum(abs(X_R(KB)-X(KB)).^2));
    
end
Sig_r=Sig_r(WN+1:N+WN);
end
