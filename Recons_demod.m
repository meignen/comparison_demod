function [snr_demod2_o,snr_demod3_o,snr_demod4_o] = Recons_demod(cas,SNR,ra)
 
 %input
 %cas : type of signal 
 %SNR : input SNR
 %ra  : ratio N/Nfft
 %output 
 %snr_demod2_o : output SNR for mode reconstruction using SST2-demod 
 %snr_demod3_o : output SNR for mode reconstruction using SST3-demod
 %snr_demod4_o : output SNR for mode reconstruction using SST4-demod
 
 
 if (cas == 1)
    N = 1024;
    t = (0:N-1)/N;
    a  = 2;
    s2 = a.*exp(2*pi*1i*(250*t+50*t.^3));
    s1 = a.*exp(2*pi*1i*(130*t+100*t.^2));
    s  = s1+s2;
    signal = [s1 ; s2];
    s = s(:);
    nr = 2;
    clwin = round(10/ra);
    sigma_opt = 0.05;
 elseif (cas == 2)
  %% Test signal 2
    N = 1024;
    t = (0:N-1)/N;
    a  = 2;
    s2 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
    s1 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
    s  = s1+s2;
    signal = [s1 ; s2];
    s  =s(:);
    nr = 2;
    clwin = round(30/ra);
    sigma_opt = 0.04;
 else
    N = 1024;
    t = (0:N-1)/N;
    a =  1+ 7*(1-t).^4;
    phi = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
    s = a.*exp(2*pi*1i*(phi));
    signal = s;
    s = s(:);
    clwin = round(30/ra);
    nr = 1;
    sigma_opt = 0.025;
 end

 %the window is the Gaussian window    
 
 prec  = 10^(-3);
 L     = sigma_opt*N;
 Lh    = floor(L*sqrt(-log(prec)/pi))+1;
 h     = amgauss(2*Lh+1,Lh+1,L); 
 
 Nfft = N/ra;
 
 if (Nfft < 2*Lh+1)
  error ( 'Nfft smaller than filer length');
 end
 
 n     = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s,n,SNR);
 
 gamma = 10^(-6);

 %computation of SST2, SST3 and SST4, and of the different reassignment operators
 
 [~,~,SST2,SST3,SST4,~,omega2,omega3,omega4] = sstn(sn,sigma_opt,Nfft,gamma);
 
 %demodulation using SST2, evaluation on the ridge 
 [sp2_om,integ2_om] = demod_multi_omega(sn,SST2,omega2,nr,clwin);
 
 %demodulation using SST3, evaluation on the ridge or not
 [sp3_om,integ3_om] = demod_multi_omega(sn,SST3,omega3,nr,clwin);

  %demodulation using SST4, evaluation on the ridge or not
 [sp4_om,integ4_om] = demod_multi_omega(sn,SST4,omega4,nr,clwin);

 index = Lh:N-Lh; 
 
 %signal reconstruction
 d = 0:5;
 len_d = length(d);

 snr_demod2_o = zeros(nr,len_d); %using demodulated signal with omega2
 snr_demod3_o = zeros(nr,len_d); %using demodulated signal with omega3
 snr_demod4_o = zeros(nr,len_d); %using demodulated signal with omega3
 
 for p = 1:nr
  
  %demodulation based on omega2 computed on the ridge
  [~,~,SST2_d_o,~,~] = sst2(sp2_om(p,:),sigma_opt,Nfft,gamma);
  TFR_int  = zeros(size(SST2_d_o));
  TFR_int(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:) =... 
  SST2_d_o(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:);
  [C2_o,~] = exridge(TFR_int,0,0,clwin);
  
  %demodulation based on omega3 computed on the ridge
  [~,~,SST3_d_o,~,~] = sst2(sp3_om(p,:),sigma_opt,Nfft,gamma);
  TFR_int  = zeros(size(SST3_d_o));
  TFR_int(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:) =... 
  SST3_d_o(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:);
  [C3_o,~] = exridge(TFR_int,0,0,clwin);
  
    
  %demodulation based on omega4 computed on the ridge
  [~,~,SST4_d_o,~,~] = sst2(sp4_om(p,:),sigma_opt,Nfft,gamma);
  TFR_int  = zeros(size(SST4_d_o));
  TFR_int(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:) =... 
  SST4_d_o(round(101*(Nfft/N))-clwin:round(101*(Nfft/N))+clwin,:);
  [C4_o,~] = exridge(TFR_int,0,0,clwin);
  
  for  k = 1:len_d
      
   %signal reconstruction with demodulation and second order
   %synchrosqueezing
   imf               = 1/Nfft*recmodes(SST2_d_o,C2_o,d(k));
   imf               = imf.*exp(2*1i*pi*(integ2_om(p,:)-100.*t));
   snr_demod2_o(p,k) = snr(signal(p,index),imf(index)-signal(p,index));

   %signal reconstruction with demodulation and third order
   %synchrosqueezing 
   imf               = 1/Nfft*recmodes(SST3_d_o,C3_o,d(k));
   imf               = imf.*exp(2*1i*pi*(integ3_om(p,:)-100.*t));
   snr_demod3_o(p,k) = snr(signal(p,index),imf(index)-signal(p,index));
    
   %signal reconstruction with demodulation and fourth order
   %synchrosqueezing 
   imf               = 1/Nfft*recmodes(SST4_d_o,C4_o,d(k));
   imf               = imf.*exp(2*1i*pi*(integ4_om(p,:)-100.*t));
   snr_demod4_o(p,k) = snr(signal(p,index),imf(index)-signal(p,index));

  end
 end 
end 