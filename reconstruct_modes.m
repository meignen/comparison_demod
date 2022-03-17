function [SNR_modes,inter_moy,coeff_util] = reconstruct_modes(cas,cas1,SNR,downsamp,ra)
 %INPUT
 %cas        : type of studied signal 
 %cas1       : if 0 then classical hard-thresholding, if 1 we use the variant 
 %SNR        : SNR corresponding to the added noise
 %downsamp   : downsmpling factor in time
 %ra         : fixes the ratio between N and Nfft
 
 %OUTPUT
 %SNR_modes  : SNR associated with signal reconstruction
 %inter_moy  : average number of coefficients used at considered time instants
 %coeff_util : total number of coefficients used for mode reconstruction
 
 
 if (cas == 1)
    a  = 2;
    N = 1024;
    t = (0:N-1)/N;
    s2 = a.*exp(2*pi*1i*(250*t+50*t.^3));
    s1 = a.*exp(2*pi*1i*(130*t+100*t.^2));
    s  = s1+s2;
    signal = [transpose(s1) transpose(s2)];
    s = s(:);
    nr = 2;
    clwin = round(30/ra);
    sigma_opt = 0.05;
elseif (cas == 2)
%% Test signal 2
    a  = 2;
    N = 1024;
    t = (0:N-1)/N;
    s2 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
    s1 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
    s  = s1+s2;
    signal = [transpose(s1) transpose(s2)];
    s = s(:);
    nr = 2;
    clwin = round(30/ra);
    sigma_opt = 0.04;
 else
    N = 1024;
    t = (0:N-1)/N;
    a =  1+ 7*(1-t).^4;
    phi = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));
    s = a.*exp(2*pi*1i*(phi));
    signal = transpose(s);
    s = s(:);
    clwin = round(30/ra);
    nr = 1;
    sigma_opt = 0.025;
 end
 
 %the window is the Gaussian window    
 prec = 10^(-3);
 L =  sigma_opt*N;
 Lh = floor(L*sqrt(-log(prec)/pi))+1;
 h = amgauss(2*Lh+1,Lh+1,L); 
 
 %Lh
 n = randn(N,1)+1i*randn(N,1);
 [sn] = sigmerge(s,n,SNR);
 
 Nfft = N/ra;
 [tfr,~] = tfrstft_three_case_down(sn,Nfft,2,h,Lh,downsamp,0); 
 Abstfr = abs(tfr);
 
 %estimation of the noise level
 Y2 = real(tfr);
 gamma_estime = median(abs(Y2(:)))/0.6745;
 
 %ridge extraction
 [Cs]  = exridge_mult(tfr,nr,0,0,clwin);
 Cs    = Cs';
 modes = zeros(N,nr); 
 index = Lh:N-Lh;

 interval  = zeros(N,nr);
 SNR_modes = zeros(nr,1); 
 
 for j=1:nr
  %construction of the TF mask
  tfr_int = zeros(Nfft,N/downsamp);  
  for r = 1:N/downsamp,  
   
   val = 3*gamma_estime; %threshold for the transform depending on the noise level
  
   if (Abstfr(Cs(r,j),r) > val)
    k1 = 0;
    k2 = 0;
     
    eta1 = - 1;
    while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > val)
     if (k1 ~= Cs(r,j)-1)
      k1 = k1+1;
     else
      eta1 = k1;   
     end
    end
    if (eta1 < 0)
     eta1 = k1-1;
    end
     
    eta2 = -1;
    while (eta2 < 0) && (Abstfr(Cs(r,j)+min(k2,Nfft-Cs(r,j)),r) > val)
     if (k2 ~= Nfft-Cs(r,j))
      k2 = k2+1;   
     else
      eta2 = k2;   
     end
    end
    if (eta2 < 0)
     eta2 = k2;
    end
    
    if (cas1 == 0)
     X = max(1,Cs(r,j)-eta1):min(Nfft,Cs(r,j)+eta2);
     interval(r,j) = length(X);
     tfr_int(X,r)  = tfr(X,r);
    else
     X = max([1 Cs(r,j)-eta1 Cs(r,j)-clwin]):min([Nfft Cs(r,j)+eta2 Cs(r,j)+clwin]);  
     interval(r,j) = length(X);
     tfr_int(X,r)  = tfr(X,r);
    end
   end
  end
  modes(:,j) = itfrstft_three_case_down(tfr_int,2,N,h,0);
  SNR_modes(j) = snr(signal(index,j),modes(index,j)-signal(index,j));
 end
 
 inter_moy   = mean(interval);
 coeff_util  = sum(interval);  
  
 