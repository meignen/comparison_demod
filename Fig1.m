 sigma_opt = 0.05;
 close all
 %the window is the Gaussian window    
 N = 1024;
 prec = 10^(-3);
 L =  sigma_opt*N;
 Lh = floor(L*sqrt(-log(prec)/pi))+1;
 h = amgauss(2*Lh+1,Lh+1,L); 
 
 Nfft = N;
 % the first type of signal
 t = (0:N-1)/N;
 a  = 2;
 s2 = a.*exp(2*pi*1i*(250*t+50*t.^3));
 s1 = a.*exp(2*pi*1i*(130*t+100*t.^2));
 s  = s1+s2;
 s = s(:); 
 
 [tfr,~] = tfrstft_three_case_down(s,Nfft,1,h,Lh,1,0); 
 Abstfr = abs(tfr);
 imagesc(t,(0:Nfft/2-1)*N/Nfft,Abstfr(1:Nfft/2,:));
 set(gca,'ydir','normal');
 xlabel('time');
 ylabel('frequency');
 set(gca,'fontsize',24);
 
 figure () 
 %the second type of signal
 s2 = a.*exp(2*pi*1i*(330*t+16*cos(3*pi*t)));
 s1 = a.*exp(2*pi*1i*(190*t+9*cos(3*pi*t)));
 s  = s1+s2;
 s  = s(:);
 sigma_opt = 0.04;
 
 %the window is the Gaussian window    
 N = 1024;
 prec = 10^(-3);
 L =  sigma_opt*N;
 Lh = floor(L*sqrt(-log(prec)/pi))+1;
 h = amgauss(2*Lh+1,Lh+1,L); 

 
 [tfr,~] = tfrstft_three_case_down(s,Nfft,1,h,Lh,1,0); 
 Abstfr = abs(tfr);
 imagesc(t,(0:Nfft/2-1)*N/Nfft,Abstfr(1:Nfft/2,:));
 set(gca,'ydir','normal');
 xlabel('time');
 ylabel('frequency');
 set(gca,'fontsize',24);
 
 %the third type of signal
 figure() 
 %a1 = exp(2*(1-t).^3);%exp(2*(1-t).^3 + 1.5*t.^4);
 a2 =  1+ 7*(1-t).^4;%1+ 5*t.^3 + 7*(1-t).^6;

 %phi1 = 50*t+30*t.^3-20*(1-t).^4;
 phi2 = 340*t-2.*exp(-2*(t-0.2)).*sin(14*pi.*(t-0.2));

 %s1 = a1.*exp(2*pi*1i*(phi1));
 s2 = a2.*exp(2*pi*1i*(phi2));
 s  = s2;
 s  = s(:);
 
 sigma_opt = 0.025;
 
 %the window is the Gaussian window    
 N = 1024;
 prec = 10^(-3);
 L =  sigma_opt*N;
 Lh = floor(L*sqrt(-log(prec)/pi))+1;
 h = amgauss(2*Lh+1,Lh+1,L); 
 
 
 [tfr,~] = tfrstft_three_case_down(s,Nfft,1,h,Lh,1,0); 
 Abstfr = abs(tfr);
 imagesc(t,(0:Nfft/2-1)*N/Nfft,Abstfr(1:Nfft/2,:));
 set(gca,'ydir','normal');
 xlabel('time');
 ylabel('frequency');
 set(gca,'fontsize',24);
 s = s1+s2;
 s = s(:);

 