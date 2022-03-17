function Fig_compare_downsample_demod(cas,ra)
 
 %cas: corresponds to the type of signal we are studying
 %ra : proportion between Nfft and N
 
 nbreal = 5;

 d =0:5;
 len_d = length(d);
 
 if (cas <= 2)
  nr = 2;
 else
  nr = 1;
 end
 
 snr_demod2_o     = zeros(nbreal,nr,len_d);
 snr_demod3_o     = zeros(nbreal,nr,len_d);
 snr_demod4_o     = zeros(nbreal,nr,len_d);

 snr_demod2_o_10  = zeros(nbreal,nr,len_d);
 snr_demod3_o_10  = zeros(nbreal,nr,len_d);
 snr_demod4_o_10  = zeros(nbreal,nr,len_d);

 snr_demod2_o_30  = zeros(nbreal,nr,len_d);
 snr_demod3_o_30  = zeros(nbreal,nr,len_d);
 snr_demod4_o_30  = zeros(nbreal,nr,len_d);

 % SNR for HT reconstruction
 SNR_modes1_0 = zeros(nbreal,nr);
 SNR_modes2_0 = zeros(nbreal,nr);
 SNR_modes1_10 = zeros(nbreal,nr);
 SNR_modes2_10 = zeros(nbreal,nr);
 SNR_modes1_30 = zeros(nbreal,nr);
 SNR_modes2_30 = zeros(nbreal,nr);

 if (cas == 3)
  SNR_modes3_0  = zeros(nbreal,nr);
  SNR_modes3_10 = zeros(nbreal,nr);
  SNR_modes3_30 = zeros(nbreal,nr);
  
  SNR_modes4_0  = zeros(nbreal,nr);
  SNR_modes4_10 = zeros(nbreal,nr);
  SNR_modes4_30 = zeros(nbreal,nr);

 end
 
 % SNR for HT reconstruction, with alternative technique
 SNR_modes1_1_0 = zeros(nbreal,nr);
 SNR_modes2_1_0 = zeros(nbreal,nr);
 SNR_modes1_1_10 = zeros(nbreal,nr);
 SNR_modes2_1_10 = zeros(nbreal,nr);
 SNR_modes1_1_30 = zeros(nbreal,nr);
 SNR_modes2_1_30 = zeros(nbreal,nr);

 
 %the number of coefficients used for the reconstruction of each mode
 %for reconstrcution based on STFT (first technique)
 coeff_util1_0  = zeros(nbreal,nr);
 coeff_util2_0  = zeros(nbreal,nr);
 coeff_util1_10  = zeros(nbreal,nr);
 coeff_util2_10  = zeros(nbreal,nr);
 coeff_util1_30  = zeros(nbreal,nr);
 coeff_util2_30  = zeros(nbreal,nr);

 %the number of coefficients used for the reconstruction of each mode
 %for reconstrcution based on STFT (second technique)
 coeff_util1_1_0  = zeros(nbreal,nr);
 coeff_util2_1_0  = zeros(nbreal,nr);
 coeff_util1_1_10  = zeros(nbreal,nr);
 coeff_util2_1_10  = zeros(nbreal,nr);
 coeff_util1_1_30  = zeros(nbreal,nr);
 coeff_util2_1_30  = zeros(nbreal,nr);
 
 for k=1:nbreal
  k
  %mode reconstruction using the demodulation procedure
  [snr_demod2_o(k,:,:),snr_demod3_o(k,:,:),snr_demod4_o(k,:,:)]...
                       = Recons_demod(cas,0,ra);
  [snr_demod2_o_10(k,:,:),snr_demod3_o_10(k,:,:),snr_demod4_o_10(k,:,:)]...
                       = Recons_demod(cas,10,ra);
  [snr_demod2_o_30(k,:,:),snr_demod3_o_30(k,:,:),snr_demod4_o_30(k,:,:)]...
                       = Recons_demod(cas,30,ra);
                 
  %mode reconstruction based downsampled STFT using optimized Hamming filter
  %with different downsampling values

  %reconstruction with simple hard-thresholding
  %0 dB
  [SNR_modes1_0(k,:),~,coeff_util1_0(k,:)] = reconstruct_modes(cas,0,0,32,ra);
  [SNR_modes2_0(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,0,16,ra);
  if (cas == 3)
   [SNR_modes3_0(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,0,8,ra); 
   [SNR_modes4_0(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,0,4,ra); 
  end
  
  %10 dB
  [SNR_modes1_10(k,:),~,coeff_util1_10(k,:)] = reconstruct_modes(cas,0,10,32,ra);
  [SNR_modes2_10(k,:),~,coeff_util2_10(k,:)] = reconstruct_modes(cas,0,10,16,ra);
  if (cas == 3)
   [SNR_modes3_10(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,10,8,ra); 
   [SNR_modes3_10(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,10,4,ra); 
  end     

  %30 dB
  [SNR_modes1_30(k,:),~,coeff_util1_30(k,:)] = reconstruct_modes(cas,0,30,32,ra);
  [SNR_modes2_30(k,:),~,coeff_util2_30(k,:)] = reconstruct_modes(cas,0,30,16,ra);
  if (cas == 3)
   [SNR_modes3_30(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,30,8,ra); 
   [SNR_modes4_30(k,:),~,coeff_util2_0(k,:)] = reconstruct_modes(cas,0,30,4,ra); 
  end     

  %reconstruction using the variant of hard-thresholding
  %0 dB
  [SNR_modes1_1_0(k,:),~,coeff_util1_1_0(k,:)] = reconstruct_modes(cas,1,0,32,ra);
  [SNR_modes2_1_0(k,:),~,coeff_util2_1_0(k,:)] = reconstruct_modes(cas,1,0,16,ra);
 
  %10 dB
  [SNR_modes1_1_10(k,:),~,coeff_util1_1_10(k,:)] = reconstruct_modes(cas,1,10,32,ra);
  [SNR_modes2_1_10(k,:),~,coeff_util2_1_10(k,:)] = reconstruct_modes(cas,1,10,16,ra);
 
  %30 dB
  [SNR_modes1_1_30(k,:),~,coeff_util1_1_30(k,:)] = reconstruct_modes(cas,1,30,32,ra);
  [SNR_modes2_1_30(k,:),~,coeff_util2_1_30(k,:)] = reconstruct_modes(cas,1,30,16,ra);
 end
 
 X2_o = zeros(nr,len_d);
 X3_o = zeros(nr,len_d);
 X4_o = zeros(nr,len_d);
 
 X2_o(:,:) = mean(snr_demod2_o);
 X3_o(:,:) = mean(snr_demod3_o);
 X4_o(:,:) = mean(snr_demod4_o);

 Xdown1   = zeros(nr,1);
 Xdown2   = zeros(nr,1);
 
 Xdown1(:,:)   = mean(SNR_modes1_0);
 Xdown2(:,:)   = mean(SNR_modes2_0);
 if (cas == 3)
  Xdown3   = zeros(nr,1);  
  Xdown3(:,:)   = mean(SNR_modes3_0);
  Xdown4   = zeros(nr,1);  
  Xdown4(:,:) = mean(SNR_modes4_0);
 end
 
 Xdown1_1   = zeros(nr,1);
 Xdown2_1   = zeros(nr,1);
 
 Xdown1_1(:,:)   = mean(SNR_modes1_1_0);
 Xdown2_1(:,:)   = mean(SNR_modes2_1_0);

 if (cas == 1)
  %we only compute the demodulation with second order synchrosqueezing, higher orders 
  %synchrosqueezing transforms lead to the same results
  
  %SNR = 0 dB 
  
  figure()
  plot(d,X2_o(1,:),...
       d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
       d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
       'linewidth',2,'markersize',20);
  hold on;
  plot(d,X2_o(2,:),'--',...
       d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
       d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
       'linewidth',2,'markersize',20);
  xlabel('d');
  ylabel('output SNR');
  legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');
      
  set(gca,'fontsize',24);
  hold off;
 
   %SNR = 10 dB
   figure()
   X2_o(:,:)   = mean(snr_demod2_o_10);
   Xdown1(:,:) = mean(SNR_modes1_10);
   Xdown2(:,:) = mean(SNR_modes2_10);
  
   Xdown1_1(:,:) = mean(SNR_modes1_1_10);
   Xdown2_1(:,:) = mean(SNR_modes2_1_10);
  
   plot(d,X2_o(1,:),...
       d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
       d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
       'linewidth',2,'markersize',20);
   hold on;
   plot(d,X2_o(2,:),'--',...
       d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
       d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
       'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');

   set(gca,'fontsize',24);
   hold off;
   
   %SNR = 30 dB
   figure()
   X2_o(:,:)   = mean(snr_demod2_o_30);
   Xdown1(:,:) = mean(SNR_modes1_30);
   Xdown2(:,:) = mean(SNR_modes2_30);
   
   Xdown1_1(:,:) = mean(SNR_modes1_1_30);
   Xdown2_1(:,:) = mean(SNR_modes2_1_30);
 
   plot(d,X2_o(1,:),...
       d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
       d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
       'linewidth',2,'markersize',20);
   hold on;
   plot(d,X2_o(2,:),'--',...
       d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
       d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
       'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');
   set(gca,'fontsize',24);
   hold off;
  elseif (cas == 2),
   % again we do not consider SST3-demod and SST4-demod
 
   % SNR = 0 dB
   figure()
   plot(d,X2_o(1,:),...
       d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
       d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
       'linewidth',2,'markersize',20);
   hold on;
   plot(d,X2_o(2,:),'--',...
        d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
        d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
       'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');

   set(gca,'fontsize',24);
   hold off;
 
   % SNR = 10 dB
   figure()
   X2_o(:,:)   = mean(snr_demod2_o_10);
   Xdown1(:,:) = mean(SNR_modes1_10);
   Xdown2(:,:) = mean(SNR_modes2_10);
   
   Xdown1_1(:,:) = mean(SNR_modes1_1_10);
   Xdown2_1(:,:) = mean(SNR_modes2_1_10);
  
   plot(d,X2_o(1,:),...
        d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
        d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
        'linewidth',2,'markersize',20);
   hold on;
   plot(d,X2_o(2,:),'--',...
        d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
        d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
       'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');

   set(gca,'fontsize',24);
   hold off;
   
   % SNR = 30 dB
 
   X2_o(:,:)   = mean(snr_demod2_o_30);
   Xdown1(:,:) = mean(SNR_modes1_30);
   Xdown2(:,:) = mean(SNR_modes2_30);
 
   Xdown1_1(:,:) = mean(SNR_modes1_1_30);
   Xdown2_1(:,:) = mean(SNR_modes2_1_30);
 
   figure()
   
   plot(d,X2_o(1,:),...
       d,Xdown1(1)*ones(1,len_d),'-<',d,Xdown2(1)*ones(1,len_d),'->',... 
       d,Xdown1_1(1)*ones(1,len_d),'--<',d,Xdown2_1(1)*ones(1,len_d),'-->',...
       'linewidth',2,'markersize',20);
   hold on;
   plot(d,X2_o(2,:),'--',...
        d,Xdown1(2)*ones(1,len_d),'-d',d,Xdown2(2)*ones(1,len_d),'-o',...
        d,Xdown1_1(2)*ones(1,len_d),'--d',d,Xdown2_1(2)*ones(1,len_d),'--o',...
        'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod, f_1$','$STFT-M_1, R = 32, f_1$','$STFT-M_1, R = 16, f_1$',...
      '$STFT-M_2, R = 32, f_1$','$STFT-M_2, R = 16, f_1$',...
      '$FSST2-demod, f_2$','$STFT-M_1, R = 32, f_2$','$STFT-M_1, R = 16, f_2$',...
      '$STFT-M_2, R = 32, f_2$','$STFT-M_2, R = 16, f_2$'},'Interpreter','latex');

   set(gca,'fontsize',24);
   hold off;
 
 elseif (cas == 3),
   % In that case, we consider SST3-demod and SST4-demod
 
   % SNR = 0 dB
   figure()
   plot(d,X2_o(1,:),d,X3_o(1,:),'--',d,X4_o(1,:),'-.',...
        d,Xdown2(1)*ones(1,len_d),'->',...
        d,Xdown3(1)*ones(1,len_d),'-o', d,Xdown4(1)*ones(1,len_d),'-d',...
        'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod$','$FSST3-demod$','$FSST4-demod$',...
       '$STFT-M_1, R = 16$','$STFT-M_1, R = 8$','$STFT-M_1, R = 4$'},...
       'Interpreter','latex');
   set(gca,'fontsize',24);
  
   %SNR = 10 dB

   figure()
   X2_o(:,:) = mean(snr_demod2_o_10);
   X3_o(:,:) = mean(snr_demod3_o_10);
   X4_o(:,:) = mean(snr_demod4_o_10);

   Xdown1(:,:)   = mean(SNR_modes1_10);
   Xdown2(:,:)   = mean(SNR_modes2_10);
   Xdown3(:,:)   = mean(SNR_modes3_10);
   Xdown4(:,:)   = mean(SNR_modes3_10);

   plot(d,X2_o(1,:),d,X3_o(1,:),'--',d,X4_o(1,:),'-.',...
        d,Xdown2(1)*ones(1,len_d),'->',...
        d,Xdown3(1)*ones(1,len_d),'-o',d,Xdown4(1)*ones(1,len_d),'-d', 'linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod$','$FSST3-demod$','$FSST4-demod$',...
       '$STFT-M_1, R = 16$','$STFT-M_1, R = 8$','$STFT-M_1, R = 4$'},...
       'Interpreter','latex');

   set(gca,'fontsize',24);

   
   %SNR = 30 dB
   figure()
   X2_o(:,:)   = mean(snr_demod2_o_30);
   X3_o(:,:)   = mean(snr_demod3_o_30);
   X4_o(:,:)   = mean(snr_demod4_o_30);
   Xdown1(:,:) = mean(SNR_modes1_30);
   Xdown2(:,:) = mean(SNR_modes2_30);
   Xdown3(:,:) = mean(SNR_modes3_30);
   Xdown4(:,:) = mean(SNR_modes4_30);
   
   plot(d,X2_o(1,:),d,X3_o(1,:),'--',d,X4_o(1,:),'-.',...
        d,Xdown2(1)*ones(1,len_d),'->',...
        d,Xdown3(1)*ones(1,len_d),'-o',d,Xdown4(1)*ones(1,len_d),'-d','linewidth',2,'markersize',20);
   xlabel('d');
   ylabel('output SNR');
   legend({'$FSST2-demod$','$FSST3-demod$','$FSST4-demod$',...
       '$STFT-M_1, R = 16$','$STFT-M_1, R = 8$','$STFT-M_1, R = 4$'},...
       'Interpreter','latex');

   set(gca,'fontsize',24);
 end    
end

