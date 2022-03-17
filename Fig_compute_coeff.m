function [coeff_util0,coeff_util10,coeff_util30] = Fig_compute_coeff(cas)
 
 %cas: corresponds to the type of signal we are studying
 %ra : portion between Nfft and N
 
 nbreal = 5;

 if (cas <= 2)
  nr = 2;
 else
  nr = 1;
 end
 
 %the number of coefficients used for the reconstruction of each mode
 %for reconstrcution based on STFT (first technique)
 coeff_util_0  = zeros(nbreal,nr);
 coeff_util_10  = zeros(nbreal,nr);
 coeff_util_30  = zeros(nbreal,nr);
  
 for k=1:nbreal
  k
  %reconstruction using the variant of hard-thresholding
  if (cas <= 2)
   [~,~,coeff_util_0(k,:)]  = reconstruct_modes(cas,1,0,16,4);
   [~,~,coeff_util_10(k,:)] = reconstruct_modes(cas,1,10,16,4);
   [~,~,coeff_util_30(k,:)] = reconstruct_modes(cas,1,30,16,4);
  else
   [~,~,coeff_util_0(k,:)]  = reconstruct_modes(cas,1,0,8,4);
   [~,~,coeff_util_10(k,:)] = reconstruct_modes(cas,1,10,8,4);
   [~,~,coeff_util_30(k,:)] = reconstruct_modes(cas,1,30,8,4);
  end
 end 
 coeff_util0 = mean(coeff_util_0); 
 coeff_util10 = mean(coeff_util_10);
 coeff_util30 = mean(coeff_util_30);
end