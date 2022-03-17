function [sp1_s,integ1] = demod_multi_omega(s,VSST,omega2,nr,jump)   

 %computation of the different ridges
 [Cs21,~] = exridge_mult(VSST,nr,0,0,jump);

 A = size(VSST);
 
 integ1 = zeros(size(Cs21));
 
 sp1_s  =  zeros(nr,A(2));
 
 t = (0:A(2)-1)/A(2);
 
 for k = 1:nr
  YY = zeros(1,A(2));
  for kk = 1: A(2)
   YY(kk) = omega2(Cs21(k,kk),kk);
  end
  
  integ1(k,:) = cumtrapz(t,YY);
  sp1_s(k,:)  = transpose(s).*exp(-2*1i*pi*(integ1(k,:)-100.*t)); 
 end 
 