[coeff_util0,coeff_util10,coeff_util30]       = Fig_compute_coeff(1); 
[coeff_util0_1,coeff_util10_1,coeff_util30_1] = Fig_compute_coeff(2);
[coeff_util0_2,coeff_util10_2,coeff_util30_2] = Fig_compute_coeff(3);
 
N = 1024;

[sum(coeff_util0/N) 7*2 sum(coeff_util10/N) 9*2 sum(coeff_util30/N) 11*2] 
[sum(coeff_util0_1/N) 7*2 sum(coeff_util10_1/N) 11*2 sum(coeff_util30_1/N) 11*2] 
[sum(coeff_util0_2/N) 11 sum(coeff_util10_2/N) 11 sum(coeff_util30_2/N) 11] 