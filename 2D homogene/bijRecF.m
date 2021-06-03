function [i,j]=bijRecF(k,nx)
  j = floor((k-1)/nx);
  i = k - nx*j -1;
endfunction
