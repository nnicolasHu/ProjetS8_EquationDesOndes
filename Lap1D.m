function K=Lap1D(d)
  K=-2*eye(d)+diag(ones(1,d-1),1)+diag(ones(1,d-1),-1);
endfunction
