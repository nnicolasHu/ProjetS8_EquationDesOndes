function K=Lap1D(d)
  I=[2:d 1:d 1:d-1];
  J=[1:d-1 1:d 2:d];
  K=sparse(I,J,[ones(1,d-1) -2*ones(1,d) ones(1,d-1)],d,d);
endfunction
