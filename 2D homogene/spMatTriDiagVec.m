function A=spMatTriDiagVec(u,v,w)
  d=length(v);
  I=[2:d 1:d 1:d-1];
  J=[1:d-1 1:d 2:d];
  A=sparse(I,J,[u v w],d,d);
endfunction
