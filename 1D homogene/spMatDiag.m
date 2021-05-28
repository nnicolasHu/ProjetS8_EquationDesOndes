function D=spMatDiag(v)
  d=length(v);
  D=sparse(1:d,1:d,v,d,d);
endfunction
