function C=spMatKron(A,B)
  n=size(A,1);
  m=size(B,1);
  C=sparse(n*m,n*m);
  
  [I,J,K]=find(A);
  for k=1:length(K)
    C((I(k)-1)*m+1:I(k)*m,(J(k)-1)*m+1:J(k)*m)=K(k)*B;
  endfor
endfunction
