function F=EvalSource2DVec(f,x,y,nx,ny,t)
  N=nx*ny;
  [I,J]=bijRecF(1:N,nx);
  X=x(I+1);
  Y=y(J+1);
  F=(f(t,X,Y))';
endfunction