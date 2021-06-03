function F=EvalFun2DVec(f,x,y,nx,ny)
  N=nx*ny;
  [I,J]=bijRecF(1:N,nx);
  X=x(I+1);
  Y=y(J+1);
  F=f(X,Y);
endfunction