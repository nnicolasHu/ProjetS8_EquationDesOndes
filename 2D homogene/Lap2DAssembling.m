function Axy=Lap2DAssembling(nx,ny,hx,hy)
  Jnx=sparse(2:nx-1,2:nx-1,ones(1,nx-2),nx,nx);
  Jny=sparse(2:ny-1,2:ny-1,ones(1,ny-2),ny,ny);
  Ax=(1/(hx*hx))*spMatTriDiagVec([ones(1,nx-2) 0],[0 -2*ones(1,nx-2) 0],[0 ones(1,nx-2)]);
  Ay=(1/(hy*hy))*spMatTriDiagVec([ones(1,ny-2) 0],[0 -2*ones(1,ny-2) 0],[0 ones(1,ny-2)]);
  Axy=spMatKron(Jny,Ax)+spMatKron(Ay,Jnx);

endfunction
