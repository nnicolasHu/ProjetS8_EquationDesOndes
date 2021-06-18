function Axy=Lap2DAssembling(nx,ny,hx,hy)
  Jnx=sparse(1:nx,1:nx,ones(1,nx),nx,nx);
  Jny=sparse(1:ny,1:ny,ones(1,ny),ny,ny);
  Ax=(1/(hx*hx))*spMatTriDiagVec([ones(1,nx-1)],[-2*ones(1,nx)],[ones(1,nx-1)]);
  Ay=(1/(hy*hy))*spMatTriDiagVec([ones(1,ny-1)],[-2*ones(1,ny)],[ones(1,ny-1)]);
  Axy=spMatKron(Jny,Ax)+spMatKron(Ay,Jnx);
endfunction
