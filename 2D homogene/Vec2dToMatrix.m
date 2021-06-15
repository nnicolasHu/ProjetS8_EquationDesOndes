function M=Vec2dToMatrix(u,nx,ny)
  M = zeros(nx,ny);
  for i=0:ny-1
    indice = i*nx;
    M(:,i+1) = u(indice+1:indice+nx);
  endfor
  
end