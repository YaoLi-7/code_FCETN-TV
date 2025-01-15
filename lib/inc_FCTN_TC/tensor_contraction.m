
function Out = tensor_contraction(X,Y,n,m)
Lx = size(X);      Ly = size(Y);   
Nx = ndims(X);     Ny = ndims(Y);
indexx = 1:Nx;     indexy = 1:Ny;
indexx(n) = [];    indexy(m) = [];

tempX = permute(X,[indexx,n]);  tempXX=reshape(tempX,prod(Lx(indexx)),prod(Lx(n)));
tempY = permute(Y,[m,indexy]);  tempYY=reshape(tempY,prod(Ly(m)),prod(Ly(indexy)));
tempOut = tempXX*tempYY;
Out     = reshape(tempOut,[Lx(indexx),Ly(indexy)]);
