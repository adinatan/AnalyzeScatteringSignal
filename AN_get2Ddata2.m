function d = AN_get2Ddata2(x, y, vec,N)
%This function loads an input of the format [x y int] and
%generates a matrix of int spanned by x and y. 
 
%Resetting zeros
x = x - (min(x(:)));
y = y - (min(y(:)));


%Dividing by N
x = ceil(x/N)+1;
y = ceil(y/N)+1;

 

minx = min(x(:));
miny = min(y(:));
maxx = max(x(:));
maxy = max(y(:));


  sd = [maxx-minx+1 maxy-miny+1];
 
  d = zeros(sd);

  for n=1:numel(vec)
  d(x(n),y(n) )=vec(n);
  end


