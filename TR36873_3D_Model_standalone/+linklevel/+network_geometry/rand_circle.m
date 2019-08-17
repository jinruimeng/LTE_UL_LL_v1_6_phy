function pos = rand_circle(N, x, y, r)
% RAND_CIRC(N,x,y,r) generates N random points in a circle with radius r 
% and center at (x,y)

Ns = round(4/pi*N + 2.5*sqrt(N) + 100);
X = rand(Ns,1)*(2*r) - r;
Y = rand(Ns,1)*(2*r) - r;
I = find(sqrt(X.^2 + Y.^2) <= r);
X = X(I(1:N)) + x;
Y = Y(I(1:N)) + y;
pos = [X, Y];