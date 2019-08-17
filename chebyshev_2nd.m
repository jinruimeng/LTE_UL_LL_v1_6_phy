function [ y ] = chebyshev_2nd( x,n )
% Chebyshev polynomial of second kind

if (x==1)
    y = n+1;
elseif (x==-1)
    y = (n+1)*(-1)^n;
else
    y = sin((n+1)*acos(x))/sin(acos(x));
end

end

