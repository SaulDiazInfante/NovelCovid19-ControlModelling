function [y] = phiab(a,b,k,x)
y = phia(a,k,x) * phib(b,k,x);
end

