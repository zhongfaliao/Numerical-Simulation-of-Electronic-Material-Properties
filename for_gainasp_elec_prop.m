clc;
clear all;
% syms y;

y =0:0.01:1;
x =0:0.01:1;
z = NaN(length(x),length(y));
for m =1:1:length(x);
    for n =1:1:length(y);
        z(m,n)=1.35 +0.668*x(m) -1.068*y(n) +0.758*x(m)^2 +0.078*y(n)^2 -0.069*x(m)*y(n) -0.332*x(m)^2*y(n) +0.03*x(m)*y(n)^2;
    end
end
surf(y,x,z)

%  1.35 + (0.758*x + 0.642)*x + (0.101*y - 1.101)*y-(0.28*x - 0.109*y + 0.159)*x*y