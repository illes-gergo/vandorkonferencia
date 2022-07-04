%LN kristály csoport törésmutatója

function [ ng ] = ngp( lambda )
syms l;
lambda1 = lambda*1e6;
a1 = 1.39;
a2 = 0.172;
b1 = 4.131;
b2 = 0.234;
c1 = 2.57;
c2 = 0.345;
d1 = 2.056;
d2 = 27.52;

n0 = sqrt(1+a1*l.^2./(l.^2-a2^2)+b1*l.^2./(l.^2-b2^2)+c1*l.^2./(l.^2-c2^2)+d1*l.^2./(l.^2-d2^2));
        
        
        a = n0-l*diff(n0);
        l = lambda1;
        ng = eval(a);
end

