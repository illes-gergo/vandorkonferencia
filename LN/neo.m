%LN kristály törésmutatója Sellmeier módszerrel

function [ n ] = neo(lambda,T)
lambda1 = lambda*1e6;
a1 = 5.078;
a2 = 0.0964;
a3 = 0.2065;
a4 = 61.16;
a5 = 10.55;
a6 = 1.59e-2;
b1 = 4.677e-7;
b2 = 7.822e-8;
b3 = -2.653e-8;
b4 = 1.096e-4;
T = T-273;
f = (T-24.5)*(T+570.82);

n = sqrt(a1+b1*f+(a2+b2*f)./(lambda1.^2-(a3+b3*f)^2)+(a4+b4*f)./(lambda1.^2-a5^2)-a6*lambda1.^2);

end

