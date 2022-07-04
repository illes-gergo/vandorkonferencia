%relatív dielektromos állandó THz-en
function [er_] = er(omega, Nc)
tsc = 180e-15;
meff = 0.25*9.109e-31;
q = 1.602e-19;
e0 = 8.8541878e-12;

op2 = q^2*Nc/e0/meff;

er_ = 9.09+2.06*363.4^2./(363.4^2-omega.^2*0.01^2/(2*pi*3e8)^2-2*1i*0.55*omega*0.01/(2*pi*3e8))+...
    -op2./(omega.^2+1i*omega/tsc);
end

