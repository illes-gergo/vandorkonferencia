%LN kirstály fázistörésmutatója

function [ nTHz ] = nTHzo(omega,Nc)
%nu = omega/2/pi/3e8*1e-2; %[nu] = 1/cm
%nTHz = real(9.09+2.06*363.4^2./(363.4^2-nu.^2-2*1i*0.55*nu));

nTHz = real(sqrt(er(omega,Nc)));
end

