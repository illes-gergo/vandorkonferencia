%LN kirstály fázistörésmutatója

function [ nTHz ] = nTHzo( omega, T)
if T == 100
    nTHz = 4.69732-0.03006e-12*omega/2/pi+0.04066e-24*(omega/2/pi).^2;
elseif T == 300
    nTHz = 4.91372-0.01747e-12*omega/2/pi+0.04004e-24*(omega/2/pi).^2;
end;
end

