function [ alpha ] = aTHzo( omega, T)
if T == 100
    [a b] = min(abs(omega-1.5708e13));
    alpha(1:b) = 100*(6.20304e-12*omega(1:b)/2/pi-5.03352e-24*(omega(1:b)/2/pi).^2+2.63295e-36*(omega(1:b)/2/pi).^3);
    alpha(b:length(omega)) = 100*(6.20304e-12*omega(b)/2/pi-5.03352e-24*(omega(b)/2/pi).^2+2.63295e-36*(omega(b)/2/pi).^3);    
elseif T == 300
    alpha = 100*(2.16411e-12*omega/2/pi+10.81e-24*(omega/2/pi).^2);
end;
end

