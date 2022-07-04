function out = V5_fgv(z,A_kompozit,omega,T,k_omega,k_OMEGA,khi_eff,dnu,domega,omega0,dt,k_omega0,k_OMEGA0,beta4,simp,dz,w0)

% T = 100;    %K
c = 3e8;    %m/s
hv = 6.626e-34/2/pi; %J*s
%beta2 = 0.05e-13; %m^2/W
e0 = 8.8541878e-12; %As/Vm
% lambda0 = 1031.8e-9;   %m
% N = 2*1e4;    %db
% tau = 150e-15;  %s
% I0 = 20e1/tau;    %GW/cm^2
% khi_eff =  360e-12; %pm/V;
% e0 = 8.854e-12;  %F*m^2
% nu0 = 0.5e12;
%
% omega0 = 2*pi*c/lambda0;
% omegaMAX = 5e14*2*pi;
% domega = omegaMAX/N;
% dnu = domega/2/pi;
%
% omega = (0:N-1)*domega;
%
% deltaOmega =2*sqrt(2*log(2))/tau;
%
% lambda = 2*pi*c./omega;
% lambda(1) = lambda(2);
% ngp0 = ngp(lambda0,T);
% np0 = neo(lambda0,T);
%
% nTHz = nTHzo(lambda0,T);
% vfTHz = c./nTHz;
%
% gamma = acos(ngp0/nTHzo(2*pi*nu0,T));
%
% A0 = sqrt(2*I0/neo(lambda0,T)/e0/c)*tau*sqrt(pi/log(2)); % ???
% A0 = sqrt(2*I0/neo(lambda0,T)/e0/c)*tau/(2*sqrt(2*pi*log(2)));
% Aop = A0*exp(-((omega-omega0).^2/deltaOmega.^2));
% plot(Aop(:,1));
% return;
% n_omega = neo(lambda,T);
% k_OMEGA = real(omega.*nTHzo(omega,T)/c);%+1e5;
% ddk_omega = -ngp0.^2/omega0/c/np0*tan(gamma)^2;
% k_omega = real(1/cos(gamma).*(omega.*n_omega/c+(omega-omega0).^2/2.*ddk_omega));%+1e5;
%
%
%
% if z_index == 1
%

% end

ATHz = A_kompozit(1,:,1);
Aop = A_kompozit(1,:,2);

%TPA
global Nc;
global kdz;
%Nc = [Nc beta2*e0^2*c^2*neo(2*pi*c/omega0)^2*sum(abs(ifft(Aop.*exp(-1i*k_omega*z)*2*pi/dt)).^4)/2/h/omega0*dt]
k_OMEGA = real(omega.*nTHzo(omega,Nc(end))/c);


It = e0/2*c*neo(2*pi*c/omega0)*abs(ifftshift(ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z)*2*pi/dt))).^2;
Nt = beta4*cumsum(It.^4)*dt/4/hv/omega0;
%At = ifftshift(ifft(Aop.*exp(-1i*(1*(k_omega-k_omega0)*z))))*2*pi/dt;
ITHzt = abs(ifftshift(ifft(ATHz.*exp(-1i*(kdz+k_OMEGA*dz)+1i*k_OMEGA0*z)))).^2;
ITHzt = ITHzt/max(ITHzt);
THzint = sum(ITHzt);
%ITHzint = cumsum(ITHzt);
%ITHzint = ITHzint/max(ITHzint);
%if ATHzint(end)>0
%    ATHzint = ATHzint/max(ATHzint);
%    Neff = (1-ATHzint).*Nt;
%else
%    Neff = Nc(end);
%end;
if THzint>0
    Neff = sum(Nt.*ITHzt)./THzint;
else
    Neff = Nc(end);
end;
%ATHzt = ATHzt/max(ATHzt);
%[a i1] = min(abs(ITHzint-exp(-2)));
%[b i2] = min(abs(ITHzint-(1-exp(-2))));
%if i2>i1
%    N_ = sum(Neff(i1:i2))/(i2-i1);
%else
%    N_ = Nc(end);
%end;
%return;
%THz_dur = dt*(i2-i1);
%sulyozott = sum(ATHzt(i1:i2).*It(i1:i2).^4)/(i2-i1)
Nc = [Nc Neff+8e20];

%At2 = ifftshift(ifft(Aop))*omegaMAX;
%    It2 = e0/2*c*abs(At2).^2;
%    Aop= Aop-fft(fftshift(beta4/2*It2.^3.*At2*dz))/omegaMAX;

%MFA = 1*fft(fftshift(It.^3.*At.*beta4/2))*dt/2/pi.*exp(1i*k_omega*z)+0*1i*k_omega.*Aop;

%Nc = [Nc beta4*e0^4*c^4*neo(2*pi*c/omega0)^4/2^4*sum(abs(ifft(Aop.*exp(-1i*k_omega*z)*2*pi/dt)).^8)/4/hv/omega0*dt]
%Imax = e0*c*neo(2*pi*c/omega0)*max(abs(ifft(Aop.*exp(-1i*k_omega*z)*2*pi/dt)).^2)/2
%Nc = [Nc beta4*Imax.^4/4/hv/omega0*150e-15];

%k_OMEGA = real(omega.*nTHzo(omega,Nc(end))/c);

abszorpcio = 2*omega/c.*imag(sqrt(er(omega,Nc(end))));

temp1 = zeros(size(ATHz));
[~,I] = max(abs(Aop));
for nagy_omega = 2:ceil(10e12/dnu)
    temp1(nagy_omega) = -1*abszorpcio(nagy_omega)/2*ATHz(nagy_omega)-1i*khi_eff*omega(nagy_omega).^2/2/c^2/k_OMEGA(nagy_omega)...
        *sum(Aop(nagy_omega:end-1).*conj(Aop(1:end-nagy_omega))...
        .*exp(-1i*(k_omega(nagy_omega:end-1)-k_omega(1:end-nagy_omega)-k_OMEGA(nagy_omega))*z))*domega;
end;
temp2 = zeros(size(temp1));
%/lambda;

for kis_omega = I-simp:I+simp
    zR = k_omega(kis_omega)*w0^2/2;%pi*w0^2*omega(kis_omega)/2/pi/c;
    q = sqrt(1+z^2/zR^2);
    
    %eA = 1*(1i*k_omega(kis_omega)*z+(3/4*z^2-zR^2/2)/(z^2+zR^2))/(2*1i*k_omega(kis_omega)*(z^2+zR^2)-z);
    eA = (5*z^2+4*1i*k_omega(kis_omega)*z*q^2*zR^2-4*q^2*zR^2)/...
        (4*z*q^2*zR^2+8*1i*k_omega(kis_omega)*q^4*zR^4);
    %eP = 1*(1+z^2/zR^2)^(1/4)/(2*1i*k_omega(kis_omega)+z/(z^2+zR^2));
    temp2(kis_omega) = -1*eA*Aop(kis_omega)...
        -1*1i*khi_eff*omega(kis_omega).^2/c^2/2/k_omega(kis_omega)*(sum(Aop(kis_omega:end-1).*conj(ATHz(1:end-kis_omega))...
        .*exp(-1i*(k_omega(kis_omega:end-1)-k_omega(kis_omega)-k_OMEGA(1:end-kis_omega))*z))...
        +sum(Aop(kis_omega:-1:1).*(ATHz(1:kis_omega))...
        .*exp(-1i*(k_omega(kis_omega:-1:1)-k_omega(kis_omega)+k_OMEGA(1:kis_omega))*z)))*domega;
%temp2(kis_omega) = 0;
end;

%tpa


out = zeros(1,length(omega),2);
out(1,:,1) = temp1;
out(1,:,2) = temp2;
end