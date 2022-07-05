%% tisztogatÃ¡s, kÃ¶nyvtÃ¡r kivÃ¡lasztÃ¡s
clc;clear all;close all;


%workdir = matlab.desktop.editor.getActiveFilename;
%while workdir(end) ~= '\'
%    workdir = workdir(1:end-1);
%end
%cd(workdir(1:end-1));-


tic;
dir_n = strrep(strrep(datestr(datetime), ' ', '_'),':','.');
%mkdir(dir_n);
%copyfile('initial.xlsx', dir_n);

%D = xlsread('initial.xlsx');

fluence = 20*1e1; %mj/cm^2
tau = 200*1e-15;
I0 = fluence./tau;
lambda0 = 1030*1e-9;
elochirp = 0*1e-3;
nu0 = 0.5*1e12;
z_vegso = 6*1e-3;
T = 300;    %K
beta4 = 0*2.6e-4*1e-37;
dz = 10*1e-6;
omegaMAX = 2*pi*1e12*800;
N = 2e4;
simp = 700;
utem = 20;
w0 = 20*3.0404e-6;%30e-6; %beamlet szélessége

c = 3e8;    %m/s
c0 = 3e8;
khi_eff = 360e-12; %pm/V;
e0 = 8.854e-12;  %F*m^2

deg = pi/180;

global Nc;
Nc = 0*8e20;
global kdz;
omega0 = 2*pi*c/lambda0;

dt = 2*pi/omegaMAX;
t = (0:N-1)*dt;
t = t-t(end)/2;
domega = omegaMAX/N;
dnu = domega/2/pi;

omega = (0:N-1)*domega;
nu=omega/2/pi;

deltaOmega =2*sqrt(2*log(2))/tau;

lambda = 2*pi*c./omega;
lambda(1) = lambda(2);
ngp0 = ngp(lambda0,T);
np0 = neo(lambda0,T);

nTHz0 = nTHzo(2*pi*nu0,T);
vfTHz = c./nTHz0;

nTHz = nTHzo(omega,T);
nTHz(1) = 1;

gamma = acos(ngp0/nTHz0);
w0 = w0*cos(gamma);  %tényleges w0
z_period = w0/sin(gamma);
ngVPHG = 1.5221;
gammaNM = 0*atan(ngVPHG/ngp0*tan(gamma));
GDz = w0*tan(gamma-gammaNM)/cos(gammaNM);
%return;
if 10*dz>z_period
    dz =z_period/10;
    LI = 10;
else
    dz = Lambda/round(z_period/dz);
    LI = fix(Lambda/dz);
end;
%dz = 50e-6;
z = 0:dz:z_vegso;


%return;
n_omega = neo(lambda,T);
%n_THz = nTHzo(omega,0);
k_OMEGA = real(omega.*nTHz/c);%+1e5;
k_OMEGA0 = real(omega.*nTHzo(2*pi*nu0,T)/c);
kdz = k_OMEGA*0;
ddk_omega = -ngp0.^2/omega0/c/np0*tan(gammaNM)^2;
k_omega = real(1/cos(gammaNM).*(omega.*n_omega/c+(omega-omega0).^2/2.*ddk_omega));%+1e5;
k_omega0 = real(1/cos(gamma).*omega.*ngp0/c);

GDz = LI*dz*real(1/cos(gammaNM).*omega*ngp0/c-1/cos(gamma)*omega*ngp0/c);


A0 = sqrt(2*I0/neo(lambda0,T)/e0/c)*tau/(2*sqrt(2*pi*log(2)));
Aop = A0*exp(-((omega-omega0).^2/deltaOmega.^2)).*exp(1i*(k_omega-k_omega0)*elochirp);


Aop0 = Aop;

ATHz = zeros(size(Aop));

%return;

FI = 4*nTHz.^2./(1+nTHz).^2;
FA = 2*nTHz./(1+nTHz);
pF = sum(abs(Aop).^2)*np0;


%% kezdÅ‘ komplex tÃ©rerÅ‘ssÃ©gek Ã¶sszeillesztÃ©se
A_komp(1,:,1) = ATHz;
A_komp(1,:,2) = Aop0;

%% diffegyenlet lÃ©trehozÃ¡sa
v6_fgv =@(z,A_kompozit) diffegy(z,A_kompozit,omega,T,k_omega,k_OMEGA,khi_eff,dnu,domega,omega0,dt,k_omega0,k_OMEGA0,beta4,simp,dz,w0/2);


% fileID = fopen(strcat(dir_n,'\data.txt'),'w');
% fprintf(fileID,'GaAs-ben történõ THz-keltés 1 pumpáló impulzussal beamletekkel\n\n');
% fprintf(fileID,'Pumpa hullámhossza: %6.2f \n',lambda0*1e9);
% fprintf(fileID,'Pumpa döntése: %6.2f° \n',gamma/deg);
% fprintf(fileID,'Beamlet döntése: %6.2f° \n',gammaNM/deg);
% fprintf(fileID,'Pumpa FL intenzitása: %6.2f GW/cm^2\n',I0*1e-13);
% %fprintf(fileID,'Pumpa intenzitása csörpölve: %6.2f GW/cm^2\n',chI0*1e-13);
% fprintf(fileID,'Pumpa FL félértékszélessége: %6.2f fs\n',tau*1e15);
% %fprintf(fileID,'Pumpa csörpölt félértékszélessége: %6.2f fs\n',tau*1e15);
% 
% fprintf(fileID,'Kristály anyaga: GaP\n');
% fprintf(fileID,'Fázisillesztési frekvencia (periodikus polarizációval): %6.2f THz\n',nu0*1e-12);
% %fprintf(fileID,'Periódus hossza: %6.2f um\n',Lambda*1e6);
% fprintf(fileID,'Kristályhossz: %6.2f mm\n',z_vegso*1e3);
% %fprintf(fileID,'Kristály hõmérséklete: %3.0f K\n',T);
% fprintf(fileID,'Pumpa fázis törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0*1e9,np0);
% fprintf(fileID,'Pumpa csoport törésmutatója %6.2f nm-es hullámhosszon: %6.4f\n',lambda0*1e9,ngp0);
% fprintf(fileID,'THz fázis törésmutatója %6.2f THz-es frekvencián: %6.4f\n',nu0*1e-12,nTHz0);
% fprintf(fileID,'Nemlineáris szuszceptibilitás (Chi_eff): %6.2f pm/V\n\n',khi_eff*1e12);
% 
% fprintf(fileID,'Térbeli lépések száma: %d (dz = %f6.2 um)\n',length(z), dz*1e6);
% fprintf(fileID,'Spektrális felbontás: %d (d_nu = %f6.2 GHz)\n',N, dnu*1e-9);
% fprintf(fileID,'Spektrális tartomány: (Nu_max=) %6.2 THz',omegaMAX/2/pi*1e-12);
% fprintf(fileID,'Idõbeli felbontás: %d (d_t = %f6.2 fs)\n',N, dt*1e15);
% fprintf(fileID,'Idõbeli tartomány: (T_max=) %6.2 ps\n\n',N*dt*1e12);
% 
% fprintf(fileID,'Számítás idõpontja: %s\n',dir_n);
% 
% fclose(fileID);
% 
% 
%    % Aop = Aop.*exp(-1i*GDz*ngp0/c0/2*(omega-0*omega0));
% 
%     Aop = Aop.*exp(1i*GDz/2);
A_contain = zeros([length(z),length(omega),2]);

for ii = 1:length(z)
A_komp(1,:,1) = ATHz;
A_komp(1,:,2) = Aop;

%% diffegyenlet megoldÃ¡sa
[z2,A_komp] = RK4_M(v6_fgv,dz,(ii-1)*dz,A_komp,(ii+0.1)*dz);

Nc2(ii) = Nc(end); %lépésenként a szabad töltéshordozók száma
%kdz = kdz%+real(omega.*nTHzo(omega,Nc(end))/c)*dz;
ATHz = A_komp(2,:,1).';
Aop = A_komp(2,:,2).';
clear A_komp;
A_contain(ii,:,1) = ATHz;
A_contain(ii,:,2) = Aop;
effic(ii) = sum(abs(ATHz).^2.*FI.')/pF ;
%% kirajzolÃ¡s
subplot(2,3,1)
plot(lambda,abs(Aop).^2);
xlim([lambda0-200e-9,lambda0+200e-9]);
title('Pumpa intensity spectrum');
xlabel('Frequency (Hz)');

% subplot(2,3,2)
% plot(z(1:ii)*1e3,Nc2);
% title('Nc');
% xlabel('Crystal length (mm)');

abszorpcio = 2*omega/c.*imag(sqrt(er(omega,Nc2(end))));
ab0 = 2*omega/c.*imag(sqrt(er(omega,0)));

subplot(2,3,3)
plot(t*1e12, 1e-13*np0*e0*c0/2*abs(ifftshift(ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2);
xlim([-2 2]);
title('Pump Intensity');
xlabel('time (ps)');
ylabel('Intensity (GW/cm^2)');
%return;
% 
% %subplot(2,4,5);
% %plot(nu*1e-12,[nTHzo(2*pi*nu,Nc2(end)); nTHzo(2*pi*nu,0)]);
% xlim([0,8]);
% title('Instantaneous THz refr. ind.');
% xlabel('Frequency (THz)');
% legend('Nc inst.','Nc = 0');

subplot(2,3,6);
plot(nu*1e-12,[abszorpcio; ab0]/100);
xlim([0,8]);
title('Instantaneous THz absorption [1/cm]');
xlabel('Frequency (THz)');
legend('Nc inst.','Nc = 0');

subplot(2,3,4);
b = plot(z(1:ii)*1e3,effic);
title('Efficiency');
xlabel('Crystal length(mm)');


subplot(2,3,5);
plot(nu*1e-12,abs(ATHz));
xlim([0,8]);
title('THz electric field spectrum');
xlabel('Frequency (THz)');


subplot(2,3,2)
plot(t*1e12,1e-5*real(ifftshift(ifft(FA.'.*ATHz.*exp(-1i*kdz+1i*k_OMEGA0*z(ii)).'))*omegaMAX));
title('THz electric field (out)');
xlabel('Time (ps)');
ylabel('Electric Field (KV/cm)');
xlim([-10, 5])

%return;
drawnow;

Iev(ii) = max(np0*e0*c0/2*abs(ifftshift(ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2);

% if mod(ii,utem)==0
%     %if D(21) == 1    
%         dlmwrite(strcat(dir_n,'\PumpInt-',int2str(fix(ii*dz*1e6)),'_um.txt'),[(t).' np0*e0*c0/2*abs(ifftshift(ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2]);
%     %end;
%     
%     %if D(20) == 1
%         dlmwrite(strcat(dir_n,'\PumpSpec-',int2str(fix(ii*dz*1e6)),'_um.txt'),[lambda.' abs(Aop).^2]);
%     %end;
%     
%     %if D(19) == 1
%         dlmwrite(strcat(dir_n,'\THzt-',int2str(fix(ii*dz*1e6)),'_um.txt'),[t.' real(ifftshift(ifft(FA.'.*ATHz.*exp(-1i*kdz-1i*k_OMEGA0*z(ii)).')))*omegaMAX]);
%     %end;
%     %if D(18) == 1
%         dlmwrite(strcat(dir_n,'\THzSpec-',int2str(fix(ii*dz*1e6)),'_um.txt'),[nu.' abs(ATHz)]);
%     %end;
% end;
% 
 if mod(ii,(LI)) == 0
    Aop = Aop.*exp(-1i*GDz*ngp0/c0/cos(gamma)*(omega-0*omega0).');
    
     Aop = Aop.*exp(1i*GDz.');

 end;

end;

hossz = toc;
c0 = 3e8;

h = fix(hossz/3600);
min = fix((hossz-h*3600)/60);
sec = hossz-h*3600-min*60;
% fileID = fopen(strcat(dir_n,'\data.txt'),'a+');
% fprintf(fileID,'\nCalculation time: %2.0f h %2.0f min %2.0f sec', h, min, sec);
% fclose(fileID);
% 
% 
%     dlmwrite(strcat(dir_n,'\PumpInt.txt'),[(t).' np0*e0*c0/2*abs(ifftshift(ifft(Aop.*exp(-1i*(k_omega-k_omega0)*z(ii)).'))*omegaMAX).^2]);
%     dlmwrite(strcat(dir_n,'\PumpSpec.txt'),[nu.' abs(Aop).^2]);
%     dlmwrite(strcat(dir_n,'\THzt.txt'),[t.' real(ifftshift(ifft(FA.'.*ATHz.*exp(-1i*kdz+1i*k_OMEGA0*z(ii)).')))*omegaMAX]);
%     dlmwrite(strcat(dir_n,'\THzSpec.txt'),[nu.' abs(ATHz)]);
%     dlmwrite(strcat(dir_n,'\Effic.txt'),[z.' effic.']);
%     dlmwrite(strcat(dir_n,'\Nc.txt'),[z.' Nc2.']);
% 
%  clearvars -except JJ Iev;
%  close all;
%  clc;





