clc;clear;close all;
load("vandorkonferencia.mat");
writematrix([z.',effic.'],"hatásfok2");
[M,I] = max(effic)
AOP = ATHz(:,I);
EOPO = AOP.*exp(-1i.*k_OMEGA.'.*z(I));
EOPT = real(ifft(EOPO));
[MM,II] = max(abs(EOPT))
offset = round((II-length(t)/2));
EOPT = circshift(EOPT,-offset)*omegaMAX/1e5;
plot(t,EOPT);
writematrix([t.'-max(t)/2,EOPT],"TPF_tere")

clear;

load("szamitas.mat")
ATHz = A_contain(:,:,1).';
[M,I] = max(effic)
AOP = ATHz(:,I);
EOPO = AOP.*exp(-1i.*k_OMEGA.'.*z(I));
EOPT = real(ifft(EOPO));
[MM,II] = max(abs(EOPT))
offset = round((II-length(t)/2));
EOPT = circshift(EOPT,-offset)*omegaMAX/1e5;
plot(t,EOPT);
writematrix([t.',EOPT],"echelon_tere")