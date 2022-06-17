clc;
clear;
close all;
%========================= Initializing Parameters ========================
I0_alpha = 6e-6 ;
k = 0.9;
i_dark = 0.6e-9;
Vsn = 36e-3 ;
Ip = k * (i_dark + I0_alpha) ;
Vth = 5e-3 ;
Vin = 1.58e-6 ;
Ven = 78e-6 ;
Rf = 500e3 ;
% Vo = Rf * I  + Vc ;
L = 40e3 ;
c = 3e8 ;
n = 1.444 ;
Fs = 2e4 ;
t = 0:1/Fs:1 ;
size = length(t) ;
num_samp = 10 ;
del_est = zeros(1,num_samp);
D_vals = zeros(1,num_samp) ;
f = @(x) (0.5*pi + sin(2*pi*300.*x).*exp(-x./0.05)); %.*exp(-x./0.05)
for i=1:num_samp
    x = zeros(1,size)+randi(L);
    D = ((2*L-x)/(c/n))- (x/(c/n)) ;
%     disp(x(1,10));
    %====================== Defining the Noises ======================
    %Shot Noise
    shot_noise1 = (poissrnd(Vsn*1e3,1,size)*5e-7 -0.00002)+ Vsn;
    shot_noise2 = (poissrnd(Vsn*1e3,1,size)*5e-7 -0.00002)+ Vsn;

    %Thermal Noise
    Vth1 = normrnd(0,Vth, 1,size);
    Vth2 = normrnd(0,Vth, 1,size);

    %Relative Intensity Noise
    rltv_noise = normrnd(0,0.1e-9, 1,size);

    n_visible1 = normrnd(0,1e-1, 1,size);
    n_visible2 = normrnd(0,1e-1, 1,size);

    % Detecting at Reveiver
    i1_c = normrnd(0.6e-9,1e-10, 1,size);
    i2_c = normrnd(0.6e-9,1e-10, 1,size);

    I1 = Ip .* (1+rltv_noise).*(1+n_visible1.*cos(f(t))) + i1_c;
    I2 = Ip .* (1+rltv_noise).*(1+n_visible2.*cos(f(t-D))) + i2_c;
    
    Vin1 = normrnd(0,Vin,1,size);
    Vin2 = normrnd(0,Vin,1,size);
    
    Ven1 = normrnd(0,Ven,1,size) ;
    Ven2 = normrnd(0,Ven,1,size) ;
    
    V1_t = Rf * I1 + sqrt(Ven1.^2 + Vin1.^2 + Vth1.^2 + shot_noise1.^2) ;
    V2_t = Rf * I2 + sqrt(Ven2.^2 + Vin2.^2 + Vth2.^2 + shot_noise2.^2) ;
    
    %======================== CC TDE =======================
    [Rvv,lag] = xcorr(V1_t,V2_t , 'biased');
    [~,In] = max(abs(Rvv));
    del_est(1,i) = abs(lag(In)/Fs -mean(D)) ;  
%     disp(L-del_est(1,i)*c/(2*n)) ;
    D_vals(1,i) = D(10) ;
end
figure(1); plot(t,V1_t,t,V2_t);grid on;xlim([0,0.01]);
fprintf('Mean of Time Delay Estimation Error = %f \n',mean(del_est)) ;
fprintf('Standard Deviation of TDE Error = %f \n' ,std(del_est)) ;
% L-mean(del_est)*c/(2*n)
% D_vals(1,:)