clear all;
close all;
format compact;
clc;

%Global constants
K = 1.38 * 10^(-23);            %Boltzmann constant
T = 273.15;                     %0 degrees Celsius in absolute scale 
q = 1.6 * 10^(-19);             %electron charge

%Catalog at STC (Standard Test Conditions)
Pmpr = 270;
Vmpr = 31.1;
Impr = 8.67;
NOCT = 45;                      %Normal Operating Cell Temperature
miu_Voc = -0.1184;
miu_Isc = 0.0049;
Vocr = 38.2;
Iscr = 9.19;

%-------------------Environment variables-------------------------------%
%ta = 20;                        %atmospheric temperature (20 for NOCT)
tc = 65;  %<---                 %cell temperature (25 for STC)
Tc = T+tc;                      %absolute cell temperature
Tr = T + 25;                    %absolute reference temperature
G = 1000; %<---                 %irradiance (1000 for STC, 800 for NOCT)
Gr = 1000;                      %reference irradiance
Vtr = K*Tc/q;                   %thermal voltage equivalent
%-----------------------------------------------------------------------%

%---------------------Circuit variables---------------------------------%
Vo = 100;                       %Load voltage
dIl = 0.2;                      %delta I_L
%-----------------------------------------------------------------------%


%Parameters of the mathematical model (m, Rsh, Rs, Io, Is)
syms  m Rsh Rs
V = [m, Rsh, Rs];
F = @(V) [Impr - Iscr + (Vmpr + V(3)*Impr + V(3)*Iscr)/ V(2) + (Iscr - (Vocr - V(3)*Iscr)/ V(2)) * exp((Vmpr + V(3)*Impr - Vocr)/ (V(1)*Vtr));
            Impr + ((-(V(2)*Iscr - Vocr + V(3)*Iscr) * exp((Vmpr + V(3)*Impr - Vocr)/ (V(1)*Vtr))/ (V(2) * V(1)*Vtr) - 1/V(2))/ (1 + V(3) * (V(2)*Iscr - Vocr + V(3)*Iscr) * exp((Vmpr + V(3)*Impr - Vocr)/ (V(1)*Vtr))/ (V(2) * V(1)*Vtr) + V(3)/V(2))) * Vmpr;
            1/V(2) + (-(V(2)*Iscr - Vocr + V(3)*Iscr) * exp((V(3)*Iscr - Vocr)/ (V(1)*Vtr))/ (V(2) * V(1)*Vtr) - 1/V(2))/ (1 + V(3) * (V(2)*Iscr - Vocr + V(3)*Iscr) * exp((V(3)*Iscr - Vocr)/ (V(1)*Vtr))/ (V(2) * V(1)*Vtr) + V(3)/V(2))];

InitialGuess = [60; 7000; 0.2];
% options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt,'MaxFunctionEvaluations',2000);
options = optimoptions('fsolve','MaxFunctionEvaluations',2000);
sol = fsolve(F, InitialGuess,options);
% ShouldBeZero = F(sol)
m = sol(1);
Rsh = sol(2);
Rs = sol(3);
Ior = (Iscr - (Vocr - Rs*Iscr)/ Rsh) * exp(-Vocr/ (m*Vtr));
Isr = Ior * exp(Vocr/ (m*Vtr)) + Vocr/Rsh;


%Simulation of the model
%Tc = Tc;ta + G*(NOCT-20)/ 800 + 273.15;
Vt = K*Tc/ q;
Isc = G/ Gr * (Iscr + miu_Isc*(Tc - Tr));
Voc = Vocr + miu_Voc*(Tc - Tr) + m*Vt*log(G/Gr);
Io = (Isc - (Voc - Rs*Isc)/ Rsh) * exp(-Voc/ (m*Vt));
Is = Io * exp(Voc/ (m*Vt)) + Voc/ Rsh;

Vd_vetor = linspace(0,50);      %Adjust max value for better plots
I = Is - Io * (exp(Vd_vetor/ (m*Vt)) - 1) - Vd_vetor/ Rsh;
V = Vd_vetor - Rs*I;
P = V .* I;
[Pmpp, ind] = max(P);           %Returns the maximum power and where it occurs
Imax = I(ind);
Vmax = V(ind);

%MATLAB results of the simulation
% figure
% plot(V,I)
% axis([0 Voc*1.1 0 Isc*1.1])
% xlabel('Voltage [V]')
% ylabel('Current [I]')
% 
% figure
% plot(V,P)
% axis([0 Voc*1.1 0 Pmpp*1.1])
% xlabel('Voltage [V]')
% ylabel('Power [W]')


%Simulink results of the simulation
Load = Pmpp/I(ind)^2;
Imax = I(ind);
sim('PVArray.slx')




