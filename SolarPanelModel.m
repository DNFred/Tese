clear all;
close all;
format compact;
clc;

%Global constants
K = 1.38 * 10^(-23);            %Boltzmann constant
T = 273.15;                     %0 degrees Celsius in absolute scale 
q = 1.6 * 10^(-19);             %electron charge

%Catalog values at STC (Standard Test Conditions)
Pmpr = 270;
Vmpr = 31.1;
Impr = 8.67;
NOCT = 45;                      %Normal Operating Cell Temperature
miu_Voc = -0.1184;
miu_Isc = 0.0049;
Vocr = 38.2;
Iscr = 9.19;
Pmin = 20;

%-------------------Environment variables-------------------------------%
ta = 20;  %<------------        %atmospheric temperature (20 for NOCT)
tc = 25;                        %cell temperature (25 for STC)
Tc = T + tc;                    %absolute cell temperature
Tr = T + 25;                    %absolute reference temperature
G = 1000; %<------------        %irradiance (1000 for STC, 800 for NOCT)
Gr = 1000;                      %reference irradiance
Vtr = K*Tr/q;                   %thermal voltage equivalent
%-----------------------------------------------------------------------%

%---------------------Circuit variables---------------------------------%
N = 1;    %<------------        %Number of panels
Vo = 50;                        %Load voltage
t_PWM = 20e-6;                  %PWM signal period
k_v = 1/ (t_PWM*5000);           %voltage gain for v_ref

timeintegrator1 = 1;
timeintegrator2 = 2;
timeintegrator3 = 3;
timeintegrator4 = 4;
timeintegrator5 = 5;
timeintegrator6 = 6;
timeintegrator7 = 7;
timeintegrator8 = 8;
timeintegrator9 = 9;
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
Tc = ta + G*(NOCT-20)/ 800 + 273.15;
Vt = K*Tc/ q;
Isc = G/ Gr * (Iscr + miu_Isc*(Tc - Tr));
Voc = Vocr + miu_Voc*(Tc - Tr) + m*Vt*log(G/Gr);
Io = (Isc - (Voc - Rs*Isc)/ Rsh) * exp(-Voc/ (m*Vt));
Is = Io * exp(Voc/ (m*Vt)) + Voc/ Rsh;

%Plot of the model
Vd_vetor = linspace(0,50);      %Adjust max value for better plots
I = Is - Io * (exp(Vd_vetor/ (m*Vt)) - 1) - Vd_vetor/ Rsh;
V = Vd_vetor - Rs*I;
P = V .* I;
[Pmpp, ind] = max(P);           %Returns the maximum power and where it occurs
Impp = I(ind);
Vmpp = V(ind);


%Simulink parameters for N panels
Load = (N*Voc)^2/(N*Pmin);
delta_Il = (N*Impp)*0.1/ 2;
L_inductor = (N*Vo) * t_PWM/ (4 * (N*Impp)*0.01);
k_e = 1/ L_inductor * 4;
C1 = t_PWM * (N*Impp)*0.1/ (8 * (N*Vmpp)*0.0001);
C2 = (N*Vo) * t_PWM/ (Load * (N*Vo)*0.001);
k_vc = 1/ C1;                   %v_pv error gain
k_l = 0.01/ C2;                 %v_o error gain (variable load)

%Parameters for Lyapunov integral term for Il calculation (reduce following
%error)
wn = 200;
csi = 1.25;
ki = wn^2;
kv = 2 * 1.25 * wn;


% %Parameters for integrated gama
% a1 = 3;
% wn_ig = 200;
% syms kI kv_ig ki_ig
% [SkI,Skv,Ski] = solve (1/(kI*(kv_ig+ki_ig))==1/(wn_ig^3),(2*kv_ig+ki_ig)/(kI*(kv_ig+ki_ig))==a1/(wn_ig^2), 1/(C1^2*kI*(kv_ig+ki_ig))+kv_ig/kI==a1/wn_ig, kI, kv_ig, ki_ig);
% kv_ig = max(double(Skv));
% ki_ig = max(double(Ski));
% ke_ig = ki_ig + kv_ig;


