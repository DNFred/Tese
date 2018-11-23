clear all;
close all;
format compact;
clc;
pause('on');

%Global constants
K = 1.38 * 10^(-23);            %Boltzmann constant
T = 273.15;                     %0 degrees Celsius in absolute scale 
q = 1.6 * 10^(-19);             %electron charge

%Catalog values at STC (Standard Test 400^2Conditions)
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
Npanels = 3;    %<------------        %Number of panels
Vo = 400;                       %Load voltage
t_PWM = 50e-6;                  %PWM signal period
k_v = 1/ (t_PWM*500);%10000);           %voltage gain for v_ref
k_vc = 1/ (t_PWM*50);%);

t_int = 0.5;
trigger_signal = 1e-3;
GREEN = 0;
ORANGE = 1;
RED = 2;
ON = 1;
OFF = 0;
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
Cp = 3.5e-6;
Ior = (Iscr - (Vocr - Rs*Iscr)/ Rsh) * exp(-Vocr/ (m*Vtr));
Isr = Ior * exp(Vocr/ (m*Vtr)) + Vocr/Rsh;


% figure       %only active if doing IV and PV curve demo
% while G>200


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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %IV-curve
% plot(V,I);
% xlim([0 40]);
% ylim([0 15]);
% xlabel('Voltage');
% ylabel('Current');
% hold on
% 
% % %PV-curve
% % plot(V,P);
% % xlim([0 40]);
% % ylim([0 300]);
% % xlabel('Voltage');
% % ylabel('Power');
% % hold on
% G = G-200;
% end
% 
% axp = gca;
% axp.XAxisLocation = 'origin';
% axp.YAxisLocation = 'origin';
% % determine startpoint and endpoint for the arrows 
% xs=0.1305;
% xe=0.95;
% ys=0.1098;
% ye=0.95;  
% % make the arrows
% annotation('arrow', [xs xe],[ys ys]);
% annotation('arrow', [xs xs],[ys ye]);
% % remove old box and axes
% box off
% set(gca,'YTick',[])
% set(gca,'yticklabel',({}))
% set(gca,'XTick',[])
% set(gca,'xticklabel',({}))
% legend('1200 W/m^2', '1000 W/m^2', '800 W/m^2', '600 W/m^2', '400 W/m^2','Location', 'northeast');
% set(gca, 'FontSize', 16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %MPPT problem
% figure
% plot(simout)
% hold on
% plot(20*simout1)
% ylim([0 400]);
% xlim([0 1.03]);
% xlabel('Voltage');
% axp = gca;
% axp.XAxisLocation = 'origin';
% axp.YAxisLocation = 'origin';
% xs=0.1305;
% xt=0.9;
% xe=0.95;
% ys=0.1898;
% ye=0.95;  
% annotation('arrow', [xs xe],[ys ys]);
% annotation('arrow', [xs xs],[ys ye]);
% %annotation('arrow', [xt xt],[ys ye]);
% box off
% set(gca,'YTick',[])
% set(gca,'yticklabel',({}))
% set(gca,'XTick',[Vmpp/(3*Voc), (Voc*0.9+Vmpp)/(3*Voc), (1.92*Voc+Vmpp)/(3*Voc)])
% set(gca,'xticklabel',({'V_1','V_2','V_3'}))
% lgd = legend('$P_{PV}$', '$I_{PV}$', 'Location', 'northwest');
% set(lgd,'Interpreter','latex');
% set(lgd,'FontSize',17);
% set(gca,'FontSize',17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Demonstration of the IV curve of solar cell
% Id = Io * (exp(Vd_vetor/ (m*Vt)) - 1);
% figure
% plot(V,Id,'k');
% hold on
% Ip = Id - Is;
% Vp = plot(V,Ip,'k');
% set(Vp,'LineWidth',2)
% xlim([0 40]);
% ylim([-10 10]);
% lgd = legend('$I_{dark}$', '$I_{cell}$','Location','northwest');
% set(lgd,'Interpreter','latex');
% set(gca, 'FontSize', 28)
% axp = gca;
% axp.XAxisLocation = 'origin';
% axp.YAxisLocation = 'origin';
% % determine startpoint and endpoint for the arrows 
% xs=0.1305;
% xe=0.95;
% ys=0.05;
% ym=0.516;
% ye=0.95;    
% % make the arrows
% annotation('arrow', [xs xe],[ym ym]);
% annotation('arrow', [xs xs],[ye ys]); 
% box off
% Ipe = -Isc;
% set(gca,'YTick',[Ipe,0])
% set(gca,'yticklabel',({'I_{pe}','0'}))
% set(gca,'XTick',[0])
% set(gca,'xticklabel',({'0'}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PV-curve for demonstration (Tc = ta + G*(NOCT-20)/ 800 + 273.15;)
% figure
% VP = plot(V,P*0.7,'k');
% set(VP,'LineWidth',2)
% hold on
% VI = plot(V(1:71),I(1:71)*10,'k');
% xlim([0 Voc*1.1])
% xlabel('Voltage')
% ylim([-50 250*0.7])
% ylabel('Current')
% dP = diff(P)./diff(V);
% plot(V(2:end)-0.52,dP*6,'--k');
% lgd = legend('$P_{pan}$', '$I_{pan}$', '$\frac{dP_{pan}}{dV_{pan}}$','Location','northwest');
% set(lgd,'Interpreter','latex');
% set(lgd,'FontSize',17);
% VM = plot([0,Vmpp],[Impp*10,Impp*10],':k');
% IM = plot([Vmpp,Vmpp],[0,Pmpp*0.7],':k');
% PM = plot([Vmpp, Voc*1.1],[Pmpp*0.7,Pmpp*0.7],':k');% 
% %%%%%%%%%%Put arrows in the Axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % determine position of the axes
% axp = gca;
% axp.XAxisLocation = 'origin';
% axp.YAxisLocation = 'origin';
% % determine startpoint and endpoint for the arrows 
% xs=0.1305;
% xt=0.9;
% xe=0.95;
% ys=0.291;
% ye=0.95;  
% annotation('arrow', [xs xe],[ys ys]);
% annotation('arrow', [xs xs],[ys ye]);
% annotation('arrow', [xt xt],[ys ye]);
% % remove old box and axes
% box off
% set(gca,'YTick',[Impp*10, Pmpp*0.7])
% set(gca,'yticklabel',({'I_{MPP}','P_{MPP}'}))
% set(gca,'XTick',[Vmpp])
% set(gca,'xticklabel',({'V_{MPP}'}))
% % set(gca,'YTick',[Impp*10, Isc*10, Pmpp*0.7])
% % set(gca,'yticklabel',({'I_{MPP}','I_{SC}','P_{MPP}'}))
% % set(gca,'XTick',[Vmpp Voc])
% % set(gca,'xticklabel',({'V_{MPP}','V_{OC}'}))
% set(gca, 'FontSize', 17)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PV-curve for voltage control demostrantion (Tc = ta + G*(NOCT-20)/ 800 + 273.15;)
figure
VP = plot(V,P,'k','LineWidth',2);
% set(VP,'LineWidth',2)
hold on
% VI = plot(V(1:71),I(1:71)*20,'k');
% lgd = legend('$P_{PV}$','$I_{PV}$','Location','northwest');
% set(lgd,'Interpreter','latex');
%set(lgd,'FontSize',30);
xlim([0 Voc*1.1])
ylim([0 Pmpp*1.1])
VM = plot([Vmpp,Vmpp],[0,Pmpp*1.1],':k');
left = '$$\rightarrow\frac{dP_{PV}}{dV_{PV}} < 0$$';
text(Vmpp*1.05,Pmpp*1.005,left,'Interpreter','latex','FontSize',30)
right = '$$\frac{dP_{PV}}{dV_{PV}} > 0 \leftarrow$$';
text(Vmpp*0.737,Pmpp*1.005,right,'Interpreter','latex','FontSize',30)
xlabel('Voltage');
set(gca,'XTick',[Vmpp])
set(gca,'xticklabel',({'V_{MPP}'}))
set(gca,'YTick',[])
set(gca,'yticklabel',({''}))
set(gca, 'FontSize', 28)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Hysteresis loop for converter control
% x = [-2 1 1 2 -1 -1 -2];
% y = [1 1 0 0 0 1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% q.ShowArrowHead = 'off';
% hold on;
% 
% headWidth = 5;
% headLength = 5;
% x = [0 0.5];
% y = [1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [1 1];
% y = [1 0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [2 1.5];
% y = [0 0];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [0 -0.5];
% y = [0 0];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-1 -1];
% y = [0 0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-2 -1.5];
% y = [1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% set(gca,'YTick',[0, 1])
% set(gca,'yticklabel',({'0','1'}))
% set(gca,'XTick',[-1 0 1])
% set(gca,'xticklabel',({'-\DeltaI_L/2','0','\DeltaI_L/2'}))
% xlabel('e_I_L');
% ylabel('\gamma');
% set(gca, 'FontSize', 17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Bode plot
% H = tf([(2*pi*100)^2 0],[1 2*0.7*2*pi*100 (2*pi*100)^2]);
% bode(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Hysteresis loop for inverter control
% x = [2 -0.5 -0.5 1 -1 -1 -2 0.5 0.5 1 1 2];
% y = [1 1 0 0 0 -1 -1 -1 0 0 1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% q.ShowArrowHead = 'off';
% hold on;
% 
% headWidth = 5;
% headLength = 5;
% x = [0.5 0];
% y = [1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-0.5 -0.5];
% y = [1 0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-1 -1];
% y = [0 -0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-2 -1.5];
% y = [-1 -1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [-0.5 0];
% y = [-1 -1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [0.5 0.5];
% y = [-1 -0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [1 1];
% y = [0 0.5];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% x = [2 1.5];
% y = [1 1];
% q = quiver(x(1:end-1),y(1:end-1),diff(x),diff(y),0,'color',[0 0 0]);
% U = q.UData;
% V = q.VData;
% X = q.XData;
% Y = q.YData;
% ah = annotation('arrow','headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
% set(ah,'parent',gca);
% set(ah,'position',[X Y U V]);
%         
% set(gca,'YTick',[-1, 0, 1])
% set(gca,'yticklabel',({'-1','0','1'}))
% set(gca,'XTick',[-1 -0.5 0 0.5 1])
% set(gca,'xticklabel',({'-\DeltaI_L/2','-\DeltaI_L/4','0','\DeltaI_L/4','\DeltaI_L/2'}))
% xlabel('e_{inv}');
% ylabel('\gamma_i');
% set(gca, 'FontSize', 17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Input Capacitor C_1 current
% figure
% %Capacitor current
% x=[0 2 3 7 8];
% y = [0 1 -1 1 -1];
% plot(x,y)
% hold on
% %Switch position
% x=[0 2 2 3 3 7 7 8];
% y=[1.5 1.5 0 0 1.5 1.5 0 0];
% plot(x,y)
% x=[0 2 2.5];
% y=[0 1 0];
% h=area(x,y);
% h.FaceColor = [0 1 0];
% ylim([-2 2]);
% ylabel('I_C');
% xlim([0 8]);
% axp = gca;
% axp.XAxisLocation = 'origin';
% axp.YAxisLocation = 'origin';
% xs=0.1305;
% xe=0.95;
% ys=0.1098;
% yt=0.5165;
% ye=0.95;  
% annotation('arrow', [xs xe],[yt yt]);
% annotation('arrow', [xs xs],[ys ye]);
% box off
% set(gca,'YTick',[-1 0 1])
% set(gca,'yticklabel',({'-\Deltai_C/2', '0', '\Deltai_C/2'}))
% set(gca,'XTick',[2.5 5])
% set(gca,'xticklabel',({'T/2', 'T'}))
% ylim = get(gca,'XLim');
% h = line([2 2],[-max(ylim)*2/3 max(ylim)*2/3]);
% set(h,'LineStyle',':')
% h = line([3 3],[-max(ylim)*2/3 max(ylim)*2/3]);
% set(h,'LineStyle',':')
% h = line([7 7],[-max(ylim)*2/3 max(ylim)*2/3]);
% set(h,'LineStyle',':')
% on = 'S ON';
% text(0.7,1.2,on)
% text(4.7,1.2,on)
% off = 'S OFF';
% text(2.15,1.2,off)
% text(7.15,1.2,off)
% q = '$$\Delta Q$$';
% text(1.5,0.4,q,'Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Boost parameters for N panels
Np = 1;     %<------------      %1 for all
delta_Il = (Np*Impp)*0.1/ 2;
L_inductor = Vo*t_PWM/ (4 * 2*delta_Il);
Ns = 1;     %<------------      %1 for all
C1 = t_PWM * (Np*Impp)*0.1/ (8 * (Ns*Vmpp)*0.0001);
Nsc = 3;    %<------------      %3 for all
RLoad = (Vo/Nsc)^2/Pmpp;
C2 = Vo/Nsc * t_PWM/ (RLoad * Vo/Nsc*0.001);%N*Impp/(2*pi*50*0.01*Vo);

%Inverter parameters
delta_Ii = 0.05*Npanels*Impp;
Npc = 1;    %<------------      %1 for 1MPPT and (series)/ 3 for (parallel)
Li = Vo*t_PWM/ (4*delta_Ii);
Cdc_inv = Npc*Impp/(2*pi*50*0.01*Vo);
k_inv = 0.005/Cdc_inv;
k_A = 1/Cdc_inv;
Vgrid = 230;
GAMA = Vgrid/400;

syms kIi kvi kil wn_i ai Ci
[SkIi,Skvi,Skii] = solve (1/(kIi*kil)==1/(wn_i^3),(kvi+kil)/(kIi*kil)==ai/(wn_i^2), GAMA^2/(Cdc_inv^2*kIi*kil)+(kIi+kvi)/(kIi*kil)==ai/wn_i, kIi, kvi, kil);
 
tdi = 100*t_PWM;
ai = 3;
wn_i = 1/(ai*tdi);
 
kIc1i = simplify (subs (SkIi), 'Criterion','preferReal', 'Steps', 100);
kvc1i = simplify (subs (Skvi), 'Criterion','preferReal', 'Steps', 100);
kic1i = simplify (subs (Skii), 'Criterion','preferReal', 'Steps', 100);
 
kIc2i = eval (kIc1i);
kvc2i = eval (kvc1i);
kic2i = eval (kic1i);

k_inv = kvc2i(2)/10;%*10
kIi = kIc2i(2);
k_A = kIi/10;%*10

Td = 1/50;
a = 2.7;    
Gi = -GAMA;
Tzv = a^2*Td;
Tpv = a^3*Gi*Td^2/Cdc_inv;

%Parameters for voltage control
wn = 2;
csi = 1.25;
k1 = wn^2;
k2 = 2 * 1.25 * wn;


%Parameters for current control
%(versão professor corrigida)
syms kI kv ki wn_ig a1 Cq
[SkI,Skv,Ski] = solve (1/(kI*ki)==1/(wn_ig^3),(kv+ki)/(kI*ki)==a1/(wn_ig^2), 1/(C1^2*kI*ki)+(kI+kv)/(kI*ki)==a1/wn_ig, kI, kv, ki);
 
td = 10*t_PWM;
a1 = 3;
wn_ig = 1/(a1*td);
 
kIc1 = simplify (subs (SkI), 'Criterion','preferReal', 'Steps', 100);
kvc1 = simplify (subs (Skv), 'Criterion','preferReal', 'Steps', 100);
kic1 = simplify (subs (Ski), 'Criterion','preferReal', 'Steps', 100);
 
kIc2 = eval (kIc1);
kvc2 = eval (kvc1);
kic2 = eval (kic1);

kv_ig = kvc2(2);
kI = kIc2(2);
ki_ig = kI;

%Iluminance tests
t1 = 6;
t2 = 11;
t3 = 16;
t4 = 8;
t5 = 22;
t6 = 20;

g1 = 1000;
g2 = 1000/2.25;
g3 = 1000/3.5;

%Tests - Results
% figure
% plot(Po)
% xlim([4 16]);
% ylabel('Power [W]');