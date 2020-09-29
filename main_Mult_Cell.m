% Simulink simulator for Interval Observer for cells in series

clear;clc;close all;
n = 5; % Number of cells in series

%% Set model parameters

load('data\r_rc_rc_scale_cap_A123_26650.mat')
% With temperature effect - Extract all values
Rs2_SOC_c = Rs2_SOC_c*10;
Rs2_SOC_d = Rs2_SOC_d*10;

% Cc = C12_SOC_c(:,6:7);
% R1c = Rs2_SOC_c(:,6:7);
% R2c = R12_SOC_c(:,6:7);
% Cd = C12_SOC_d(:,6:7);
% R1d = Rs2_SOC_d(:,6:7);
% R2d = R12_SOC_d(:,6:7);

C_c = C12_SOC_c(:,6);
R1_c = Rs2_SOC_c(:,6);
R2_c = R12_SOC_c(:,6);
C_d = C12_SOC_d(:,6);
R1_d = Rs2_SOC_d(:,6);
R2_d = R12_SOC_d(:,6);

tauc = 1./R2_c./C_c;
taud = 1./R2_d./C_d;
C_max = max(max(max(C_c)),max(max(C_d)));
C_min = min(min(min(C_c)),min(min(C_d)));
% C0 = mean(mean(C));
R2_max = max(max(max(R2_c)),max(max(R2_d)));
R2_min = min(min(min(R2_c)),min(min(R2_d)));
% R20 = mean(mean(R2));
tau0 = mean([mean(mean(tauc)),mean(mean(taud))]);
delta_tauc = tauc - tau0;
delta_taud = taud - tau0;
delta_tau_max = max(max(max(delta_tauc)),max(max(delta_taud)));
delta_tau_min = min(min(min(delta_taud)),min(min(delta_taud)));
R1_max = max(max(max(R1_c)),max(max(R1_d)));
R1_min = min(min(min(R1_c)),min(min(R1_d)));

% Thermal model parameters
Vdot = 9.5e-3;
Rho = 1.205;
Cp = 1005;
Cf = Vdot*Rho*Cp;
Cc = 62.7;
Cs = 4.5;
Rc = 1.94;
Ru = 15;
Rcc = 1000/1000;

a = -1/Rc/Cc;
b = 1/Rc/Cc;
c = 1/Rc/Cs;
d = -(1/Ru/Cs+1/Rc/Cs);
e = 1/Ru/Cs;

% Electrical model parameters
eta = 1; % Coulombic efficiency
gamma = 1; % decay rate
M = 0.05; % maximum polarization due to hysteresis

%% Load input profile
% UDDS
load('data\UDDSx2_DFN_Co_1sec')
t_UDDS = out.time + 2;
I_UDDS = out.cur/15; % Scale down magnitude of input current
Tend = 1000;
index = find(t_UDDS == Tend);
I_UDDS = I_UDDS(1:index);
t_UDDS = t_UDDS(1:index);
delta_t = 0.1; % Reduce the step size from 1 to delta_t
t = (0:delta_t:Tend)';
Cur = interp1(t_UDDS,I_UDDS,t); % Interpolate
Tf1 = 25*ones(size(t)); % Coolant flow temperature of cell #1

% % Constant current
% delta_t = 0.01;
% t = (0:delta_t:99)';
% Cur = -5*ones(size(t));
% Tf1 = 25*ones(size(t)); % Coolant flow temperature of cell #1
% Tend = t(end);

%% Set initial conditions for dynamical states
% Initial conditions for SOC, Vc, h, Tc, Ts
SOC0 = [0.28; 0.3; 0.32; 0.34; 0.36]; %0.35*ones(n,1); %
Vc0 = zeros(n,1);
h0 = zeros(n,1);
Tc0 = 25*ones(n,1);
Ts0 = 25*ones(n,1);

%% Determine consistent initial conditions for algebraic states
H = 0; % Indicate if hysteresis is included
if(H == 1)
    load('data\LFP_OCV_Data.mat','SOC_LFP') % OCV for LFP
    SOC_data = SOC_LFP;
    OCV_data = OCV_LFP_func(SOC_data);
    clear SOC_LFP
    [R1_init,R2_init,C_init,I_init] = Alg_init_hyst(Ts0,SOC0,Vc0,h0,Cur);
    Q = 3600*[2.3 2.4 2.2 2.15 2.35]; % [Coulombs] LFP = 2.3Ah
    Q_max = max(Q); Q_min = min(Q);
    OCV_min = OCV_LFP_func(0.1); OCV_max = OCV_LFP_func(0.9); 
else
    load('data\NMC_OCV_Data.mat') % OCV for LFP
    SOC_data = SOC_NMC;
    OCV_data = OCV_NMC_func(SOC_data);
    clear SOC_NMC
    [R1_init,R2_init,C_init,I_init] = Alg_init(Ts0,SOC0,Vc0,Cur);
    Q = 3600*[2.8 2.9 2.7 2.65 2.85]; % [Coulombs] NMC = 2.8Ah
    Q_max = max(Q); Q_min = min(Q);
    OCV_min = OCV_NMC_func(0.1); OCV_max = OCV_NMC_func(0.9); 
end

%% Initial conditions for interval observer

SOC_upper0 = max(SOC0) + 0.05;
SOC_lower0 = min(SOC0) - 0.05;
Vc_upper0 = max(Vc0) + 0.01;
Vc_lower0 = min(Vc0) - 0.008;
Ts_lower0 = min(Ts0) - 2;
Ts_upper0 = max(Ts0) + 2;
Tc_lower0 = min(Ts0) - 3;
Tc_upper0 = max(Ts0) + 3;

z_1_upper0 = OCV_NMC_func(SOC_upper0) + Vc_lower0 - OCV_NMC_func(0);
z_1_lower0 = OCV_NMC_func(SOC_lower0) + Vc_upper0 - OCV_NMC_func(0);
z_2_upper0 = -tau0*Vc_lower0;
z_2_lower0 = -tau0*Vc_upper0;

L1_lower = 10;
L1_upper = 10;
L2_lower = -0.00002;
L2_upper = -0.00002;
L1t_lower = b/1.01;
L1t_upper = b/1.01;
L2t_lower = 1;
L2t_upper = 1;

load('OCV_Map')

%% Simulate and parse out states

sim('Mult_Cell.slx')

% Parse out states
SOC_smlk = nan(length(tout),n);
Vc_smlk = nan(length(tout),n);
V_smlk = nan(length(tout),n);
Tc_smlk = nan(length(tout),n);
Ts_smlk = nan(length(tout),n);
R1_smlk = nan(length(tout),n);
R2_smlk = nan(length(tout),n);
C_smlk = nan(length(tout),n);
Tf_smlk = nan(length(tout),n); Tf_smlk(:,1) = Tf1;

for i = 1:(7*n-1)
    if(i <= 5)
        SOC_smlk(:,i) = yout{i}.Values.Data;
    elseif(i <= 10)
        Vc_smlk(:,i-n) = yout{i}.Values.Data;
    elseif(i <= 15)
        V_smlk(:,i-2*n) = yout{i}.Values.Data;
    elseif(i <= 20)
        Tc_smlk(:,i-3*n) = yout{i}.Values.Data;
    elseif(i <= 25)
        Ts_smlk(:,i-4*n) = yout{i}.Values.Data;
    elseif(i <= 30)
        R1_smlk(:,i-5*n) = yout{i}.Values.Data(:,1);
        R2_smlk(:,i-5*n) = yout{i}.Values.Data(:,2);
        C_smlk(:,i-5*n) = yout{i}.Values.Data(:,3);
    else
        Tf_smlk(:,i-6*n+1) = yout{i}.Values.Data;
    end
end

SOC_upper = yout{35}.Values.Data;
SOC_lower = yout{36}.Values.Data;
Vc_upper = yout{37}.Values.Data;
Vc_lower = yout{38}.Values.Data;

z1_lower = yout{39}.Values.Data;
z2_lower = yout{40}.Values.Data;
z1_upper = yout{41}.Values.Data;
z2_upper = yout{42}.Values.Data;

Tc_lower = yout{43}.Values.Data;
Ts_lower = yout{44}.Values.Data;
Tc_upper = yout{45}.Values.Data;
Ts_upper = yout{46}.Values.Data;

%%
close all;
figure()
plot(tout,z1a.Data,tout,z1b.Data,tout,z1c.Data,tout,z1d.Data,tout,z1e.Data,'LineWidth',1.5)
hold on
plot(tout,z1_lower,'k--',tout,z1_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$z_1(t)$','Interpreter','latex')
grid on
set(gca,'FontSize',14)
xlim([0 t(end)])

figure()
plot(tout,z2a.Data,tout,z2b.Data,tout,z2c.Data,tout,z2d.Data,tout,z2e.Data,'LineWidth',1.5)
hold on
plot(tout,z2_lower,'k--',tout,z2_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$z_2(t)$','Interpreter','latex')
grid on
set(gca,'FontSize',14)
xlim([0 t(end)])

figure()
plot(tout,SOC_smlk(:,1),tout,SOC_smlk(:,2),tout,SOC_smlk(:,3),tout,SOC_smlk(:,4),tout,SOC_smlk(:,5),'LineWidth',1.5)
hold on
plot(tout,SOC_lower,'k--',tout,SOC_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$SOC(t)$ [A]','Interpreter','latex')
grid on
set(gca,'FontSize',14)
leg = legend('$z_1$','$z_2$','$z_3$','$z_4$','$z_5$','$\underline{z}$','$\overline{z}$');%,);
set(leg,'Interpreter','latex','FontSize',16);
xlim([0 t(end)])

figure()
plot(tout,Vc_smlk(:,1),tout,Vc_smlk(:,2),tout,Vc_smlk(:,3),tout,Vc_smlk(:,4),tout,Vc_smlk(:,5),'LineWidth',1.5)
hold on
plot(tout,Vc_lower,'k--',tout,Vc_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$V_c(t)$ [V]','Interpreter','latex')
grid on
set(gca,'FontSize',14)
xlim([0 t(end)])

figure()
plot(tout,Tc_smlk(:,1),tout,Tc_smlk(:,2),tout,Tc_smlk(:,3),tout,Tc_smlk(:,4),tout,Tc_smlk(:,5),'LineWidth',1.5)
hold on
plot(tout,Tc_lower,'k--',tout,Tc_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$T_c(t)$ [$^o$C]','Interpreter','latex')
grid on
set(gca,'FontSize',14)
xlim([0 t(end)])

figure()
plot(tout,Ts_smlk(:,1),tout,Ts_smlk(:,2),tout,Ts_smlk(:,3),tout,Ts_smlk(:,4),tout,Ts_smlk(:,5),'LineWidth',1.5)
hold on
plot(tout,Ts_lower,'k--',tout,Ts_upper,'r--','LineWidth',1.5)
xlabel('$t$ [Seconds]','Interpreter','latex')
ylabel('$T_s(t)$ [$^o$C]','Interpreter','latex')
grid on
set(gca,'FontSize',14)
xlim([0 t(end)])
