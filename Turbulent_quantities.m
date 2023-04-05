clc
clear
close all

%Data
L = 2;                                                      % Bump chord [m]
T_ref = 300;                                                % Total temperature [K]
rho_ref = 1.174;                                            % Total pressure [Pa]
gamma = 1.4;                                                % Heat capacity ratio
Ma = 0.2;                                                   % Mach number
R = 287;                                                    % Air constant [J/Kg*K]
P_ref = rho_ref*(R*T_ref);                                  % Total density [Kg/m^3]

eta = 1+((Ma^2)*(gamma-1)/2);

%T_t = T_ref*1.008;                                          % Static temperature [K]
%P_t = P_ref*1.02828;                                        % Static pressure [Pa]
%rho_t =  P_t/(R*T_t);                                       % Static denisty [Kg/m^3]

T_t = T_ref*eta;                                            % Static temperature [K]
P_t = P_ref*(eta^(gamma/(gamma-1)));                        % Static pressure [Pa]
rho_t =  rho_ref*(eta^(1/gamma-1));                         % Static denisty [Kg/m^3]

%From NASA - Sandia experiment
ratio_p = 0.5825;                                           % Exit wall pressure ratio P_exit/P_t
%P_exit = P_t*ratio_p;                                       % Exit pressure [Pa]

c = sqrt(gamma*R*T_ref);                                    % Speed of sound [m/s]

U = c*Ma;                                                   % Free stream velocity [m/s]
u_square = U^2;

Re =5e+6;                                                   % Reynolds number

nu = (U*L)/Re;                                              % Kinematic viscosity [m/s^2]                          
mu = rho_ref*nu;                                            % Dynamic viscosisty [Kg/m*s] Obtained from mu and rho_ref

%Sutherland lsqrt(T_ref)/(1+(Ts/T_ref))aw
Ts = 110.4;                                                 % Sutherland temperature [K]
As = 1.458e-06;                                             % Sutherland constant [kg/m*s*K^0.5]

mu_s = As * sqrt(T_ref)./(1+(Ts/T_ref));                    % Sutherland Dynamic viscosisty [Kg/m*s]
nu_s = mu_s/rho_ref;                                        % Sutherland Kinematic viscosity [m^2/s]
nut_s = 3*nu_s;                                             % Sutherland turbulent kinematic viscosity [m^2/s]
Re_s = (U*L)/nu_s;                                          % Sutherland Reynolds number 


%Turbulent data for Spalart-Allmaras method
Prt = 0.85;                                                 % Turbulent Prandtl number
nut = 3*nu;                                                 % Turbulent kinematic viscosity [m^2/s]
nuTilda = nut;                                              % Spalart-Allmaras kinematic viscosity [m^2/s]
alphat = nut/Prt;                                           % Kinematic turbulent thermal conductivity [m^2/s] 

%Turbulent data for kOmega method
I = 1/100;                                                  % Turbulence intensity
C_mu = 0.09;                                                % Constant

k = 1.5*(I*U)^2;                                            % Turbulent kinetic energy [m^2/s^2]
omega = (sqrt(k))/((C_mu^0.25)*L);                          % Turbulence specific dissipation rate [1/s]

%Wall quantities
yPlus = 6;                                                % Adimensional y at inner layer
Cf = 0.026/(Re)^(1/7);                                      % Friction coefficient at the wall
tau_w = 0.5*Cf*rho_ref*U^2;                                 % Wall shear stress [Pa]
U_tau = sqrt(tau_w/rho_ref);                                % Friction velocity [m/s]
Delta_s = (yPlus*nu)/U_tau;                                 % First cell height at the wall to obtain yPlus=1
