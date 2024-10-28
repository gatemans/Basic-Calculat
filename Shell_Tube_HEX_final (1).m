clear all
close all
clc


%% Cold/Hot-Side Fluid Thermal Data:

T_c_in = input('Enter Cold-Side Inlet Temperature (C): ');
T_c_out = 1.3 .* T_c_in;
T_c_m = (T_c_in + T_c_out) / 2;

T_h_in = input('Enter Hot-Side Inlet Temperature (C): ');
T_h_out = input('Enter Hot-Side Outlet Temperature (C): ');
T_h_m = (T_h_in + T_h_out) / 2;


%% Water Properties: Density [kg/m3] (-30:150 C), Viscosity [Pa.s] (273:373 K), Thermal Conductivity [W/m.K] (274:370 K)

a = -2.8054253e-10; b = 1.0556302e-7; c = -4.6170461e-5; d = -0.0079870401;
e = 16.945176; f = 999.83952; g = 0.01687985;
rho_c_w = (f + e * T_c_m + d * T_c_m ^ 2 + c * T_c_m ^ 3 + b * T_c_m ^ 4 + a * T_c_m ^ 5) / (1 + g * T_c_m);
rho_h_w = (f + e * T_h_m + d * T_h_m ^ 2 + c * T_h_m ^ 3 + b * T_h_m ^ 4 + a * T_h_m ^ 5) / (1 + g * T_h_m);

A = -3.7188; B = 578.919; C = -137.546;
miu_c_w = (exp(A + (B / (C + T_c_m + 273.15)))) / 1000;
miu_h_w = (exp(A + (B / (C + T_h_m + 273.15)))) / 1000;

k_c_w = 0.6065 * (-1.48445 + 4.12292 * ((T_c_m + 273.15) / 298.15) - 1.63866 * ((T_c_m + 273.15) / 298.15) ^ 2);
k_h_w = 0.6065 * (-1.48445 + 4.12292 * ((T_h_m + 273.15) / 298.15) - 1.63866 * ((T_h_m + 273.15) / 298.15) ^ 2);

cp_w = 4193;                                % Specific Heat (Cold/Hot-Side) (J/kg.K)

alpha_c_w = k_c_w / (rho_c_w * cp_w);       % Thermal Diffusivity (Cold-Side)
alpha_h_w = k_h_w / (rho_h_w * cp_w);       % Thermal Diffusivity (Hot-Side)

Pr_c_w = miu_c_w / (rho_c_w * alpha_c_w);     % Prandtl Number (Cold-Side)
Pr_h_w = miu_h_w / (rho_h_w * alpha_h_w);     % Prandtl Number (Hot-Side)


%% Shell & Tube Thermal & Geometrical Data: (Assumptions: No fouling - Square arrangement - Material (Carbon steel))
k = 60;            % Shell & tubes thermal conductivity [W/m.k]
Q = input('Enter Heat Transfer Rate (W): ');    
dt = ((T_h_in - T_c_out) - (T_h_out - T_c_in)) / (log((T_h_in - T_c_out) / (T_h_out - T_c_in))) ;   % LMTD (C)
m_dot_c = Q / (cp_w * (T_c_out - T_c_in));                      % Cold flow mass flow rate (kg/s)
m_dot_h = Q / (cp_w * (T_h_in - T_h_out));                      % Hot flow mass flow rate (kg/s)

CL = 1;                                % Tube layout constant (90 degree)
CTP = 0.9;                             % Tube count calculation constant (Fixed tube sheet - Two passes)
PR = input('Enter Tube Pitch Ratio; '); %  (1.25 < PR < 1.5)

if PR < 1.25 || PR > 1.5
    disp('Error on PR value.')
else
end

OD = [0.25,0.25,0.25,0.375,0.375,0.375,0.375,0.5,0.5,0.5,0.5,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,...
    0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,1,1,1,1,1,1,...
    1,1,1,1,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,2,2,2.5] .* 25.4; % Tube outer diameter (mm)

t = 25.4 .* [0.028,0.022,0.018,0.049,0.035,0.028,0.022,0.065,0.049,0.035,0.028,0.109,0.095,0.083,0.072,0.065,...
    0.058,0.049,0.042,0.035,0.134,0.120,0.109,0.095,0.083,0.072,0.065,0.058,0.049,0.035,0.134,0.120,...
    0.109,0.095,0.083,0.065,0.049,0.035,0.165,0.134,0.120,0.109,0.095,0.083,0.072,0.065,0.049,0.035,...
    0.180,0.165,0.134,0.120,0.109,0.095,0.083,0.065,0.049,0.035,0.134,0.109,0.083,0.065,0.120,0.095,0.148]; % Tube thickness (mm)

ID = OD - 2 .* t;               % Tube inner diameter (mm)

PT = PR .* OD;                  % Tube pitch (mm)

A_f = (pi * ID .^ 2) ./ 4;      % Internal flow area (mm2)

Ds = 500;                       % Shell diameter (mm)

L = [0.5:0.5:10];

B = 200;        % Baffle Spacing (mm)

C = PT - OD;    % Clearance 

for i = 1 : 65
    Nt(i) = ceil(0.785 * (CTP / CL) * ((Ds ^ 2) / (PR ^ 2 * OD(i) ^ 2)));       % No. of tubes
end

De = 4 .* (PT .^ 2 - (pi .* OD() .^ 2) ./ 4) ./ (pi .* OD());             % Equivalent shell diameter (mm)

As = (Ds * B .* C) ./ PT;

for j = 1:20

A_o = pi .* OD .* L(j) .* Nt;

A_i = pi .* ID .* L(j) .* Nt;

end

Gs = (m_dot_h ./ As) .* 10e6;         % Shell-side mass velocity [kg/s.m2]

Re_s = Gs .* De ./ (1000 * miu_h_w);

h_s = (360 .* k_h_w .* ((De .* Gs ./ (1000 .* miu_h_w)) .^ 0.55) .* ((cp_w * miu_h_w / k_h_w) ^ 0.33)) ./ De; % Shell-side Convection (W/m2.K)

A_tp = ((pi .* ID .^ 2) .* Nt) ./ (8 * 10e6);   % Total tube inside area (m2)

u_t = m_dot_c ./ (rho_c_w .* A_tp);     % Tube velocity (m/s)

Re_t = rho_c_w .* u_t .* ID ./ (1000 * miu_c_w);

f = (1.58 .* log(Re_t) - 3.28) .^ (-2);

Nu_t = (f .* (Re_t - 1000) .* Pr_c_w) ./ (2 + 25.4 .* (f ./ 2) .^ 0.5 .* (Pr_c_w .^ 0.667 - 1));

h_t = (Nu_t .* ID) ./ k_c_w;

[x,y] = max(abs((1 ./ ((1 ./ h_s) + (1 ./ h_t)))));

fprintf('Overall Heat Transfer Coefficient = %6.4f \n',x)

fprintf('Thickness = %4.2f \n',t(y))

fprintf('Inner Diamter = %4.2f \n',ID(y))

fprintf('Outer Diamter = %4.2f \n',OD(y))

fprintf('No. of Tubes = %4.2f \n',Nt(y))