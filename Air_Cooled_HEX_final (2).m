clear all
close all
clc


%% Cold/Hot-Side Fluid Thermal Data:

% Air Data Input
T_a_in = input('Enter Air Inlet Temperature (C): ');
T_a_out = 0.7 .* T_a_in;
T_a_m = (T_a_in + T_a_out) / 2;
u = input('Enter Air Velocity (m/s): ');

% Water Data Input
T_w_in = input('Enter Water Inlet Temperature (C): ');
T_w_out = input('Enter Water Outlet Temperature (C): ');
T_w_m = (T_w_in + T_w_out) / 2;


%% Water Properties:

% Density [kg/m3] (-30:150 C)
a = -2.8054253e-10; b = 1.0556302e-7; c = -4.6170461e-5; d = -0.0079870401;
e = 16.945176; f = 999.83952; g = 0.01687985;
rho_w = (f + e * T_w_m + d * T_w_m ^ 2 + c * T_w_m ^ 3 + b * T_w_m ^ 4 + a * T_w_m ^ 5) / (1 + g * T_w_m);

% Viscosity [Pa.s] (273:373 K)
A = -3.7188; B = 578.919; C = -137.546;
miu_w = (exp(A + (B / (C + T_w_m + 273.15)))) / 1000;

% Thermal Conductivity [W/m.K] (274:370 K)
k_w = 0.6065 * (-1.48445 + 4.12292 * ((T_w_m + 273.15) / 298.15) - 1.63866 * ((T_w_m + 273.15) / 298.15) ^ 2);

% Specific Heat [J/kg.K]
cp_w = 4193;

% Thermal Diffusivity [m2/s]
alpha_w = k_w / (rho_w * cp_w);

% Prandtl Number
Pr_w = miu_w / (rho_w * alpha_w);


%% Air Properties: (Ideal Gas Assumption - Atmospheric Pressure)

% Density [kg/m3]
rho_a = 101325 / (287.058 * (T_a_m + 273.15));

% Viscosity [Pa.s]
miu_a = (1.458 * (10 ^ -6) * (T_a_m + 273.15) ^ 1.5) / (110.4 + 273.15 + T_a_m);

% Thermal Conductivity [W/m.K]
k_a = (1.5207 * (10 ^ -11) * (T_a_m + 273.15) ^ 3) - (4.8574 * (10 ^ -8) * (T_a_m + 273.15) ^ 2) + (1.0184 * (10 ^ -4) * (T_a_m + 273.15)) - (3.9333 * (10 ^ -4));

% Specific Heat [J/kg.K]
cp_a = 1005;

% Thermal Diffusivity [m2/s]
alpha_a = k_a / (rho_a * cp_a);

% Prandtl Number
Pr_w = miu_a / (rho_a * alpha_a);


%% Heat Exchanger Geometrical Data: (Finned-Tube)

T_lmtd = ((T_w_in - T_a_in) - (T_w_out - T_a_out)) / log(((T_w_in - T_a_in)) / (T_w_out - T_a_out));

Q = input('Enter Heat Transfer Rate (W): ');

m_dot_a = Q / (cp_a * (T_a_out - T_a_in));      % Air mass flow rate (kg/s)

m_dot_w = Q / (cp_w * (T_w_in - T_w_out));      % Water mass flow rate (kg/s)

C_a = m_dot_a * cp_a;

C_w = m_dot_w * cp_w;

C_min = min(C_a,C_w);

C_max = max(C_a,C_w);

Q_max = C_max * (T_w_in * T_a_in);

e = Q / Q_max;

C_ratio = C_min / C_max;

if C_a < C_w
    NTU = (-log(C_ratio * (log(1 - e)) + 1)) / C_ratio;
else
    NTU = -log(1 + (log(1 + C_ratio * e)) / C_ratio);
end

NR = [3, 4, 5, 6];              % No. of rows

sigma = [315, 394, 394];        % Fin density (fins/m)

P = [0.0603, 0.0603, 0.0635];   % Pitch (m)

FV = [3.3, 3.12, 2.97, 2.84, 3.18, 3.05, 2.92, 2.79, 3.56, 3.35, 3.18, 3.05];            % Face velocity (m/s)

OD = [0.25,0.25,0.25,0.375,0.375,0.375,0.375,0.5,0.5,0.5,0.5,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,...
    0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,1,1,1,1,1,1,...
    1,1,1,1,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,2,2,2.5] .* 25.4 ./ 1000; % Tube outer diameter (m)

t = 25.4 .* [0.028,0.022,0.018,0.049,0.035,0.028,0.022,0.065,0.049,0.035,0.028,0.109,0.095,0.083,0.072,0.065,...
    0.058,0.049,0.042,0.035,0.134,0.120,0.109,0.095,0.083,0.072,0.065,0.058,0.049,0.035,0.134,0.120,...
    0.109,0.095,0.083,0.065,0.049,0.035,0.165,0.134,0.120,0.109,0.095,0.083,0.072,0.065,0.049,0.035,...
    0.180,0.165,0.134,0.120,0.109,0.095,0.083,0.065,0.049,0.035,0.134,0.109,0.083,0.065,0.120,0.095,0.148] ./ 1000; % Tube thickness (m)

ID = OD - 2 .* t;               % Tube inner diameter (m)

for i = 1:12

if i < 4

h(i) = (6.75 .* (FV(i) .* 196.85) .^ 0.5) .* 5.678263341;

else

h(i-3) = (8 .* (FV(i) .* 196.85) .^ 0.5) .* 5.678263341;

end

end

FA = ((Q .* 3.4144) ./ (FV .* (T_w_in - T_w_out) .* 1.95 .* 196.85)) .* 0.0929;   % Face area (m2)

L = 1:1:14;         % Length (m)

Y = 0;

for j = 1:14

for z = 1:12

    Y(j,z) = (FA(z) ./ L(j)).*1000;         % Width (m)

end

end

Nt_1 = ceil(Y ./ P(1));

Nt_2 = ceil(Y ./ P(3));

NDT_1 = zeros(1,4);

NDT_2 = zeros(1,4);

for zz = 1:4

NDT_1(zz) = max(max(max(ceil(Nt_1 ./ NR(zz)))));

NDT_2(zz) = max(max(max(ceil(Nt_2 ./ NR(zz)))));

end

AA = max(max(NDT_1))

BB = max(max(NDT_2))

CC = max(AA,BB)

[x,y] = max(abs((1./h).^-1));

fprintf('Overall Heat transfer = %6.4f \n',x)

if i < 4

fprintf('Pitch = %6.4f m \n',P(1))

fprintf('Face area = %6.4f m2 \n',FA(1))

fprintf('No. of fins per meter = %4.2f \n',sigma(1))

else

fprintf('Pitch = %6.4f m \n',P(3))

fprintf('Face area = %6.4f m2 \n',FA(3))

fprintf('No. of fins per meter = %4.2f \n',sigma(3))

end

if y < 4

fprintf('No. of rows = %2d \n',NR(1))

elseif 3 < y < 7

fprintf('No. of rows = %2d \n',NR(2))

elseif 6 < y < 10

fprintf('No. of rows = %2d \n',NR(3))

else

fprintf('No. of rows = %2d \n',NR(4))

end

fprintf('No. of tubes per row= %4d \n',CC)