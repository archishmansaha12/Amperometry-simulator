%% AMPEROMETRY SIMULATOR %%

clear, clc;

%% INDEPENDENT VARIABLES %%
C      = 1.0E-3;  % [=] mol/L, initial concentration of O. Default = 1.0
D      = 1E-5;    % [=] cm^2/s, O & R diffusion coefficient. Default = 1E-5
V      = -600E-3; % [=] V, Voltage applied. Default = 1E-3
n      = 1.0;     % [=] number of electrons transfered. Default = 1
alpha  = 0.5;     % [=] dimensionless charge-transfer coefficient. Default = 0.5
k0     = 10E-2;   % [=] cm/s, electrochemical rate constant. Default = 1E-2
kc     = 1E-3;    % [=] 1/s, chemical rate constant. Default = 1E-3
T      = 298.15;  % [=] K, temperature. Default = 298.15

tk     = 100;    % [=] s, total time of measurement
t_app  = 30;     % [=] s, when reactant is added to solution
Dt     = 0.5;    % [=] s, intervals of measurement

%% PHYSICAL CONSTANTS %%
F      = 96485;   % [=] C/mol, Faraday's constant
R      = 8.3145;  % [=] J/mol-K, ideal gas constant
f      = F/(R*T); % [=] 1/V, normalized Faraday's constant at room temperature

%% SIMULATION VARIABLES %%
L      = tk/Dt;    % [=] number of instances in measurement
L1     = t_app/Dt; % [=] instances after which reactant added
DM     = 0.45;   % [=] model diffusion coefficient

%% DERIVED CONSTANTS %%

Dx  = sqrt(D*Dt/DM);      % [=] cm, delta x
j   = ceil(4.2*L^0.5)+5;  % number of boxes. If L~200, j=65


%% PRE-INITIALIZATION %%
C = C / 1000;           % Convert C from mol/L to mol/cm3
k = 0:L;                % time index vector
t = Dt * k;             % time vector
kf = k0*exp(  -alpha *n*f*V) % [=] cm/s, fwd rate constant
kb = k0*exp((1-alpha)*n*f*V) % [=] cm/s, rev rate constant

O = C*ones(L+1,j); % [=] mol/cm^3, concentration of O
R = zeros(L+1,j);  % [=] mol/cm^3, concentration of R
JO = zeros(1,L+1); % [=] mol/cm^2-s, flux of O at the surface

%% START SIMULATION %%
% i1 = time index. i2 = distance index
for i1 = L1:L
    % Update bulk concentrations of O and R
    for i2 = 2:j-1
        O(i1+1,i2) = O(i1,i2) + DM*(O(i1,i2+1)+O(i1,i2-1)-2*O(i1,i2));

        R(i1+1,i2) = R(i1,i2) + DM*(R(i1,i2+1)+R(i1,i2-1)-2*R(i1,i2));% ...
            %- km * R(i1,i2);
    end

    % Update flux
    JO(i1+1)   = ( kf*O(i1+1,2) - kb*R(i1+1,2) ) ./ (1 + Dx/D*(kf + kb) );

    % Update surface concentrations
    O(i1+1,1) = O(i1+1,2) - JO(i1+1)*(Dx/D);
    R(i1+1,1) = R(i1+1,2) + JO(i1+1)*(Dx/D);% - km*R(i1+1,1);
end

% Calculate current density, Z, from flux of O
Z = n.*F.*JO * 1000; % [=] A/cm^2 -> mA/cm^2, current density

%% PLOT RESULTS %%

figure;
plot(t,Z);
grid on;
xlabel('Time (s)'); 
ylabel('Current density (mA/cm^2)');
title('Amperometric readings for given initial concentration of M^+ at V = -600mV');