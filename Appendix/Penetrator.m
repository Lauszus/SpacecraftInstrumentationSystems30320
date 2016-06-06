%% Outer pressure of Penetrator. Solving for a specifc pressure
clear all 
close all 
penetrator = 1;
submarine = 0; 

%Choose metal used for the penetrator 
Titanium = 1; 
Inconel = 0; 
Aluminium = 0; 

if Titanium == 1
    ys = 827e6; %yield strength Pa
    rho = 4430; % density, kg/m^3
    young = 114e9; %Youngs modolus, Pa
elseif  Inconel == 1
    ys = 951e6; % yeild strength pa
    rho = 7850; % density kg/m^3
elseif Aluminium == 1
    ys = 241e6; % yeild strength pa
    rho = 2700; % density kg/m^3
end    
fos = 2; %factor of safety

% Setting up equation for pressure on a cylinder 
syms p_o p_i ro ri %variables
sigma_max(p_o,p_i,ro,ri) = -2* p_o*ro^2/(ro^2-ri^2) + p_i*(ro^2 + ri^2)/(ro^2-ri^2);
%p_o = outer pressure, p_i = inner pressure, ro = outer radius, ri = inner radius, sigma_max = maximum stress

%Solving equation for the inner radius 
ri(p_o,p_i,ro) = solve(sigma_max(p_o,p_i,ro,ri)==-ys/fos,ri); 

%Defining variables for penetrator and submarine
if penetrator == 1
    p_o = 12e6; %Outer pressure set to 12MPa
    p_i = 0; %Inner pressure set to zero
    ro = 0.1; %Outer radius
    
    % Calculating inner radius
    ri = ri(p_o,p_i,ro); 
    ri = double(ri(1));
    
    %Thickness of penetrator calculated when know inner and outer radius
    thickness = ro-ri;

    %Volume and mass of the cylindrical part of the penetrator
    h = 1.4;
    V_shell_cyl = pi*(ro^2-ri^2)*h;
    Vi_cyl = pi*ri^2*h;
    m_shell_cyl=V_shell_cyl*rho;

    % Volume and mass of tip of penetrator if assuming it is an sphere
    V_shell_sphere = 4/3*pi*(ro^3-ri^3);
    Vi_sphere = 4/3*pi*ri^3;
    m_shell_sphere = V_shell_sphere*rho;

    %Total mass and volume of penetrator 
    V_shell_total = V_shell_sphere + V_shell_cyl; 
    Vi_total = Vi_sphere + Vi_cyl; 
    m_shell_total = m_shell_cyl+m_shell_sphere;

    V_total = V_shell_total + Vi_total;
    
    strain = (ys/2)/young*100; % deformation in percentage
    
elseif submarine == 1
    p_o = 100e6; %Outer pressure set to 12MPa
    p_i = 0; %Inner pressure set to zero
    ro = 0.03; %Outer radius
    
    % Calculating inner radius
    ri = ri(p_o,p_i,ro); 
    ri = double(ri(1));

    %Thickness of penetrator calculated when know inner and outer radius
    thickness = ro-ri;

    %Volume and mass of the cylindrical part of the penetrator
    h = 0.15;
    V_shell_cyl = pi*(ro^2-ri^2)*h;
    Vi_cyl = pi*ri^2*h;
    m_shell_cyl=V_shell_cyl*rho;

    % Volume and mass of tip of penetrator if assuming it is an sphere
    V_shell_sphere = 4/3*pi*(ro^3-ri^3)*2;
    Vi_sphere = 4/3*pi*ri^3*2;
    m_shell_sphere = V_shell_sphere*rho;

    %Total mass and volume of penetrator 
    V_shell_total = V_shell_sphere + V_shell_cyl; 
    Vi_total = Vi_sphere + Vi_cyl; 
    m_shell_total = m_shell_cyl+m_shell_sphere;

    V_total = V_shell_total + Vi_total;
    
    strain = (ys/2)/young*100; % deformation in percentage
end 
%% Calculation of buoyancy 
 
% Calculating buoyancy for penetrator without balloon
M_all = 180; %Maximum mass of the penetrator [kg] 
g_europa = 1.314;% gravitational constant [m s^-2]
rho_water = 1000; % density of water [kg m^-3]
F_bouyancy = rho_water*g_europa*V_total; % Buoyancy for penetrator [N]
F_gravity = g_europa * M_all; % Gravitational force [N]
F_diff = M_all - (rho_water * V_total); 

%Bouynacy needed with a balloon (calculated by setting F_gravity equal
%F_buoyancy
% V_balloon = F_diff/(rho_water*g_europa); %m3
V_needed = M_all/rho_water;
V_balloon = V_needed - V_total;

%Example where CO2 is used as gas by using ideal gas law, PV=nRT
T_water = 273; %Temperature of water [K] 
molecularweight_co2 = 44; % Molecular Weight of CO2 [g mole^-1] 
R = 0.0831; %Ideal gas constant [L atm mole^-1 K^-1]
P_water = 0.35e6*9.87e-6; %Pressure of water in ballon [atm]

m_co2_liter = molecularweight_co2*P_water/(R*T_water);% Mass of CO2 per liter [kg liter^-1]
m_co2_total = V_balloon * 1000 * m_co2_liter; %Total mass of CO2 [kg}

%% Radiation 
%Constants for Plutonium decaying as alpha decay
c = 299792458;
Units = 1.66053904020e-27; % Atomic unit 
Pu238 = 238.049559894*Units; %Mass of plutonium238 [kg] 
U234 = 234.040952088*Units; %Mass og Uranium234 [kg]
Alpha = 4.00150646649*Units; %Mass og alpha particle [kg]

%Energy for each decaying atom
E_atom = (Pu238-(U234+Alpha)) * c^2; %[J]

% Energy produced for 1 kg of plutonium each second
WantedPlutonium = 1; % [kg]
A0 = WantedPlutonium/Pu238; %Number of atoms
t = 1; %Time [s]
T = 87.7*365.25*24*3600; %Half-life of Pu238 [s]
A = A0*2^(-t/T); %Number og particles left after the time, t, has passed 

% Energy each second for 1 kg of plutonium-238
E_sec = E_atom*(A0-A); %[J]

% Mass needed to produce 2kW after 7 years in space including 8% of 
%thermal energy used to produce electricity and assuming PuO2
%Ratio of Pu238 in PuO2
pu_2_puo2 = (238-32)/238;

%Extra energy due to electric energy generation
elec_ener = (2e3-150)/2e3;

%Energy for 1 atom when compensating for the above
E_atom2 = E_atom*pu_2_puo2*elec_ener; 

%Calculate thermal energy generation after 7 years for 1 kg
t2 = 7*365*24*3600; % 7 years calculated to seconds
A02 = A0*2^(-t2/T); % Number of atoms after 7 years
A2 = A02*2^(-t/T); % Number of atoms decaying each second after 7 years
E_sec2 = E_atom2*(A02-A2); %Energy produced for each kg

%Mass needed to produce 2kW
P = 2e3; %thermal energy wanted to melt [W]
m_pu = P/E_sec2; %mass needed [kg]

%Dimensions of the PuO2
h_rtg = 0.3; %Height of PuO2 [m]
rho_PuO2 = 11.46e3; %Density of PuO2 [kg/m^3]
v_rtg = m_pu/rho_PuO2; % volume of PuO2[m^3]

%Cylindrical holes where water flows to take up the energy produced by the PyO2
r_holes = 0.001; %radius of holes [m]
v_holes = pi*r_holes^2*h_rtg; %volume of holes [m^3]

%Radius of PuO2 when 20 holes are inside the volume of the RTG
n_hot = 20; %number of holes
r_rtg = sqrt((v_rtg+v_holes*n_hot)/(pi*h_rtg)); %radius of PuO2 [m]

%Heating of water at RTG
% A_circumference_rtg = 2*pi*r_rtg*h_rtg + n*2*pi*r_holes*h_rtg; %[m^2]
A_surface_rtg = n_hot*r_holes^2*pi + (r_rtg*1.2-(r_rtg))^2*pi; %Area of flowing water [m^2]
rho_water = 1000;% Density of water [kg m^-3]
vel_rtg = 0.03; %Velocity of waterflow[m s^-1]
Cp_water = 4.219e3; %Heat capaity for water [j kg^-1 K^-1]

flow_water = A_surface_rtg * vel_rtg * rho_water; %Flow of water [kg s^-1]

dT_hot = P/(flow_water*Cp_water) %Change in temperature [K]

% Cooling of water
%Formula for convective: q = h * A * dT, (A = C*l)

h_water = 3000;  %Convective Heat Transfer Coefficient [W m-^2 K^-1]

D_pipe = 2*r_rtg*pi/20; %Diameter of pipe [m]
C_pipe = pi*D_pipe; % Circumference of pipe [m]
dT_cool = (dT_hot)/2; % Average temperature of water while cooling
n_cool = 6; %number of cooling pipes
length_pipe = (P/n_cool)/(h_water*C_pipe*dT_cool)% length of pipe for all energy is dissipated

T_0 = 273; 
x = 0.001:0.001:0.4;
vel_pipe = flow_water/((D_pipe/2)^2*pi*n_cool*rho_water);
T =  exp(-((2.*h_water.*x)./(D_pipe/2*Cp_water*rho_water*vel_pipe))).*dT_hot; 

plot(x,T)
grid on
ylabel('Temperature difference [K]','fontweight','bold','fontsize',12)
xlabel('Length of pipe [m]','fontweight','bold','fontsize',12) 

%% Thermocouple
% Dimension of pellets: 1.4mm*1.4mm*1.12mm 
%Power from 127 pellets
P_127 = 7; % W 

%Total area of pellets: 
P_elec = 150; % Electric energy needed
A_pellet = P_elec/P_127*127*0.0014*0.0014 %Area needed for thermo couples [m^2]
A_rtg = r_rtg*1.2*2*pi*h_rtg %Total circumference of rtg [m^2]

A_pellet/A_rtg;

%% Calculation of the mass of the heat pipes
%Inconel
ys = 170e6; % yield strength Mpa
rho = 8470; % density kg/m^3

fos = 2; %factor of safety

syms p_o p_i ro ri
sigma_max(p_o,p_i,ro,ri) = -2* p_o*ro^2/(ro^2-ri^2) + p_i*(ro^2 + ri^2)/(ro^2-ri^2);

ri(p_o,p_i,ro) = solve(sigma_max(p_o,p_i,ro,ri)==ys/fos,ri);
p_o = 0;
p_i = 5e6;
ro = D_pipe/2; 

ri = ri(p_o,p_i,ro);
ri = double(ri(1));
thickness_pipe = ro-ri;

%Design of heat pipes: 18 vertical pipes: 6 high, 12 low. 17 toros needed
h_high = 1.4;

V_pipe = pi*(ro^2-ri^2)*h_high*6;
Vi_pipe = pi*ri^2*h_high*6;

V_total = V_pipe+Vi_pipe;
m_pipe=V_pipe*rho

%% Given two inner radii calculate the stree
% clear all 
%  
% p0=80e6;
% ro = 0.2; 
% ri = [0.1 0.15];
% for i = 1:length(ri);
%     r(:,i) = ri(i):(ro-ri(i))/9:ro;
%     for j = 1:10;
%         sigma_r(i,j) = p0*ro^2/(ro^2-ri(i)^2)*(ri(i)^2/r(j,i)^2-1);
%         sigma_theta(i,j) = -p0*ro^2/(ro^2-ri(i)^2)*(ri(i)^2/r(j,i)^2+1);
%     end 
% end
% 
% sigma_total = sigma_r + sigma_theta;
% for k = 1:length(ri)
%     figure
%     plot(r(:,k),sigma_r(k,:),r(:,k),sigma_theta(k,:),r(:,k),sigma_total(k,:))
%     legend('\sigma_r','\sigma_\theta','\sigma_{total}')
%     title(['ro = 0.2m and ri =' num2str(ri(k)) 'm and Pressure =' num2str(p0) 'Pa'])
%     xlabel('Radius')  
%     ylabel('Stress')
% end 
% 
% %% External pressure for one value
% clear all 
% 
% p0=15e6;
% ro = 0.1; 
% ri = 0.09; 
% r = ri:(ro-ri)/9:ro;
% sigma_r_ex = p0*ro^2/(ro^2-ri^2)*(ri^2./r.^2-1);
% sigma_theta_ex = -p0*ro^2/(ro^2-ri^2)*(ri^2./r.^2+1);
% sigma_total_ex = sigma_theta_ex+sigma_r_ex;
% 
% p_int = 10e6;
% sigma_r_int = p_int*ri^2*(1-ro^2./r.^2)/(ro^2-ri^2);
% sigma_theta_int = p_int*ri^2*(1+ro^2./r.^2)/(ro^2-ri^2);
% sigma_total_int = sigma_theta_int+sigma_r_int;
% 
% sigma_theta= sigma_theta_int+sigma_theta_ex;
% sigma_r = sigma_r_int + sigma_r_ex;
% sigma_total = sigma_theta+sigma_r;
% 
% subplot(1,2,1)
% plot(r,sigma_r_ex,r,sigma_theta_ex,r,sigma_total_ex)
% legend('\sigma_r','\sigma_\theta','\sigma_{total}')
% title(['ro = 0.2m and ri =0.95m and Pressure =' num2str(p0) 'Pa'])
% xlabel('Radius')  
% ylabel('Stress')
% subplot(1,2,2)
% plot(r,sigma_r_int,r,sigma_theta_int,r,sigma_total_int)
% legend('\sigma_r','\sigma_\theta','\sigma_{total}')
% title(['ro = 0.1m and ri =0.09m and Pressure =' num2str(p_int) 'Pa'])
% xlabel('Radius')  
% ylabel('Stress')
% 
% figure
% plot(r,sigma_r,r,sigma_theta,r,sigma_total)
% legend('\sigma_r','\sigma_\theta','\sigma_{total}')
% title(['ro = 0.1m and ri =0.09m and Pressure =' num2str(p_int) 'Pa'])
% xlabel('Radius')  
% ylabel('Stress')
% 
% h = 0.4;
% V = pi*(ro^2-ri^2)*h;
% m=V*4.5*1000;

