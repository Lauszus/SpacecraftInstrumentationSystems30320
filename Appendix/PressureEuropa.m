%% Pressure on Europa 
clear all 
close all

rho_ice = 917; % Desity of ice [kg m^-3]
G = 6.67e-11; % Gravitational Constant [m^3 kg^-1 s^-2]
M =  4.8e22; % Mass of Europa [kg]
r = 1569e3; % Radius of Europa [m]
g = G*M/r^2; % Gravitatioanl acceleration on Europa's surface [m s^-2]

h = 0.1e3:0.1e3:112e3; % Depth [m]
P = rho_ice * g * h; % Pressure as funtion of depth [Pa]

%Pressure, when radius of Europa decreases with depth
g2 = (G*M)./(r-h).^2;
P2 = rho_ice.*g2.*h;

%Pressure, when radius and mass of Europa decrese with depth
g3 = (G*(M-(4/3*pi*(r^3-(r-h).^3)*1000)))./(r-h).^2;
P3 = rho_ice.*g3.*h;

%Plot pressure for all three scenarios
plot(h/1e3,P/1e6,h/1e3,P3/1e6)
xlabel('Depth below Ice[km]','fontweight','bold','fontsize',12)
ylabel('Pressure [MPa]','fontweight','bold','fontsize',12)
grid on
legend('g Constant','g varying with r')
