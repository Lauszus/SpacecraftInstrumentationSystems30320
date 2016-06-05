clear; clc; close all;

differentialDisplay = 0;
%% Initialize constants
thermalData = dlmread('thermalData.txt');
C_S = 2050; C_L = 4224.5; % Heat capacities of ice and water at 273 K
rho_S = 916.2; rho_L = 999.9684; % Densities of ice and water at 273 K

rows = 51; columns = 60; layers = columns; % Grid size
maxIterations = 5e6; lastUpdate = 1;

iceThickness = 100; % [m]
Tm = 273; Ts = 100; % Melting and surface temperatures [K]
L = 333.5e3; % Latent Heat of Fusion [J/kg]
couplingCoeff = 0.99;
dz = .5; dr = .1; % GridSpacing [m]

% The different heater powers
heatings = [1000 1500 2000 3000 4500 6750 10000]; 
heaterWidth = 2; % Heater width in grid cells
heaterRow = ceil(rows/4*3); % Heater location in the grid
heaterColumn = round(columns/2) - ceil(heaterWidth / 2);
heaterLayer = round(layers/2) - ceil(heaterWidth / 2);

times = zeros(size(heatings));  % Penetration times
losses = zeros(size(heatings)); % Loss percentages

tic
for mode = 1%:2
    %% Temperature Profile data
    switch mode % 1 = linear; 2 = 2-layer
        case 1
            profileT = [Ts-1 Ts-1 Ts Tm Tm+1 Tm+1]; 
            profileD = [-1 -0.0001 0 1 1.0001 2]; % Linear
        case 2
            profileT = [Ts-1 Ts-1 Ts 240 Tm Tm+1 Tm+1]; 
            profileD = [-1 -0.0001 0 0.3 1 1.0001 2]; % 2-layer
    end
    
    for kk=3%1:length(heatings)       
        fprintf([num2str(heatings(kk)) ...
                 ' W Heater, Profile ' ...
                 num2str(mode) '\n'])
        
        heatingDensity = heatings(kk) / (heaterWidth^2 * dr^2 * dz); %W/m^3       
        dt = floor(3e6 / heatingDensity) * 500 * dz * dr^2; % Time step [s]
        time = 0; depth = 0;

        %% Initialize variables
        tempProfile = interp1(profileD,...
                              profileT,...
                              ((1:rows) - rows*3/4) * dz / iceThickness);
        T = repmat(tempProfile',1,columns,layers);
        newTemp = min(tempProfile(end),274);
        
        rho = interp1(thermalData(:,5),thermalData(:,6),T,...
                      'linear','extrap');
        k   = interp1(thermalData(:,5),thermalData(:,7),T,...
                      'linear','extrap');
        C   = interp1(thermalData(:,5),thermalData(:,8),T,...
                      'linear','extrap');
        
        LH      = zeros(rows,columns,layers);
        state   = zeros(rows,columns,layers);
        
        P       = zeros(rows,columns,layers);
        P(heaterRow-1,... % Place heater voxels
          heaterColumn:heaterColumn+heaterWidth-1,...
          heaterLayer:heaterLayer+heaterWidth-1)...
                = heatingDensity * (1-couplingCoeff);
        P(heaterRow,...
          heaterColumn:heaterColumn+heaterWidth-1,...
          heaterLayer:heaterLayer+heaterWidth-1)...
                = heatingDensity * couplingCoeff;    
        
        %% Create Plots
        fig = figure;
        f = imagesc(T(:,:,heaterLayer));%imshow(T,[100 273]);
        c = colorbar; cm = parula(1000);
        cm(end,:) = [1 0 0];
        cm(1,:)   = [0 0 0];
        colormap(cm);
        ylabel(c,'Temperature [K]')
        if differentialDisplay
            caxis([0 45]); 
        else
            caxis([Ts-1 Tm+1]); 
        end
        ylabel('Depth [m]');
        ax = gca; ax.XTick = linspace(1,columns,5);
        ax.XTickLabel = num2cell(abs(linspace(0,columns,5) - ...
                                 round(columns/2)) * dr);
        ax.YTickLabel = num2cell((str2double(ax.YTickLabel) - ...
                                 floor(rows*3/4)) * dz);
        xlabel('Radial Distance [m]'); axis equal tight;    
        suptitle('Time: 0.00 Days   Depth: 0 m   Ambient Temp: 100 K');
        
        %  waitforbuttonpress;
        %% Simulation Loop
        counter = 1;
        for i=1:maxIterations
            % Calculate temperature change. Note that the factor 6 is there
            % because del2(T) retuns nabla^2*T/6 for whatever reason
            T  = T + (6 * k .* del2(T,dr,dz,dr) + P) * dt ./ (rho .* C);
            
            % Handle latent heat
            idx = T > Tm & state == 0;
            LH(idx) = LH(idx) + (T(idx) - Tm) .* C(idx) .* rho(idx);
            T(idx) = Tm;
            
            % Handle melting
            idx = LH > 333.5e3 * rho_S;
            T(idx) = T(idx) + (LH(idx) - 333.5e3 * rho_S) / (rho_L * C_L);
            LH(idx) = 333.5e3 * rho_S;
            state(idx) = 1;
            
            % Handle freezing
            idx = state == 1 & LH < 0;
            T(idx) = T(idx) + LH(idx) / (rho_S * C_S);
            LH(idx) = 0;
            state(idx) = 0;
            
            % Get thermal properties
            % Water
            rho(logical(state)) = interp1(thermalData(:,1),...
                                          thermalData(:,2),...
                                          T(logical(state)),...
                                          'linear','extrap');
            k(logical(state))   = interp1(thermalData(:,1),...
                                          thermalData(:,3),...
                                          T(logical(state)),...
                                          'linear','extrap');
            C(logical(state))   = interp1(thermalData(:,1),...
                                          thermalData(:,4),...
                                          T(logical(state)),...
                                          'linear','extrap');
            %Ice
            rho(~logical(state)) = interp1(thermalData(:,5),...
                                           thermalData(:,6),...
                                           T(~logical(state)),...
                                           'linear','extrap');
            k(~logical(state))   = interp1(thermalData(:,5),...
                                           thermalData(:,7),...
                                           T(~logical(state)),...
                                           'linear','extrap');
            C(~logical(state))   = interp1(thermalData(:,5),...
                                           thermalData(:,8),...
                                           T(~logical(state)),...
                                           'linear','extrap');
            
            % Move window as ice melts
            if (state(heaterRow,heaterColumn,heaterLayer) == 1 && ...
                      depth <= iceThickness)
                depth = depth + dz;

                T       = circshift(T,-1,1);
                LH      = circshift(LH,-1,1);
                state   = circshift(state,-1,1);
                C       = circshift(C,-1,1);
                rho     = circshift(rho,-1,1);
                k       = circshift(k,-1,1);
                
                tempProfile = interp1(profileD,...
                                      profileT,...
                                      (((1:rows) - ceil(rows*3/4)) * dz...
                                      + depth) / iceThickness);
                                  
                newTemp = min(tempProfile(end),274);
                T(end,:,:)      = newTemp;
                LH(end,:,:)     = (333.5e3 * rho_S) * (newTemp > 273);
                state(end,:,:)  = (newTemp > 273);
                C(end,:,:)      = interp1(thermalData(:,5),...
                                          thermalData(:,8),...
                                          newTemp,'linear','extrap');
                rho(end,:,:)    = interp1(thermalData(:,5),...
                                          thermalData(:,6),...
                                          newTemp,'linear','extrap');
                k(end,:,:)      = interp1(thermalData(:,5),...
                                          thermalData(:,7),...
                                          newTemp,'linear','extrap');
                
                ax.YTickLabel = num2cell(str2double(ax.YTickLabel) + dz);
                
                % End simluation if water has been reached
                if depth >= iceThickness              
                    break;
                end                
            end
            
            %% Handle Boundaries for numerical stability
             T(:,1,:)   = repmat(tempProfile',1,columns);
             T(:,end,:) = repmat(tempProfile',1,columns);
             T(:,:,1)   = repmat(tempProfile',1,columns);
             T(:,:,end) = repmat(tempProfile',1,columns);
             T(end,:,:) = newTemp;
             
            % Update plots
            if i - lastUpdate > 100
                lastUpdate = i;
                if differentialDisplay
                    f.CData = mag2db(abs(T(:,:,heaterLayer) - ...
                                repmat(T(:,1,heaterLayer),1,50)));
                else
                    f.CData = min(T(:,:,heaterLayer),273) + ...
                                state(:,:,heaterLayer);
                end
                suptitle(['Time: ',...
                    num2str(i * dt / (24*3600),'%2.1f'),...
                    ' Days    Depth: ', ...
                    num2str(depth), ...
                    ' m    Ambient Temp: ',...
                    num2str(T(heaterRow,1,1), '%3.1f'),...
                    ' K']);
                drawnow limitrate;
            end
        end
        %% Calculate energy lost to lattice
        switch mode
            case 1
                minimumEnergyNeeded = (1/2 * (Tm-Ts)) * ...
                    iceThickness * C_S * rho_S * dr^2 * heaterWidth^2 + ...
                    L * iceThickness * heaterWidth^2 * dr^2 * rho_S;%Linear
            case 2
                minimumEnergyNeeded = 1/2 * ((Tm-Ts) * 0.3 + Tm - 240) *...
                    iceThickness * C_S * rho_S * dr^2 * heaterWidth^2 + ...
                    L * iceThickness * heaterWidth^2 * dr^2 * rho_S;%2layer
        end
        usedEnergy = heatingDensity * (heaterWidth^2 * dr^2 * dz) * i * dt;
        
        lossPercentage = (usedEnergy - minimumEnergyNeeded) / ...
                          usedEnergy * 100
        
        losses(mode, kk) = lossPercentage;
        times(mode, kk) = i * dt / (24 * 3600);
    end
end
toc

%% Plot the results
figure;
subplot(1,2,1);
hold on;
plot(heatings/(0.2^2)/1000,times*10);
hold off; grid on;
legend('Linear Profile', '2-Layer Profile');
ylabel('Melting Time [days/km]'); xlabel('Heat Generation [kW/m^2]');
title('Melting Time Profile');

subplot(1,2,2);
hold on;
plot(heatings/(0.2^2)/1000,losses);
hold off; grid on;
legend('Linear Profile', '2-Layer Profile');
ylabel('Energy Loss Percentage'); xlabel('Heat Generation [kW/m^2]');
title('Energy Loss Profile');
