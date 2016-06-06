clear all
close all


E=100000; %Mpa
rho=4.51*10^3;  %kg/m^3

%m=827;  %Mpa, yield stress (not to go over)
m=410; % Factor of safety of 2
vect=0.01:0.005:0.19;
vect=vect*100;
thickness=[];
weight=[];
innvol=[];
masshe=[];

Rg=8.31446;
A=6.02214086*10^23;
M=6.642156*A*10^-27;

Pb=80;  %Mpa, external pressure
Pb2=0.1; %Pressure in the laboratory (penetrator)
Pa=0.1; % Internal pressure in the submarine in Mpa


for k=0.01:0.005:0.19

rb=k; %m
ra=rb; %m
crit=1000;


while crit>m

ra=ra-0.0001;
R=(ra:0.001:rb);

for i=1:length(R)

%Stress at 50 km depth
sigmarr(i)=((Pa*ra^3-Pb*rb^3)/(rb^3-ra^3))-(((Pa-Pb)*rb^3*ra^3)/((rb^3-ra^3)*R(i)^3));
sigmatt(i)=((Pa*ra^3-Pb*rb^3)/(rb^3-ra^3))+(((Pa-Pb)*rb^3*ra^3)/(2*(rb^3-ra^3)*R(i)^3));
stress(i)=abs(sigmarr(i)-sigmatt(i));

%Stress in the pressured cylinder
sigmarr2(i)=((Pa*ra^3-Pb2*rb^3)/(rb^3-ra^3))-(((Pa-Pb2)*rb^3*ra^3)/((rb^3-ra^3)*R(i)^3));
sigmatt2(i)=((Pa*ra^3-Pb2*rb^3)/(rb^3-ra^3))+(((Pa-Pb2)*rb^3*ra^3)/(2*(rb^3-ra^3)*R(i)^3));
stress2(i)=abs(sigmarr2(i)-sigmatt2(i));


end

crit=max([stress stress2]);

end

thickness=[thickness rb-ra];  %in m
w=rho*4*pi*(rb^3-ra^3)/3; %in kg
weight=[weight w];
vo=4*pi*ra^3/3;
innvol=[innvol vo];

rhohe=M*Pa*10^6/(Rg*273); %kg m-3
mahe=rhohe*4*pi*ra^3/(3);

masshe=[masshe mahe];

end


% figure
% hold on
% plot(R,sigmarr,'b',R,sigmatt,'r',R,stress,'y')
% xlabel('radius (m)')
% ylabel('Stress (Mpa)')
% title(['Wall thickness of ',num2str(thickness),' m and titanium mass of ',num2str(weight),' kg'])

figure
hold on
plot(vect,thickness*1000,'b')
xlabel('Outer radius(cm)')
ylabel('Thickness (mm)')
title('Evolution of the thickness of the titanium shell with the radius')

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'Thickness','-dpdf')

figure
hold on
plot(vect,weight,'b',vect,weight+masshe,'r')
xlabel('Outer radius(cm)')
ylabel('Weight (kg)')
title('Weight of the submarine')

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'SphereWeight','-dpdf')
