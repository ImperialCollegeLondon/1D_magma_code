

x1=50.5; y1=2772.04;   x2=70.625; y2=2422.60;
rhof_2=(y2-y1)/(x2-x1)*(47-x1)+y1; %Density fluid, least evolved at 47% 
rhof_1=(y2-y1)/(x2-x1)*(74-x1)+y1; %Density fluid, most evolved AT 74%
%  rhof_2=2600; %Density fluid, least evolved 
% rhof_1=2200; %Density fluid, most evolved
%Solid
x1=51.094; y1=2959.43;  x2=70.99; y2=2660.13;
rhom_2=(y2-y1)/(x2-x1)*(47-x1)+y1; %Density solid, least evolved
rhom_1=(y2-y1)/(x2-x1)*(74-x1)+y1; %Density solid, most evolved


Fontsize=13;
figure(1316)
clf; hold on
set(gcf,'color','w')
set(gca,'fontsize',13)

N=100;
x=linspace(47,74,N);
x1=50.5; y1=2772.04;   x2=70.625; y2=2422.60;


y0=(y2-y1)/(x2-x1)*(45-x1)+y1;
y=(x-45)*(y2-y1)/(x2-x1)+y0;

Melt_d=[y(1) y(end)];
plot(x,y,'linewidth',2,'color','r')

x1=51.094; y1=2959.43;  x2=70.99; y2=2660.13;

y0=(y2-y1)/(x2-x1)*(45-x1)+y1;
y=(x-45)*(y2-y1)/(x2-x1)+y0;

plot(x,y,'linewidth',2,'color','b')


xlabel('SiO_2 (%)','FontSize',Fontsize)
ylabel('Densities (kg/m^3)','FontSize',Fontsize)



%%

syms m1VM m2VM m3VM m4VM


% 3 reference compositions
% 1: 47% SiO2, 0% H2O,    1-0.47-m1      Cop2
% 2: 47% SiO2, 2% H2O,    1-0.47-0.02-m1 Cop2
% 3: 74% SiO2, 0% H2O,    1-0.74-m2      Cop2
% 4: 74% SiO2, 7% H2O,    1-0.74-0.07-m2 Cop2
% Third component M/V coefficient VM=f(x, y)=axy+bx+cy+d

M1=[0.47  0          1-0.47]; %kg
% M2=[0.47  0.02       1-0.47-0.02]; %kg
M2=[0.47  0.05       1-0.47-0.02]; %kg
M3=[0.74  0          1-0.74];
M4=[0.74  0.07       1-0.74-0.07];

VM1=[26.86e-6/0.06009   26.27e-6/0.01802 m1VM];
VM2=[26.86e-6/0.06009   26.27e-6/0.01802 m2VM];
VM3=[26.86e-6/0.06009   26.27e-6/0.01802 m3VM];
VM4=[26.86e-6/0.06009   26.27e-6/0.01802 m4VM];

% V=[M1.*VM; M2.*VM; M3.*VM; M4.*VM];
V=[M1.*VM1; M2.*VM2; M3.*VM3; M4.*VM4];

% v1=26.86e-6*Mol1;
% v2=26.27e-6*Mol2;

% Den=1./(sum(V,2))-[2665; 2540; 2360; 2110]*1;
Den=1./(sum(V,2))-[2665; 2400; 2360; 2110]*1.02;

Jaco=jacobian(Den,[m1VM m2VM m3VM m4VM]);

x=[2e-4 3e-4 4e-4 5e-4]'; %initial guess
Maxiter=100;
err=1;
iter=1;
while(err>1e-12 && iter<Maxiter)
    F=subs(Den, [m1VM m2VM m3VM m4VM], x');
    J=subs(Jaco, [m1VM m2VM m3VM m4VM], x');
    dx=-J\F;
    x=double(x+dx);

    err=abs(dx(1))+abs(dx(2))+abs(dx(3))+abs(dx(4));
    iter=iter+1;
end

% VM1=x(1);
% VM2=x(2);
% VM3=x(3);
% VM4=x(4);


X=[M1(1) M2(1) M3(1) M4(1)]';
Y=[M1(2) M2(2) M3(2) M4(2)]';
Z=x;

A=[X.*Y X Y ones(size(X))];
coef=A\Z;

%%



N=10;
x=linspace(0.47,0.74,N);



H2O=[0, 0.012,   0.02,  0.05,  0.07,   0.10]';

VMx=x.*H2O*coef(1)+repmat(x,length(H2O),1)*coef(2)+repmat(H2O,1,length(x))*coef(3)+coef(4);

x_other=1-x-H2O;


V1=x.*26.86e-6/0.06009+H2O.*26.27e-6/0.01802+x_other.*VMx;
Den=1./V1;
%%

handel=plot(x*100,Den,'--','linewidth',2);

Legned=cell(1,length(H2O)+2);
Legned{1}='Melt old';
Legned{2}='Solid';
for i=1:length(H2O)
    Legned{i+2}=[num2str(H2O(i)*100) '%'];
end


Color=get(handel,'Color');
plot(M1(1)*100, Den(1,1),'X','MarkerSize',15,'LineWidth',2,'Color',Color{1})
plot(M2(1)*100, Den(4,1),'X','MarkerSize',15,'LineWidth',2,'Color',Color{4})
plot(M3(1)*100, Den(1,end),'X','MarkerSize',15,'LineWidth',2,'Color',Color{1})
plot(M4(1)*100, Den(5,end),'X','MarkerSize',15,'LineWidth',2,'Color',Color{5})
legend(Legned,'fontsize',Fontsize)
% Plot()
%%
% aa=[-112.528, 127.811, 112.04];
% bb=[-0.381, -1.135, -0.411];
% cc=[0.033, 0, 0];
% 
% T=linspace(600,1100,300);
% P2=linspace(3000,7000,300);
% rho_temp=zeros(300);
% for i=1:300
%     for j=1:300
%         rho_temp(j,i)=(aa(1)* max(T(i),300).^bb(1)+aa(2)*max(P2(j),400).^bb(2)+aa(3)*max(T(i),300).^bb(3).*max(P2(j),400).^cc(1))*1000;
%     end
% end
% 
% figure(1317)
% clf; hold on
% set(gcf,'color','w')
% set(gca,'fontsize',13)
% 
% [X,Y] = meshgrid(T,P2/1000);
% surf(X,Y,rho_temp,'EdgeColor','none')
% % view([0 90])
% A=colorbar;
% A.Label.String = 'Density(kg^3)';
% A.Label.FontSize = Fontsize;
% xlabel('Temperature (C)','FontSize',Fontsize)
% ylabel('Pressure (kbar)','FontSize',Fontsize)