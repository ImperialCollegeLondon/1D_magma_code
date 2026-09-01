% SiO2 TiO2 Al2O3 FeO) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% Composition = [64.6  0.53 16.5  4.47 0.01 2.39 5.23 4.49 1.54 0.04 0.01 0];
syms m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12

bb  = [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  = [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

siti    =   m1+m2;   %xmf_t(:,1) + xmf_t(:,2);
tial    =   m2+m3;   %xmf_t(:,2)+xmf_t(:,3);
fmm     =   m4+m5+m6;%xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,6);
nak     =   m8+m9;  %xmf_t(:,8) + xmf_t(:,9);
b1  =   siti;
b2  =   m3;% xmf_t(:,3);
b3  =   m4+m5+m10; %xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,10);
b4  =   m6; %xmf_t(:,6);
b5  =   m7; %xmf_t(:,7);
b6  =   m8+m11+m12; %xmf_t(:,8) + xmf_t(:,11) + xmf_t(:,12);
b7  =   m11+m12+log(1+m11); %xmf_t(:,11) + xmf_t(:,12) + log(1+xmf_t(:,11));
b12 =   siti.*fmm;
b13 =   (siti + m3 + m10).*( nak + m11 ); %(siti + xmf_t(:,3) + xmf_t(:,10)).*( nak + xmf_t(:,11) );
b14 =   m3.*nak; %xmf_t(:,3).*nak;

c1      =   m1; %xmf_t(:,1);
c2      =   tial;
c3      =   fmm;
c4      =   m7; %xmf_t(:,7);
c5      =   nak;
c6      =   log(1+m11 + m12); %log(1+xmf_t(:,11) + xmf_t(:,12));
c11     =   m3+fmm+m7-m10; %xmf_t(:,3) + fmm + xmf_t(:,7) - xmf_t(:,10);
c11     =   c11.*(nak+m11+m12); %c11.*(nak + xmf_t(:,11) + xmf_t(:,12));

bcf      =   [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf      =   [c1 c2 c3 c4 c5 c6 c11];

BT          = sum(bb.*bcf(:,:),2);
CT          = sum(cc.*ccf(:,:),2);

% Lmu=-4.55+BT/(T-CT);
% 
% dLmud1=diff(Lmu,m1);
% dLmud2=diff(Lmu,m11);
% dBd1=diff(BT,m1);
% dBd2=diff(BT,m11);

B=subs(BT,[m2 m3 m4 m5 m6 m7 m8 m9 m10 m12], [  0.53 16.5  4.47 0.01 2.39 5.23 4.49 1.54 0.04  0]);
C=subs(CT,[m2 m3 m4 m5 m6 m7 m8 m9 m10 m12], [  0.53 16.5  4.47 0.01 2.39 5.23 4.49 1.54 0.04  0]);

%%
% x: SiO%
% y: H2O%
% B=b1*x+b2*y+b3*log(y+1)+b4+b5*(x+b6)*(y+b7)
% C=c1*x+c2*y+c3*log(y+1)+c4

% Solve the system
% using A=-4.55
syms b1 b2 b3 b4 c1 c2 c3 c4 x y T
B=b1*x+b2*y+b3*log(y+1)+b4;
C=c1*x+c2*y+c3*log(y+1)+c4;
% Var=[b1 b2 b3 b4 c1 c2 c3];
Var=[b1 b2 b3 c1 c2 c3];  %variables used for fitting

logeta=-4.55+B/(T+273.15-C);
Jac=jacobian(logeta, Var); 
LHS=matlabFunction(Jac,'Vars', [b1 b2 b3 b4 c1 c2 c3 c4 x y T]);
RHS=matlabFunction(logeta,'Vars', [b1 b2 b3 b4 c1 c2 c3 c4 x y T]);
%% fitting data H2O% for least and most evolved component
% H2O_frac=[0, 0.2   0.5      1];  %from 0% to saturated
% fit_data=[1  2.33  2.66     3,...%least evolved
%           4  6.66  7.33     8];  %most evolved

P2=7000;  %bar

P3=P2/10; %Mpa
Saturation=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5);

% H2O=Saturation*H2O_frac;

%% H2O dependancy for solidus and liquidus
%(0,0)    : (57,326)
%(1200,0) : (455,326)
%(0,8000) : (57,70)
% from 8000 bar to 0 bar, water 13% to 0 %
%           13  12  11  10    9     8     7    6     5    4    3    2   1     0
Data_point=[80, 83, 86, 91, 100,  115,  135, 164,  196, 230, 267, 307, 355, 401;... %70
             0, 82, 85,90, 98, 111.5,130.5,158,  189, 220, 255, 293, 343, 388;... %102
             0,  0, 84,87.5,95.5,108.5,126, 151.5, 180, 211, 243, 280, 329, 375;... %134
             0,  0, 0,   86, 93, 105,  121, 144,   170, 199, 230, 266, 316, 362;... %166
             0,  0, 0,    0,90.5,100,  115, 136,   159, 186, 217, 252, 302, 349;... %198
             0,  0, 0,    0,   0, 95,  109, 128.5, 150, 176, 205, 241, 291, 339;... %223
             0,  0, 0,    0,   0,  0,  103, 120,   141, 165, 194, 230, 280, 328;... %248
             0,  0, 0,    0,   0,  0,    0, 113,  133.5,155.5,185,220, 270, 320;... %268
             0,  0, 0,    0,   0,  0,    0,   0,   127,146.5,176, 211.5,261,313;... %285
             0,  0, 0,    0,   0,  0,    0,   0,     0, 141, 170, 206, 256, 309;... %295
             0,  0, 0,    0,   0,  0,    0,   0,     0,   0, 164, 200, 250, 305;... %305
             0,  0, 0,    0,   0,  0,   0,    0,     0,   0,   0, 194, 244, 301;... %316
             0,  0, 0,    0,   0,  0,   0,    0,     0,    0,  0,   0, 240, 298;... %322
             0,  0,  0,   0,   0,  0,    0,   0,     0,   0,   0,   0,   0, 297];   %326   
Data_y=[70, 102,134,166,198,223,248,268,285,295,305,316,322,326];
Data_y=8000-(Data_y-70)*8000/(326-70);
Data_point=(Data_point-57)*(1200-600)/(455-57)+600;
%
% (0,0) (61,738)
% (0,1700) (838,738)
% (0,40000) (61 33)
%             20   10    5   2    0
Data_point2=[358, 428,  512, 586, 699;...%120
             321, 397,  477, 547, 653;...%208
             288, 369,  443, 508, 604;...%297
             274, 351,  412, 473, 558;...%390
               0, 340,  391, 447, 517;...%473
               0, 332,  379, 427, 477;...%561
               0, 329,  376, 420, 462;...%596
               0,   0,  374, 414, 438;...%691
               0,   0,    0, 413, 436;...%728
               0,   0,    0,   0, 436];...%738  
Data_y2=[120,208,297,390,473,561,596,691,728,738];               
Data_y2=40000-(Data_y2-120)*40000/(738-120);
Data_point2=(Data_point2-61)*(1700-800)/(838-61)+800;
%% Obtaining the temperature

A1=-99.6;  C1=1160; B1=760-C1-A1;
Ts0=A1+B1+C1;
Tl0=C1;

if length(Var)==4
    H2O=[0 Saturation 0 Saturation];
    SiO2=[47 47 74 74];
    fit_data=[2 0 7 3];  %order of melt viscosity
elseif length(Var)==6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The fitting settings!
    H2O=[0 Saturation*0.2 Saturation 0 Saturation*0.2 Saturation];
    SiO2=[47 47 47 74 74 74];
    %[least/wet, least/dry, most/wet, most/dry]
    Base_data=[0 2 4 7];
    drop_rate=0.3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mid1=Base_data(1)+(Base_data(2)-Base_data(1))*drop_rate;
    Mid2=Base_data(3)+(Base_data(4)-Base_data(3))*drop_rate;
    fit_data=[Base_data(2)  Mid1 Base_data(1) Base_data(4)  Mid2 Base_data(3)];  %order of melt viscosity
end

T=zeros(1,length(H2O));
for j=1:length(H2O)
    Ts=Cal_solidus(H2O(j),P2,Data_point,Data_y);
    Tl=Cal_liquidus(H2O(j),P2,Data_point2,Data_y2);
    a1=A1*(Tl-Ts)/(Tl0-Ts0);
    b1=B1*(Tl-Ts)/(Tl0-Ts0);
    c1=Tl;
    Cl=(SiO2(j)-47)/(74-47);
    T(j)=a1*Cl^2+b1*Cl+c1;    
end

%% Solve the parameters
N=length(Var);
% H2O=repmat(H2O,1,2);
% SiO2=[47*ones(1,4) 74*ones(1,4)];




step_size=0.02;
error=1;
iter=0;
% b1 b2 b3 b4 c1 c2 c3 c4 x y T
% initial guess

% Initial guess
b1=179.17189; b2=102.7844; b3=1182.7931; b4=0;
c1=-0.66143385;  c2=-37.328439; c3=-456.07622; c4=0;




lhs=zeros(N);
rhs=zeros(N,1);
in=[b1 b2 b3 b4 c1 c2 c3 c4]';
% x=[b1 b2 c1 c2];
while (error>1e-5 && iter <10000)
    for i=1:N
        lhs(i,:)=LHS(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),SiO2(i), H2O(i), T(i));   
        rhs(i)=  RHS(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),SiO2(i), H2O(i), T(i))-fit_data(i);   
    end
    dx=lhs\rhs;
%     b1=b1-dx(1)*step_size; b2=b2-dx(2)*step_size; b3=b3-dx(3)*step_size; b4=b4-dx(4)*step_size;
%     c1=c1-dx(5)*step_size; c2=c2-dx(6)*step_size; c3=c3-dx(7)*step_size; c4=c4-dx(8)*step_size;
    if N==4
        in([1 2 5 6])=in([1 2 5 6])-dx*step_size;
    else
        in([1 2 3 5 6 7])=in([1 2 3 5 6 7])-dx*step_size;
    end
    error=max(abs(dx(:)));
    iter=iter+1;
end


disp(in)
b1=in(1); b2=in(2);
c1=in(5); c2=in(6);
if N==6
    b3=in(3); c3=in(7);
end
%%
N=300;
Font=20;

sio2=linspace(47,74,N);
x=linspace(0,1,N);
h2o=linspace(0,10,N);

[X,Y] = meshgrid(sio2,h2o);
V=zeros(N);
for j=1:length(h2o)
    Ts=Cal_solidus(h2o(j),P2,Data_point,Data_y);
    Tl=Cal_liquidus(h2o(j),P2,Data_point2,Data_y2);
    aa1=A1*(Tl-Ts)/(Tl0-Ts0);
    bb1=B1*(Tl-Ts)/(Tl0-Ts0);
    cc1=Tl;

    T=aa1*x.^2+bb1*x+cc1;
    B=b1*sio2+b2*h2o(j)+b3*log(1+h2o(j))+b4;
    C=c1*sio2+c2*h2o(j)+c3*log(1+h2o(j))+c4;
    V(N-j+1,:)=-4.55+B./(T+273.15-C);
end
surf(X,Y,V,'EdgeColor','none')

index=find(Y(:,1)<1.21,1,'last');
plot3(X(index,:),Y(index,:),V(index,:),'color','red','linewidth',2)
text(45, 1.21, 0, 'P','fontsize', 20,'FontWeight','bold')
text(75.5, 1.41, 4.5, 'E','fontsize', 20,'FontWeight','bold')

xlabel('SiO_2(%)','fontsize',Font)
ylabel('H_2O(%)','fontsize',Font)
zlabel('Melt viscosity order (-)','fontsize',Font)
view([139.1000+180   30.4004])
% title(['(' alphabet(p) ')'] , 'fontsize',Font-1)
grid on
colorbar