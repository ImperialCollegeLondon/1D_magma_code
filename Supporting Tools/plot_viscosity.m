
Composition0 = [64.6  0.53 16.5  4.47 0.01 2.39 5.23 4.49 1.54 0.04 0.01 0];
Ratio=Composition0([2:10 12])/sum(Composition0([2:10 12]));

H2O=[0.0 0.5,1,2,3, 5];
Text=cell(1,length(H2O));
for i=1:length(H2O)
    Text{i}=['H2O=' num2str(H2O(i)) '%'];
end
SiO2=linspace(46,75,500);

T=1000;

N=length(SiO2);



Font=13;
figure(1)
clf;
subplot(2,2,1)
hold on; box off;
set(gca,'fontsize',Font)
T0=800;
T=T0+273.15;
V=Get_viscosity(H2O,SiO2,T,Ratio );

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
% legend(Text,'FontSize',Font)
title(['T=' num2str(T0) 'C'])

subplot(2,2,2)
hold on; box off;
set(gca,'fontsize',Font)
T0=900;
T=T0+273.15;
V=Get_viscosity(H2O,SiO2,T,Ratio );

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
legend(Text,'FontSize',Font)
title(['T=' num2str(T0) 'C'])

subplot(2,2,3)
hold on; box off;
set(gca,'fontsize',Font)
T0=1000;
T=T0+273.15;
V=Get_viscosity(H2O,SiO2,T,Ratio );

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
% legend(Text,'FontSize',Font)
title(['T=' num2str(T0) 'C'])

subplot(2,2,4)
hold on; box off;
set(gca,'fontsize',Font)
T0=1100;
T=T0+273.15;
V=Get_viscosity(H2O,SiO2,T,Ratio );

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
% legend(Text,'FontSize',Font)
title(['T=' num2str(T0) 'C'])
%%
H2O=linspace(0.01,15,100);

SiO2=linspace(47,74,6);
Text2=cell(1,length(SiO2));
for i=1:length(SiO2)
    Text2{i}=['SiO2=' num2str(SiO2(i)) '%'];
end


figure(4)
clf;
hold on; box off;
set(gca,'fontsize',Font)
T0=800;
T=T0+273.15;
V=Get_viscosity(H2O,SiO2,T,Ratio );
plot(H2O,V','linewidth',2)
xlabel('H_2O%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
legend(Text2,'FontSize',Font)
title(['T=' num2str(T0) 'C'])
%%
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




SiO2=linspace(47,74,500);
H2O=[0.2 0.5,1,2,3, 5];
Text=cell(1,length(H2O));
for i=1:length(H2O)
    Text{i}=['H2O=' num2str(H2O(i)) '%'];
end





figure(2)
clf; subplot(1,2,1)
hold on; box off;
set(gcf,'color','w')
set(gca,'fontsize',Font)
P2=7000;
V=Get_viscosity2(H2O,SiO2,P2, Ratio,Data_point,Data_y,Data_point2,Data_y2);

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
% legend(Text,'FontSize',Font)
title(['P=' num2str(P2/1000) 'kbar'])

subplot(1,2,2)
hold on; box off;
set(gca,'fontsize',Font)
T0=800;
P2=4000;
V=Get_viscosity2(H2O,SiO2,P2, Ratio,Data_point,Data_y,Data_point2,Data_y2);

plot(SiO2,V,'linewidth',2)
xlabel('SiO_2%','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
legend(Text,'FontSize',Font)
title(['P=' num2str(P2/1000) 'kbar'])
%%
SiO2=linspace(47,74,6);
H2O=linspace(0,12,N);
Text2=cell(1,length(SiO2));
for i=1:length(SiO2)
    Text2{i}=['SiO2=' num2str(SiO2(i)) '%'];
end

figure(3)
clf
hold on; box off;
set(gcf,'color','w')
set(gca,'fontsize',Font)
P2=7000;
V=Get_viscosity2(H2O,SiO2,P2, Ratio,Data_point,Data_y,Data_point2,Data_y2);

plot(H2O,V','linewidth',2)
xlabel('H2O %','FontSize',Font)
ylabel('Log10(\eta)','FontSize',Font)
legend(Text2,'FontSize',Font)
title(['P=' num2str(P2/1000) 'kbar'])




%%
function V=Get_viscosity2(H2O,SiO2,P2, Ratio,Data_point,Data_y,Data_point2,Data_y2)
A1=-99.6;  C1=1160; B1=760-C1-A1;

Ts0=760;
Tl0=1160;

N=length(SiO2);

V=zeros(length(H2O),N);
Composition=zeros(length(H2O),12);
for i=1:N
    Composition(:,1)=SiO2(i);
    Cl=(SiO2(i)-47)/(74-47);
    for j=1:length(H2O)
        Composition(j,11)=H2O(j);
        Composition(j,[2:10 12])=(100-sum(Composition(j,[1 11])))*Ratio; 
        
        Ts=Cal_solidus(H2O(j),P2,Data_point,Data_y);
        Tl=Cal_liquidus(H2O(j),P2,Data_point2,Data_y2);
        a1=A1*(Tl-Ts)/(Tl0-Ts0);
        b1=B1*(Tl-Ts)/(Tl0-Ts0);
        c1=Tl;
        
        T=a1*Cl^2+b1*Cl+c1;
        V(j,i)=Giordano_2008_visc(H2O(j), T, Composition(j,:));
    end    
end
end



function V=Get_viscosity(H2O,SiO2,T,Ratio )
N=length(SiO2);

V=zeros(length(H2O),N);
Composition=zeros(length(H2O),12);
for i=1:N
    Composition(:,1)=SiO2(i);
    for j=1:length(H2O)
        Composition(j,11)=H2O(j);
        Composition(j,[2:10 12])=(100-sum(Composition(j,[1 11])))*Ratio; 
        V(j,i)=Giordano_2008_visc(H2O(j), T, Composition(j,:));
    end    
end
end

function v = Giordano_2008_visc(H2Ot, T, Composition)
%Get the Giordano2008 VFT parameters from their script
%Giordano2008_Model.mat
outvalues=Giordano2008_Model(H2Ot, Composition);
At = outvalues(:,1);
Bt = outvalues(:,2);
Ct = outvalues(:,3);
v = At+(Bt./(T-Ct));
end

function outvalues=Giordano2008_Model(H2Ot, Composition)
%This is code is from Girodano 2008 (see citation below) modified
%for use by the Coumans et al., 2020 bubble growth model

% SCRIPT grdmodel08 Nov 2007
%
% MATLAB script to compute silicate melt viscosity.
%
% Citation: Giordano D, Russell JK, & Dingwell DB (2008)
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
%INPUT: Chemical compositions of silicate melts (Wt. % oxides) as:
% SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% One line for each melt composition
% _________________________________________________________________

% ________________________________________________________________
% OUTPUT: VFT Model parameters A, B and C for each melt composition
%         to model temperature dependence: log n = A + B/[T(K) - C]
% VFT_model.out contains: No. , A_value , B_value , C_value , Tg , F
% VFT_curves.out contains: Temperature Array in Kelvins &
%                          Log(n) values at T(K) values (1 line per melt)
% ________________________________________________________________

% VFT Multicomponent-Model Coefficients
% -----------------------------------------------------------------
AT      = -4.55;
bb  = [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  = [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

row = length(H2Ot);
Comp_rep = repmat(Composition,row,1);
Comp_rep(:,11) = H2Ot';

wtm=Comp_rep;

[nx,nc]=size(wtm);

bb_rep = repmat(bb,nx,1);
cc_rep = repmat(cc,nx,1);
AT_rep = repmat(AT,nx,1);

% Function molefrac_grd: converts wt % oxide bais to mole % oxide basis
%[ncomps xmf_t] = molefrac_grd(wtm);
[ncomps, xmf_t] = molepct_grd(wtm);

%outvalues=[];
% Load composition-basis matrix for multiplication against model-coefficients
% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
siti    =   xmf_t(:,1) + xmf_t(:,2);
tial    =   xmf_t(:,2)+xmf_t(:,3);
fmm     =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,6);
nak     =   xmf_t(:,8) + xmf_t(:,9);
b1  =   siti;
b2  =   xmf_t(:,3);
b3  =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,10);
b4  =   xmf_t(:,6);
b5  =   xmf_t(:,7);
b6  =   xmf_t(:,8) + xmf_t(:,11) + xmf_t(:,12);
b7  =   xmf_t(:,11) + xmf_t(:,12) + log(1+xmf_t(:,11));
b12 =   siti.*fmm;
b13 =   (siti + xmf_t(:,3) + xmf_t(:,10)).*( nak + xmf_t(:,11) );
b14 =   xmf_t(:,3).*nak;

c1      =   xmf_t(:,1);
c2      =   tial;
c3      =   fmm;
c4      =   xmf_t(:,7);
c5      =   nak;
c6      =   log(1+xmf_t(:,11) + xmf_t(:,12));
c11     =   xmf_t(:,3) + fmm + xmf_t(:,7) - xmf_t(:,10);
c11     =   c11.*(nak + xmf_t(:,11) + xmf_t(:,12));
bcf      =   [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf      =   [c1 c2 c3 c4 c5 c6 c11];

BT          = sum(bb_rep.*bcf(:,:),2);
CT          = sum(cc_rep.*ccf(:,:),2);
TG          = BT./(12-AT_rep) + CT;
F           = BT./(TG.*(1 - CT./TG).*(1 - CT./TG));

outvalues   =[AT_rep BT CT TG F];
end

function [nox xmf] = molepct_grd(wtm)
[nr nc] =   size(wtm);
% [n x]=MOLEPCT(X) Calculates mole percent oxides from wt %
% 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1
% Output: mole fractions of equivalent
mw=[60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];
mp=[];
xmf=[];

%Replicate the molecular weight row by the number of spatial nodes used
%(number of data rows)
mw_rep = repmat(mw,nr,1);

%Perform the calculation using matrix algebra
wtn = [wtm(:,1:10).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12)) wtm(:,11) 0.5.*wtm(:,12).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12))];
mp=wtn./mw_rep;
mpv= 100*(mp./sum(mp,2));
xmf=mpv;

[nr, nc]=size(xmf);
nox =nc;
end