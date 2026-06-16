eps=1e-6;
phi=linspace(0,0.2,100);
T_1=linspace(600,1600,1000);

H=1./(1+exp(1e3*(phi-0.1)));


H2=(phi+5e-2+sqrt((phi-5e-2).^2+eps))/2;

k=100;
K=k;
H3=1/k*log(exp(k*phi)+exp(k*5e-2));
H4=1./(1+exp(K*(-phi+3e-2)));

MINT=-1/K*log(exp(-K*T_1/1000)+exp(-K*(C1+50)/1000))*1000;

Cond5=linspace(0.8,1.2,1000);
SMIN=-1/K*log(exp(-K*Cond5)+exp(-K*1));


eps=1e-6;
 H5=   0.5*((phi-5e-2)-sqrt((phi-5e-2).^2+eps));

H6=1./(1+exp(K*(-phi+3e-2)));
figure(9)
clf
% plot(T_1,MINT)
plot(H2)



    % SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
    % Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));

    %%
% Tl=Tl(336);
% T_1=T(336);
% Ts=Ts(336);
% cs_1=Cs(336);
% Cb=Cb(99);


K=2e3;
% Cb_local=phi_1(670)* cl_1(670)+(1- phi_1(670))*cs_1(670);
cs_min=0;
    % Cond4=((Tl0(1)-T_1(670))/(Tl0(1)-Ts0(1)))^n_order;
    % % SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
    % % Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;
    % 
    % SMIN1=-1/K*log(exp(-K*Cb_local)+exp(-K*Cond4));
    % Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));
    % solidus=(Constrain4-cs_1(670))/1e4
    a0=760;
    b0=1160;
    AAA=(exp(cl_1.*phi_1.*-2.0e+3+cs_1.*(phi_1-1.0).*2.0e+3).*1.0./(exp(((T_1-b0)./(a0-b0)).^n_order.*-2.0e+3)+exp(cl_1.*phi_1.*-2.0e+3+cs_1.*(phi_1-1.0).*2.0e+3)).^2.*(cl_1.*2.0e+3-cs_1.*2.0e+3))./(2.0e+7./(exp(((T_1-b0)./(a0-b0)).^n_order.*-2.0e+3)+exp(cl_1.*phi_1.*-2.0e+3+cs_1.*(phi_1-1.0).*2.0e+3))+2.0e+7);
%%
figure(9)
clf
% semilogy(abs(RHS));
plot((RHS));hold on;
YLIM=get(gca,'YLim');
ylim(YLIM);

plot([N+1, N+1], YLIM,'--','color','red')
plot([(N+1)*2, (N+1)*2], YLIM,'--','color','red')
plot([(N+1)*2+N (N+1)*2+N], YLIM,'--','color','red')
text((N+1)*2+N*.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'Major','fontsize',13)
plot([(N+1)*2+N*2 (N+1)*2+N*2], YLIM,'--','color','red')
text((N+1)*2+N*1.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'Ent','fontsize',13)

plot([(N+1)*2+N*3 (N+1)*2+N*3], YLIM,'--','color','red')
text((N+1)*2+N*2.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'Solidus','fontsize',13)
plot([(N+1)*2+N*4 (N+1)*2+N*4], YLIM,'--','color','red')
text((N+1)*2+N*3.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'Liquidus','fontsize',13)
plot([(N+1)*2+N*5 (N+1)*2+N*5], YLIM,'--','color','red')
text((N+1)*2+N*4.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'Vol trans','fontsize',13)
plot([(N+1)*2+N*6 (N+1)*2+N*6], YLIM,'--','color','red')
text((N+1)*2+N*5.3,(YLIM(2)-YLIM(1))*0.02+YLIM(1),'SSat','fontsize',13)

%%
Cb=0.99;
T_1=linspace(600,1600,1e3);
Tl=1160;
Ts=760;

C1=1160;
A1=-99.6;
B1=-300.4;
K=100;
n_order=15;
Cond4=((Tl-T_1)/(Tl-Ts)).^n_order;

cs_min=-2e-3;
SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
solidus=1/K*log(exp(K*SMIN1)+exp(K*cs_min));


eps=1e-3;
% MINT=-1/K*log(exp(-K*T_1/scaling)+exp(-K*(C1+50)/scaling))*scaling;
% MINT=((T_1+C1+50)-sqrt((T_1-C1-50).^2+eps))/2;
% Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
% SMIN=-1/K*log(exp(-K*Cond5)+exp(-K*1));
% liquidus=1/K*log(exp(K*Cb)+exp(K*SMIN));

MINT=((T_1+C1+50)-sqrt((T_1-C1-50).^2+eps))/2;
Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
SMIN=(Cond5+1-sqrt((Cond5-1).^2+eps))/2;
liquidus=(Cb+SMIN+sqrt((SMIN-Cb).^2+eps))/2;

figure(10)
clf; 

plot(solidus,T_1);
hold on
plot(liquidus,T_1);
% xlim([-5e-3,2e-2])
% xlim([0.98,1+1e-3])




%%
eps=1e-3;
Cond5=linspace(0.8,1.2,1e3);
H1=(Cond5+1-sqrt((Cond5-1).^2+eps))/2;
figure(9)
clf
plot(Cond5,H1)

%%
S_2=2e-10;
S_1=1e-10;

test=1e-10% 1e-10;

eps=1e-20;


K=1e8;
0.5*(test+sqrt(test^2+eps))

1/K*log(exp(K*test)+1);

%%
T = [610 612 614 616 618 620 622 624 626 628 ...
     630 632 634 636 638 640 642 644 646 648 ...
     650 655 660 665 670 675 680 690 700 715 ...
     730 750 770 790 810 830 850 870 890 905 ...
     920 930 940 945 948 950 952 953 954 955];

P = [9500 9000 8600 8200 7800 7400 7000 6600 6200 5900 ...
     5600 5300 5000 4700 4400 4100 3800 3500 3200 3000 ...
     2800 2500 2200 1900 1600 1300 1100 900 750 600 ...
     500 420 360 310 270 240 220 205 195 190 ...
     185 180 170 160 140 120 100 80 60 40];


dist=P*1e5/2800/9.81/1e3;

Geo=dist*25+20;

Font_size=15;
figure(9)
clf
set(gcf,'color','w')

plot(dist,T)
hold on
box off
plot(dist,Geo)
legend(gca,{'Saturated Solidus', 'Geotherm'},'fontsize', Font_size)
set(gca,'tickdir','out') ; set(gca,'fontsize', Font_size);   
xlabel('Depth (km)', 'FontSize', Font_size)
ylabel('Temperature (^\circ C)', 'FontSize', Font_size)

%%
figure(10)
clf

% plot(nodez,X(N+2:2*N+2));
% hold on
% plot(nodez,X(1:N+1));

subplot(2,1,1)
plot(X(N+2:2*N+2));
hold on
plot(X(1:N+1));


subplot(2,1,2)
plot(X((1:N)  +(N+1)*2))

%%
i=737;
Cb=phi_1.*cl_1+(1-phi_1).*cs_1;
    eps=1e-8;
    Cond4=((Tl0(1)-T_1(i))/(Tl0(1)-Ts0(1)))^n_order;
    SMIN1=(Cb(i)+Cond4-sqrt((Cb(i)-Cond4)^2+eps))/2;
    Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;

    % SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
    % Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));
    solidus=(Constrain4-cs_1(i))/1e4