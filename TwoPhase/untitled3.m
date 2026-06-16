phi=linspace(0,0.2,1000);
H=1./(1+exp(1e3*(phi-0.1)));


H2=phi-5e-2+sqrt((phi-5e-2).^2+eps);

k=800;
K=k;
H3=1/k*log(exp(k*phi)+exp(k*5e-2));
H4=-1/k*log(exp(-k*phi)+exp(-k*5e-2));

eps=1e-6;
 H5=   0.5*((phi-5e-2)-sqrt((phi-5e-2).^2+eps));

 H6=1./(1+exp(K*(-phi+3e-2)));
figure(9)
clf
plot(phi,H6)



    % SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
    % Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));

    %%
% Tl=Tl(336);
% T_1=T(336);
% Ts=Ts(336);
% cs_1=Cs(336);
% Cb=Cb(336);
cs_min=1e-2;
K=1e3;
    Cond4=((Tl-T_1)/(Tl-Ts))^n_order;
    % SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
    % Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;

    SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
    Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));
    solidus=(Constrain4-cs_1)/1e4

