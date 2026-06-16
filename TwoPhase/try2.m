eps2=1e-0;

Cb=0.3;
T_1=linspace(500,800,1000);
Ts=700;
solidus=Cb/2*(1-tanh((T_1-Ts)/eps2));

figure(10)
clf
plot(solidus,T_1)
xlim([0,1])

%%
n_order=20;
cs=linspace(0,1,1000);

T1=(1-cs).^n_order;
T2=1-cs.^(1/n_order);

figure(10)
clf; hold on;
plot(cs,T1)
plot(cs,T2)
xlim([0,1])