
N=500;
Fe3ratio=linspace(0.05,0.5,N);

dFMQ=zeros(1,N);
fSO3=zeros(1,N);
ln_CS6_plus=zeros(1,N);
ln_SCSS=zeros(1,N);

for i=1:N
    results=calculate_S_redox(2.63, 5.88, 13.36, 50.93, 0.17, 10.81, 2.46, 0.25, 14.21,    Fe3ratio(i), 0, 1763, 2, 0.026, 41.8, 58);
    dFMQ(i)=results.DeltaQFM;
    fSO3(i)=results.fSO3;
    ln_CS6_plus(i)=results.ln_CS6_plus;
    ln_SCSS(i)=results.ln_SCSS;
end


%%
Font=20;
S6=exp(ln_CS6_plus).*fSO3;
figure(1)
clf;set(gcf,'color','w'); hold on
plot(dFMQ, exp(ln_SCSS),'LineWidth',2)
plot(dFMQ, S6,'LineWidth',2)
legend({'SCSS', 'S6+'},'fontsize',Font)
xlabel('dFMQ (-)','fontsize',Font)
ylabel('(ppm)', 'fontsize',Font)