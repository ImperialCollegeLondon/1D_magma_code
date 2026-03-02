

function [Jac_mom, rhs_mom, Jac_con, rhs_con, Jac_ct, rhs_ct, Jac_ent, rhs_ent, Jac_solidus, rhs_solidus, Jac_liquidus, rhs_liquidus]=Two_phase_Newton_initiator(Non_dimention)



syms phi0 phi1 umi_0 umi_1 umi_2  dzi_0 dzi_1

ufi = sym('ufi');
syms cs cl
syms cs0 cl0 cs2 cl2 phi2

%parameters
% syms mu_all c
% c=c0*phi^(2-B); c0=A/d^2
%mu_f=10^((muf2-muf1)*cl+muf1)
syms muf1 muf2 rhom_1 rhom_2 rhof_1  rhof_2 g
syms mum_0 c0 

Variables=[umi_0 umi_1 umi_2 ufi phi0 phi1 cs0 cs cl0 cl];
parameters=[muf1 muf2 mum_0 c0 dzi_0 dzi_1 rhom_1 rhom_2 rhof_1  rhof_2 g];

drhog=((cs0*rhom_1+(1-cs0)*rhom_2+cs*rhom_1+(1-cs)*rhom_2)/2-(cl0*rhof_1+(1-cl0)*rhof_2+cl*rhof_1+(1-cl)*rhof_2)/2)*g;

eps=1e-2;
momentum=(-((umi_2-umi_1)/dzi_1*(1-phi1)*mum_0*(4/3+(phi1+eps)^-0.5)-(umi_1-umi_0)/dzi_0*(1-phi0)*mum_0*(4/3+(phi0+eps)^-0.5))*2/(dzi_0+dzi_1) *2*phi0*phi1/(phi0+phi1+1e-12)...
        +phi0*phi1/(phi0+phi1+1e-12)*(2-phi1-phi0+1e-12)*drhog-c0*((phi1+phi0+eps)/2)^-1*10^((muf2-muf1)*cl+muf1)*(ufi-umi_1))/mum_0;



Jac_mom=jacobian(momentum, Variables);
Jac_mom= matlabFunction(Jac_mom, 'Vars', [Variables parameters]);
rhs_mom= matlabFunction(momentum, 'Vars', [Variables parameters]);

continuity=umi_1*(1-phi1/2-phi0/2)+ufi*(phi1/2+phi0/2);
Variables=[umi_1, ufi, phi0, phi1];
Jac_con=jacobian(continuity, Variables);
Jac_con= matlabFunction(Jac_con, 'Vars', [umi_1 ufi phi0 phi1]);
rhs_con= matlabFunction(continuity, 'Vars', [umi_1 ufi phi0 phi1]);


syms dt
syms OLD_com
syms ufi2
% syms phis0 phis1 phil0 phil1

% Variables=[umi_1 umi_2 ufi ufi2 phi1 phis0 phis1 phil0 phil1 cs0 cs cl0 cl];
% trans_comp=(phi1*cl+(1-phi1)*cs)/dt-OLD_com/dt+(ufi*phil0*cl0-ufi2*phil1*cl+umi_1*(1-phis0)*cs0-umi_2*(1-phis1)*cs)/dzi_1;

Variables=[umi_1 umi_2 ufi ufi2 phi0 phi1 phi2 cs0 cs cs2 cl0 cl cl2];
trans_comp=(phi1*cl+(1-phi1)*cs)-OLD_com-(ufi*(phi0*cl0+phi1*cl)/2-ufi2*(phi1*cl+phi2*cl2)/2+umi_1*((1-phi0)/2*cs0+(1-phi1)/2*cs)-umi_2*((1-phi1)/2*cs+(1-phi2)/2*cs2))/dzi_1*dt;
Jac_ct=jacobian(trans_comp, Variables);
Jac_ct= matlabFunction(Jac_ct, 'Vars', [Variables dt dzi_1 OLD_com]);
rhs_ct= matlabFunction(trans_comp, 'Vars', [Variables dt dzi_1 OLD_com]);



syms T0 T1 T2 dzi_2
syms cp Lf kt
syms OLD_ent
trans_enthalpy=((cp*T1+Lf*phi1)-OLD_ent-Lf*((phi0+phi1)/2*ufi-(phi1+phi2)/2*ufi2)/dzi_1*dt-kt*2/dzi_1*((T2-T1)/(dzi_2+dzi_1)-(T1-T0)/(dzi_1+dzi_0))*dt)/550e3;
Variables=[ufi ufi2 phi0 phi1 phi2 T0 T1 T2];
Jac_ent=jacobian(trans_enthalpy, Variables);
Jac_ent=matlabFunction(Jac_ent,'Vars',[Variables dzi_0 dzi_1 dzi_2 dt cp Lf kt OLD_ent]);
rhs_ent= matlabFunction(trans_enthalpy, 'Vars', [Variables dzi_0 dzi_1 dzi_2 dt cp Lf kt OLD_ent]);




% solidus and liquidus
eps=1e-12;

syms a1 b1 c1 n Ts Tl
Cb=phi1*cl+(1-phi1)*cs;

if Non_dimention==0
    MINT=((T1+1350)-sqrt((T1-1350)^2+eps))/2;
    Cond5=(-b1-sqrt(b1^2-4*a1*(c1-MINT)))/2/a1;
    SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps))/2;
    Constrain5=(Cb+SMIN+sqrt((SMIN-Cb)^2+eps))/2;
else
    MINT=((T1+1.15)-sqrt((T1-1.15)^2+eps))/2;
    Cond5=(-b1-sqrt(b1^2-4*a1*(c1-MINT)))/2/a1;
    SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps))/2;
    Constrain5=(Cb+SMIN+sqrt((SMIN-Cb)^2+eps))/2;
end

liquidus=(Constrain5-cl)/1e2;
Variables=[phi1 T1 cs cl];
Jac_liquidus=jacobian(liquidus, Variables);
Jac_liquidus=matlabFunction(Jac_liquidus, 'Vars',[Variables a1 b1 c1]);
rhs_liquidus=matlabFunction(liquidus,'Vars',[Variables a1 b1 c1 ]);

if Non_dimention==0
    Cond4=((Tl-T1)/(Tl-Ts))^n;
    SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
    Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;
else
    Cond4=(1-T1)^n;
    SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
    Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;
end

solidus=(Constrain4-cs)/1e2;
Variables=[phi1 T1 cs cl];
Jac_solidus=jacobian(solidus, Variables);
Jac_solidus=matlabFunction(Jac_solidus,'Vars',[Variables n Ts Tl]);
rhs_solidus=matlabFunction(solidus,'Vars',[Variables n Ts Tl ]);

end