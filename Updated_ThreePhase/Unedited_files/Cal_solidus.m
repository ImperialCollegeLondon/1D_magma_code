function [Ts,oversaturated, Saturation]=Cal_solidus(Vper,P,Data_point,Data_y)
%Vper: weight percent of volatile component, in 100%
%P: pressure in bar
addx=[634, 640, 648, 657, 670, 699, 730, 780, 833, 888, 953, 1020, 1095, 1160];
if P>50000
    error('Maximum allowed pressure is 50K bar');
end
if P>8000
    if Vper>13
        Saturation=13;
        oversaturated=1;
        Ts=addx(1);
    else
        z3=Vper;
        x_ind=13-floor(z3);
        y1=10000;
        y2=8000;
        y3=P;
        X1=addx(x_ind)*(z3-floor(z3))+addx(x_ind+1)*(floor(z3)-z3+1);
        X2=Data_point(1,x_ind)*(z3-floor(z3))+Data_point(1,x_ind+1)*(floor(z3)-z3+1);
        Ts=(X2*(y3-y1)+X1*(y2-y3))/(y2-y1);
    end
    return
end


y_ind=find(P>Data_y,1,'first')-1;
x1=Data_point(y_ind,y_ind);
y1=Data_y(y_ind);
z1=14-y_ind;

x2=Data_point(y_ind+1,y_ind+1);
y2=Data_y(y_ind+1);
z2=13-y_ind;

y3=P;
z3=Vper;

Saturation=0;

if (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)<0
    oversaturated=1;
    x3=(y3-y1)/(y2-y1)*(x2-x1)+x1;
    Saturation=(z2*(x3-x1)+z1*(x2-x3))/(x2-x1);
elseif z3>z2
    temp=[x1, Data_point(y_ind,y_ind+1), x2]';
    temp_A=[z1 y1 1;
        z2 y1 1;
        z2 y2 1];
    sol=temp_A\temp;
    x3=z3*sol(1)+y3*sol(2)+sol(3);
    oversaturated=0;
else
    x_ind=14-floor(Vper);
    X1=(z3-floor(z3))*Data_point(y_ind,x_ind-1)  +    (floor(z3)-z3+1)*Data_point(y_ind,x_ind);
    X2=(z3-floor(z3))*Data_point(y_ind+1,x_ind-1)+  (floor(z3)-z3+1)*Data_point(y_ind+1,x_ind);
    x3=(X2*(y3-y1)+X1*(y2-y3))/(y2-y1);
    oversaturated=0;
end
Ts=x3;