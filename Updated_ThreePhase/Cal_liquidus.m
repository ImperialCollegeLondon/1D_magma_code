function [Ts,oversaturated, Saturation]=Cal_liquidus(Vper,P,Data_point,Data_y)
%Vper: weight percent of volatile component, in 100%
%P: pressure in bar


z_val=[20,10,5,2,-1e-8];

Saturation=0;

y_ind=max(find(P>Data_y,1,'first')-1,1);
x_ind=max(find(Vper>z_val,1,'first')-1,1);

y1=Data_y(y_ind);
y2=Data_y(y_ind+1);

z3=Vper;
y3=P;
% oversaturated=0;
if Data_point(y_ind,x_ind)<800
    x_ind2=find(Data_point(y_ind,:)>800,1,'first');
    x1=Data_point(y_ind,x_ind2);
    z1=z_val(x_ind2);
    if Data_point(y_ind+1,x_ind2)<800
        x2=Data_point(y_ind+1,x_ind2+1);
        z2=z_val(x_ind2+1);
    else
        x2=Data_point(y_ind+1,x_ind2);
        z2=z_val(x_ind2);
    end
    oversaturated=1;
    x3=(y3-y1)/(y2-y1)*(x2-x1)+x1;
    Saturation=(z2*(x3-x1)+z1*(x2-x3))/(x2-x1);   
elseif Data_point(y_ind+1,x_ind)<800
    x1=Data_point(y_ind,x_ind);
    x2=Data_point(y_ind+1,x_ind+1);
    z1=z_val(x_ind);
    z2=z_val(x_ind+1);
    if (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)<0
        oversaturated=1;
        x3=(y3-y1)/(y2-y1)*(x2-x1)+x1;
        Saturation=(z2*(x3-x1)+z1*(x2-x3))/(x2-x1);
    else
        temp=[x1, Data_point(y_ind,x_ind+1), x2]';
        temp_A=[z1 y1 1;
            z2 y1 1;
            z2 y2 1];
        sol=temp_A\temp;
        x3=z3*sol(1)+y3*sol(2)+sol(3);
        oversaturated=0;
    end
else
    X1=((z3-z_val(x_ind+1))*Data_point(y_ind,x_ind)  +    (z_val(x_ind)-z3)*Data_point(y_ind,x_ind+1))/(z_val(x_ind)-z_val(x_ind+1));
    X2=((z3-z_val(x_ind+1))*Data_point(y_ind+1,x_ind)+  (z_val(x_ind)-z3)*Data_point(y_ind+1,x_ind+1))/(z_val(x_ind)-z_val(x_ind+1));
    x3=(X2*(y3-y1)+X1*(y2-y3))/(y2-y1);
    oversaturated=0;
end

Ts=x3;