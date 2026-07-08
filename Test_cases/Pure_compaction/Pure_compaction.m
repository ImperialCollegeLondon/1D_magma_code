data=readtable("output_1_CELLS_NL.txt");
data2=readtable("output_1_CELLS.txt");


McKenzie_data=readtable('MacKenzie.txt');
McKenzie_data=table2array(McKenzie_data);
g=9.81;
drho=500;
mu_f=1e3;
mu_m_mix=1e14;
d=2.75e-3;      %grain size
n=3;
b=125;

phi0=0.5;

K0=phi0^n*d^2/b;
um_0=K0/mu_f*drho*g*(1-phi0); %
sigma=(K0*mu_m_mix/mu_f)^0.5;   %Compaction length
tao0=sigma/um_0;    %time scaleA
McKenzie_data=McKenzie_data*phi0;
x_data=linspace(0,500,501)*sigma;
x_data=x_data/1000-12;




figure(9)
clf; hold on
plot(data.Depth_km_, data.MeltFrac_)
plot(data2.Depth_km_, data2.MeltFrac_)

plot(x_data, McKenzie_data(:,end))

legend({'NL',"Newton",'MC'})


x1=x_data; y1=McKenzie_data(:,end);
x2=data.Depth_km_; y2=data.MeltFrac_;
x3=data2.Depth_km_; y3=data2.MeltFrac_;

Dist=[-12 -11.5];
mask1 = (x1 > Dist(1)) & (x1 < Dist(2));
mask2 = (x2 > Dist(1)) & (x2 < Dist(2));
mask3 = (x3 > Dist(1)) & (x3 < Dist(2));

x1_range = x1(mask1);
y1_range = y1(mask1);
x2_range = x2(mask2);
y2_range = y2(mask2);
x3_range = x3(mask3);
y3_range = y3(mask3);

x_common = linspace(Dist(1)+1e-2, Dist(2)-1e-2, 1000);  % Avoiding exactly 0 and 1
y1_interp = interp1(x1_range, y1_range, x_common, 'linear','extrap');
y2_interp = interp1(x2_range, y2_range, x_common, 'linear');
y3_interp = interp1(x3_range, y3_range, x_common, 'linear','extrap');

diff_12 = y1_interp - y2_interp;
norm_L2_12 = sqrt(trapz(x_common, diff_12.^2));

diff_13 = y1_interp - y3_interp;
norm_L2_13 = sqrt(trapz(x_common, diff_13.^2));

% L-infinity norm (maximum absolute difference)
norm_Linf_12 = max(abs(diff_12));
norm_Linf_13 = max(abs(diff_13));

rms_12 = sqrt(mean(diff_12.^2));
rms_13 = sqrt(mean(diff_13.^2));



subplot(2,1,1);
plot(x1, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dataset 1');
hold on;
plot(x2, y2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Dataset 2');
plot(x3, y3, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Dataset 3');
xline([0, 1], 'k--', 'LineWidth', 1);
xlabel('x');
ylabel('y');
title('All Datasets (Region of interest highlighted)');
legend('Location', 'best');
grid on;
xlim([-0.5, 1.5]);

subplot(2,1,2);
plot(x_common, y1_interp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dataset 1 (interpolated)');
hold on;
plot(x_common, y2_interp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Dataset 2 (interpolated)');
plot(x_common, y3_interp, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Dataset 3 (interpolated)');
xlabel('x');
ylabel('y');
title('Interpolated Data in Range 0 < x < 1');
legend('Location', 'best');
grid on;

% Plot differences
figure('Position', [100, 100, 900, 400]);
plot(x_common, diff_12, 'r-', 'LineWidth', 1.5, 'DisplayName', 'y1 - y2');
hold on;
plot(x_common, diff_13, 'b-', 'LineWidth', 1.5, 'DisplayName', 'y1 - y3');
xlabel('x');
ylabel('Difference');
title('Differences between Datasets (0 < x < 1)');
legend('Location', 'best');
grid on;