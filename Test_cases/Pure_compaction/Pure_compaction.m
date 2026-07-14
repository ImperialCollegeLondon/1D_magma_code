function [Results, Err]=Pure_compaction()

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'TwoPhase')));

% test Newton's method
Newton_finished=0;
try    
    Validation=2;
    CV_Twophase8_master
    data=readtable("output_1_CELLS.txt");
    x3=data.Depth_km_; y3=data.MeltFrac_;
    disp('Newton solver finished')
    Newton_finished=1;
catch
    disp('Newton method failed')
end

clearvars -except Newton_finished x3 y3;
% test nonlinear method
Nonlinear_finished=0;
try 
    Validation=1; 
    CV_Twophase8_master
    data=readtable("output_1_CELLS.txt");
    x2=data.Depth_km_; y2=data.MeltFrac_;
    disp('Nonlinear solver finished')
    Nonlinear_finished=1;
catch
    disp('Nonlinear method failed')
end

% load the expected McKenzie data
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
x1=x_data; y1=McKenzie_data(:,end);



% figure(9)
% clf; hold on
% plot(data.Depth_km_, data.MeltFrac_)
% plot(data2.Depth_km_, data2.MeltFrac_)
% plot(x_data, McKenzie_data(:,end))
% legend({'NL',"Newton",'MC'})


Dist=[-12 -11.5];
x_common = linspace(Dist(1)+1e-2, Dist(2)-1e-2, 1000);  % Avoiding exactly 0 and 1

mask1 = (x1 > Dist(1)) & (x1 < Dist(2));
x1_range = x1(mask1); y1_range = y1(mask1);
y1_interp = interp1(x1_range, y1_range, x_common, 'linear','extrap');

err1=0; err2=0;
if Nonlinear_finished==1
    mask2 = (x2 > Dist(1)) & (x2 < Dist(2));
    x2_range = x2(mask2);
    y2_range = y2(mask2);    
    y2_interp = interp1(x2_range, y2_range, x_common, 'linear');
    diff_12 = y1_interp - y2_interp;
    % norm_L2_12 = sqrt(trapz(x_common, diff_12.^2));
    err1 = sqrt(mean(diff_12.^2));
end

if Newton_finished==1
    mask3 = (x3 > Dist(1)) & (x3 < Dist(2));
    x3_range = x3(mask3);
    y3_range = y3(mask3);
    y3_interp = interp1(x3_range, y3_range, x_common, 'linear','extrap');
    diff_13 = y1_interp - y3_interp;
    % norm_L2_13 = sqrt(trapz(x_common, diff_13.^2));
    err2 = sqrt(mean(diff_13.^2));
end

Results=[Nonlinear_finished, Newton_finished];
Err=[err1, err2];

% subplot(2,1,1);
% plot(x1, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dataset 1');
% hold on;
% plot(x2, y2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Dataset 2');
% plot(x3, y3, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Dataset 3');
% xline([0, 1], 'k--', 'LineWidth', 1);
% xlabel('x');
% ylabel('y');
% title('All Datasets (Region of interest highlighted)');
% legend('Location', 'best');
% grid on;
% xlim([-0.5, 1.5]);
% 
% subplot(2,1,2);
% plot(x_common, y1_interp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dataset 1 (interpolated)');
% hold on;
% plot(x_common, y2_interp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Dataset 2 (interpolated)');
% plot(x_common, y3_interp, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Dataset 3 (interpolated)');
% xlabel('x');
% ylabel('y');
% title('Interpolated Data in Range 0 < x < 1');
% legend('Location', 'best');
% grid on;
% 
% % Plot differences
% figure('Position', [100, 100, 900, 400]);
% plot(x_common, diff_12, 'r-', 'LineWidth', 1.5, 'DisplayName', 'y1 - y2');
% hold on;
% plot(x_common, diff_13, 'b-', 'LineWidth', 1.5, 'DisplayName', 'y1 - y3');
% xlabel('x');
% ylabel('Difference');
% title('Differences between Datasets (0 < x < 1)');
% legend('Location', 'best');
% grid on;

end