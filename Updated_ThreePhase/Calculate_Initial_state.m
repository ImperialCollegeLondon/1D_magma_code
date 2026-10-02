function x=Calculate_Initial_state(PD_range, rho_mean, SiO2_crust, H2O_crust,Precision)

% Calculate the initial amounts of the evolved component (M), refractory
% component (N) and volatile (V) required to satisfy the specified crustal
% composition and total amount constraint (based on density).

% x = [M; N; V]
% M = evolved component
% N = refractory component
% V = volatile

% PD_range = [minimum SiO2, maximum SiO2] within the model
% rho_mean = target value for the total amount of M+N+V - WORRIED ABOUT
% THIS 
% SiO2_Crust = required SiO2 content of the crust
% H2O = required H2O content of the crust
% Precision = convergence criterion for the iterative solution

% Set the initial error to 1 so that the iterative loop starts.
err=1;

% Initial guess for the amounts of the three components:
x=[1500, 1500, 20]';  %(M,N,V)

% Extract the current values of M, N and V from the solution vector.
M=x(1); N=x(2); V=x(3);

%% Initialise the 3x3 Jacobian matrix.
% The Jacobian contains the derivatives of the three equations being solved
% wth respect to the three unknowns [M, N, V].
J=zeros(3);

% Equation 1 - SiO2 composition
% This  row contains the derivatives of the SiO2 equation with respect to
% M, N and V
J(1,:)=[PD_range(2)-SiO2_crust, PD_range(1)-SiO2_crust, 0];  %(PD_range(2)*M+ PD_range(1)*N)/(M+N)=SiO2_crust

% Equation 2 - H2O composition
% This row contains the derivative of the H2O equation with respect to
% M,N,V
J(2,:)=[-H2O_crust, -H2O_crust, 100-H2O_crust]; %100*V/(M+N+V)=H2O_crust


%Equation 3: total amount constraint
%This row contains the total amount of the three components is constrained
%to equal rho_mean
J(3,:)=[1,1,1]; %M + N + V = rho_mean



% Continue iterating until the correction to M, N and V is smaller than the
% specified precision. 
while err>Precision
    
    % Calculate the residual of each of the three equations.
    % The residual indicates how far the current values of M, N and V are
    % from satisfying each equation.
    rhs=[(PD_range(2)-SiO2_crust)*M+ (PD_range(1)-SiO2_crust)*N;...
        (100-H2O_crust)*V-H2O_crust*(M+N);...
        M+N+V-rho_mean];
    
    % Solve the linear system to calculate the correction required for M, N
    % and V. 
    dx=J\rhs;

    % Update the solution by applying the correction.
    x=x-dx;

    % Find the largest absolute correction to M, N or V. 
    % The iteration stops when this is smaller than Precision.
    err=max(abs(dx(:)));

    % Update M, N and V using the new solution.
    M=x(1); N=x(2); V=x(3);
end
       
end