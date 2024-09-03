% ME 428: Final Project 
% ME 428 FINAL PROJECT: SOLUTION OF THE NON-LINEAR PENDULUM PROBLEM
% BY: ALEX PADRON (UIN: 651174065)
% Runge-Kutta 4th order method for solving a non-linear and Linear 2nd-order ODE
% Central Finite Difference scheme to calculate the first to third
% derivatives of the Linear case found using RK4

clear; % for clearing workspace
% Parameters obtained from [Table. 1]
g = 9.8;   % Gravity constant
L = 3.0;    % Length 
hValues = [1.0, 0.50, 0.10, 0.01, 0.001];  % Step sizes
duration = 10; 

% RK4 Loop for the 3 different cases 
for i = 1:length(hValues)
    
    h = hValues(i); 
    N = duration / h;  % Number of steps for the given duration and step size
    
    % Preallocate Matrices (kept giving me indexing error otherwise?)
    Y1 = zeros(1, N+1);
    Y2 = zeros(1, N+1);
    t = zeros(1, N+1);
    v = zeros(1, N+1);
    a = zeros(1, N+1);
    jerk = zeros(1, N+1);


    % Boundary Conditions:
    Y1(1) = pi/2;
    Y2(1) = 0;
    Y1Lin(1) = pi/2;
    Y2Lin(1) = 0;
    v(1) = 0;
    
    % Converted 2nd-Order ODE:
    dY1 = @(Y2) Y2; % angular position
    dY2 = @(Y1) -(g/L)*sin(Y1); % angular velocity
    % LINEAR FUNCTIONS USING SMALL-ANGLE APPROX:
    dY1Lin = @(Y2Lin) Y2Lin;
    dY2Lin = @(Y1Lin) -(g/L)*Y1Lin;

    
    % RK4 Method [K = linear, k = non-linear]
    for n = 1:N
        
        K11 = dY1Lin(Y2Lin(n));
        K12 = dY2(Y1Lin(n));
        
        k11 = dY1(Y2(n));
        k12 = dY2(Y1(n));
        
        K21 = dY1Lin(Y2Lin(n) + 0.5*h*K12);
        K22 = dY2Lin(Y1Lin(n) + 0.5*h*K11);
        
        k21 = dY1(Y2(n) + 0.5 * h * k12);
        k22 = dY2(Y1(n) + 0.5 * h * k11);
        
        K31 = dY1Lin(Y2Lin(n) + 0.5*h*K22);
        K32 = dY2Lin(Y1Lin(n) + 0.5*h*K21);
        
        k31 = dY1(Y2(n) + 0.5 * h * k22);
        k32 = dY2(Y1(n) + 0.5 * h * k21);
        
        K41 = dY1Lin(Y2Lin(n) + h * K32);
        K42 = dY2Lin(Y1Lin(n) + h * K31);
        
        k41 = dY1(Y2(n) + h * k32);
        k42 = dY2(Y1(n) + h * k31);
        
        % Linear Solution
        Y1Lin(n+1) = Y1Lin(n) + (h / 6.0) * (K11+ (2*K21) + ( 2*K31) + K41);
        Y2Lin(n+1) = Y2Lin(n) + (h / 6.0) * (K12 + (2*K22) + (2*K32) + K42);
        
        % Non-Linear solution
        Y1(n+1) = Y1(n) + (h / 6.0) * (k11 + (2*k21) + (2*k31) + k41);
        Y2(n+1) = Y2(n) + (h / 6.0) * (k12 + (2*k22) + (2*k32) + k42);
        t(n+1) = t(n) + h;
      
    end % end of RK4 Calculations (only variation is the step-size)
    
    % Finite Difference Loops:
    for n = 2:N
       
        % Central Finite Difference for first derivative (velocity)
        
        v(n) = (Y1Lin(n+1) - Y1Lin(n-1)) / (2*h);
    
    end
    
    
    for n = 2:N
        
        % Central Finite Difference for second derivative (acceleration)
        a(n) = (Y1Lin(n+1) - 2*Y1Lin(n) + Y1Lin(n-1)) / (h^2);
    
    end
    
    
    for n = 3:N-1
        
        % Central Finite Difference for third derivative (Jerk)
        jerk(n) = (Y1Lin(n+2) - 2*Y1Lin(n+1) + 2*Y1Lin(n-1) - Y1Lin(n-2)) / (2*h^3);
    
    end
%   -------------------START OF PLOTTING SECTION-------------------------------

    % Plot of Linear case for Y1
    hold on;
    figure(1)
    grid on;    
    axis([0.0, 10.0, min(Y1Lin), max(Y1Lin)])
    plot(t,Y1Lin, 'DisplayName',['h = ', num2str(h)]);
    xlabel('Time [s]');
    ylabel('Y1(Angular Position) [rad]');
    title('Linear Solution to Angular Position'); 
    hold off;
    
    % Plot results of derivatives overlapped with different step-sizes of the LINEAR case:
    
    % Velocity
    hold on;
    figure(2)
    grid on;    
    axis([0.0, 10.0, -3, 3])
    plot(t,v, 'DisplayName',['h = ', num2str(h)]);
    xlabel('Time [s]');
    ylabel('Velocity [rads/s]');
    title('First Derivative of Angular Position (velocity)'); 
    hold off;
    
    %Acceleration
    hold on;
    figure(3)
    grid on;    
    axis([0.0, 10.0, -8, 8])
    plot(t,a, 'DisplayName',['h = ', num2str(h)]);
    xlabel('Time [s]');
    ylabel('Accerlation [rads/s^2]');
    title('Second Derivative of Angular Position (acceleration)'); 
    hold off;
    
    % Jerk
    hold on;
    figure(4)
    grid on;    
    axis([0.0, 10.0, -8, 8])
    plot(t,jerk, 'DisplayName',['h = ', num2str(h)]);
    xlabel('Time [s]');
    ylabel('Jerk [rads/s^3]');
    title('Third Derivative of Angular Position (jerk)'); 
    hold off;


    % Non-Linear case: Angular Position
    hold on;
    figure(5);
    grid on;
    plot(t, Y1, 'DisplayName',['h = ', num2str(h)]);
    % sets axis limits:
    axis([0.0, 10.0, -3, 3]) % x-axis from 0 to 10, y-axis from -+0.50 offset of min and max values of function
    xlabel('Time [s]');
    ylabel('Y1(Angular Position) [rad]');
    title('Non-Linear Solution to Angular Position'); % for use in distinguishing different graphs
    hold off;
    
    % Non-Linear case: Angular Position
    hold on;
    figure(6);
    grid on;
    plot(t, Y2, 'DisplayName',['h = ', num2str(h)]);
    axis([0.0, 10.0, -3, 3]) % x-axis from 0 to 10, y-axis from -+0.50 offset of min and max values of function
    xlabel('Time [s]');
    ylabel('Y2 (Non-Linear Angular Velocity) [rad/s]');
    title('Non-Linear Solution to Angular Velocity'); % for use in distinguishing different graphs
    hold off;
    
end 
