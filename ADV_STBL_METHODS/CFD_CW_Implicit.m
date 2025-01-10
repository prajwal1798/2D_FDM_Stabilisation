clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax  = 0.05;              % Wavelength of Initial Velocity Profile  
nx    = 501;               % Number of Grid Points per wave 
dx    = xmax/(nx-1);       % Grid Size
t     = 1;              % Total simulation time
dt    = 1e-3;              % Time step size
nt    = t/dt;              % Number of time steps

x     = linspace(0,1,nx);  % Domain definition
c     = 0.1;               % Velocity of the wave 

fprintf('Grid Size: %.4fm\n', dx);
fprintf('Time Stepsize: %.3fs\n', dt);
fprintf('Wave Velocity: %.1fm/s\n', c);
fprintf('Number of Time Steps: %d\n', nt);

CFL = (c * dt) / dx;
fprintf('CFL Number: %.4f\n', CFL);

% Time instances to capture snapshots
snap_times = [0.01 , 0.1 , 0.2 , 0.3 , 0.4 , 0.45];
snapshots = cell ( size ( snap_times ) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising the Velocity Profile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = sin((pi/xmax) * x);  % Initialize u with a sine function
u(1) = 0;                % Set the value at x=0 to 0
u(x >= 0.05) = 0;        % Set the values beyond x >= 0.05 to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the initial plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2, 1, 1);
plot(x, u, 'LineWidth', 1.5, 'Color', 'r');
xlim([-0.05, 1.1]);
ylim([-0.05, 1.2]);
xlabel('X in m');
ylabel('Velocity in m/s');
title('Initial Velocity Profile');
grid on;
grid minor;
legend('Initial Velocity Profile', 'Location', 'northeast');
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit Time marching and Space marching 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the coefficient matrix A
A = diag(ones(1,nx)*(1+CFL)) + diag(-CFL*ones(1,nx-1),-1);

% Perform the time-stepping loop
tic
for n = 1:nt
    un = u;    
    u = linsolve(A, un'); % Solve using linsolve function
    
    % Reshape u to a column vector
    u = u';
    
    % Plot the convected velocity profile
    subplot(2, 1, 2);    
    plot(x, u, 'LineWidth', 1.5, 'Color', 'b');
    xlim([-0.05, 1.1]);
    ylim([-0.05, 1.2]);
    xlabel('X [m]');
    ylabel('Velocity [m/s]');
    grid on;
    grid minor;
    legend('Convected Velocity Profile', 'Location', 'northeast');
    run_time = toc;

    title_text1 = sprintf('Computation Time = %.2f s, CFL Number = %.4f', run_time, CFL);
    title({'Convected Velocity Profile';title_text1});
    pause(dt);  

    if (sum(u)<=0.001)
        break;
    end
end

% Post - processing
figure ;
for j = 1: numel (snap_times)
    plot (x , snapshots { j } , 'LineWidth ', 1.5) ;
    hold on ;
end
xlim ([ -0.05 , 1.1]) ;
ylim ([ -0.05 , 1.2]) ;
xlabel ('X [m]') ;
ylabel ('Molar Concentration [ mol /m^3] ') ;
grid on ;
grid minor ;
title ('Convected concentration Progression - Implicit Scheme ');
legend ({ 't = 0.01 s','t = 0.1 s', 't = 0.2 s', 't = 0.3 s', 't= 0.4 s', 't = 0.45 s'} , 'Location ', 'northeast ') ;