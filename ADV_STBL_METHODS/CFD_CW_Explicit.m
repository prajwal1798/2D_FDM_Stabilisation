clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing (Input Parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax  = 0.05;              % Wavelength of Initial Velocity Profile [m]
nx    = 501;               % Number of Grid Points per wave
dx    = xmax/(nx-1);       % Grid Size [m]
t     = 0.5;               % Total simulation time [s]
dt    = 1e-3;              % Time step size [s]
nt    = t/dt;              % Number of time steps

x     = linspace(0,1,nx);  % Domain definition [m]
c = 0.1;                   % Velocity of the wave [m/s] 

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
u(1) = 0;                % Initial and Boundary Condition (Set the value at x=0 to 0)
u(x >= 0.05) = 0;        % Initial Condition (Set the values beyond x >= 0.05 to 0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing (Create the initial plot) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2, 1, 1);
plot(x, u, 'LineWidth', 1.5, Color = 'r');
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
% Solving (Explicit Time marching and Space marching)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform the time-stepping loop

tic
for n = 1:nt
    un = u;

    % Perform the space loop
    for i = 2:nx
        u(i) = un(i) - CFL * (un(i) - un(i-1));
    end 

    % Post-processing (Plot the convected velocity profile)
    subplot(2, 1, 2);    
    plot(x, u, 'LineWidth', 1.5, color = 'b');
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
    pause(0.00001);    
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
