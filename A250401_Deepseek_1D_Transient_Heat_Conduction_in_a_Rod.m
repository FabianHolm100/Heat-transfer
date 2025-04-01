% Heat Transfer Problem: 1D Transient Heat Conduction in a Rod

% Problem Statement:
% Consider a 1-meter long metal rod initially at a uniform temperature of 20°C. At time 
% t=0, the left end (x=0) is suddenly heated to 100°C, while the right end
% (x=1) is kept insulated (no heat flux). 

% The rod has:
% Thermal diffusivity (α) = 1×10^-5 m^2/s
% Discretize the rod into N = 10 nodes
% Simulate the temperature distribution over time using the explicit finite difference method
% Run the simulation for t = 1000s with a time step Δt=1s

% Tasks:

% Plot the temperature distribution along the rod at t = 100,500, and 1000s
% Ensure numerical stability by checking the Fourier number condition
% (Fo=αΔt/(Δx)^2 ≤ 0.5

% Theory:

% If 𝐹𝑜=𝛼Δ𝑡/(Δ𝑥)^ > 0.5 the explicit finite difference scheme for heat conduction becomes unstable, meaning that
% small numerical errors grow exponentially with each time step, leading to oscillations or divergence in the
% solution. This happens because the time step is too large relative to the spatial step causing heat to
% propagate too quickly and violating the fundamental assumptions of the numerical method.

% The explicit finite difference method for 1D heat conduction:
% 𝑇_new,𝑖 = 𝑇_old,𝑖 + Fo(T_old,𝑖+1-2*T_old,𝑖 + T_old,𝑖-1)

% Solution

alpha = 1*10^-5;    % Thermal diffusity [m^2/s]
L = 1;              % Lenght of rod [m]
N = 10;             % Number of rodes [n]
dx = L/(N-1);       % Spatial step: A grid with 𝑁 points have N−1 intervals (steps) between them. [m]
T_Left = 100;       % Left boundary temperature [°C]
T_in = 20;          % Initial temperature of the rod [°C]
dt = 1;             % Time step [s]
t_tot = 1000;       % Total simulation time [s]
nt = 1000/dt;       % Number of time steps [s]

% Stability check (Fourier number)
Fo = alpha*dt/(dx^2)
if Fo > 0.5 
     error('Unstable: Fourier number > 0.5. Reduce dt or increase dx.');
else
    fprintf('Stable: Fourier number = %.2f\n', Fo); % displayed with two decimal places and moves to a new line for readability.
end
 
% Initialize temperature array
T = T_in * ones(N, 1);
T(1) = T_Left; % Left boundary condition

% Store results for plotting
T_records = zeros(N, 3);
record_times = [100, 500, 1000];
record_steps = record_times / dt;

% Explicit finite difference method
for k = 1:nt
    T_old = T;
    
    % Update internal nodes
    for i = 2:N-1
        T(i) = T_old(i) + Fo * (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end
     % Right boundary (insulated: dT/dx = 0)
    T(N) = T_old(N) + Fo * (T_old(N-1) - T_old(N));
    
    % Record temperatures at specified times
    if ismember(k, record_steps)
        idx = find(k == record_steps);
        T_records(:, idx) = T;
    end
end

% Plotting
x = linspace(0, L, N);
figure;
plot(x, T_records(:,1), 'r', 'DisplayName', 't = 100 s');
hold on;
plot(x, T_records(:,2), 'g', 'DisplayName', 't = 500 s');
plot(x, T_records(:,3), 'b', 'DisplayName', 't = 1000 s');
xlabel('Position (m)');
ylabel('Temperature (°C)');
title('1D Transient Heat Conduction in a Rod');
legend;
grid on;