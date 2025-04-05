lx = 0.2; % [m] (x-direction, length of plate)

ly = 0.1; % [m]  (y-direction, width of plate)

T_top = 100; % [°C] Temp. at y = 0.1

T_bottom = 20; %  [°C]  Temp. at y = 0

k = 10; % [W/(m*K)]

% Boundary condition at left edge: Insulated (dT/dx = 0)

% Boundary condition at right edge: Insulated (dT/dx = 0)


% Discretize the domain: 

nx = 21; % count of discrete points used to represent the continuous domain along the x-axis

ny = 11; % count of discrete points used to represent the continuous domain along the y-axis

dx = lx/(nx-1); % (Number of spatial steps in x-direction)

dy = ly/(ny-1); % (Number of spatial steps in y-direction)



% Finite Difference approximation: 

% Start from then governing heat equation: ∂²T/∂x² + ∂²T/∂y² = 0

% We can approximate the second-order partial derivatives using central difference approximations:

% ∂²T/∂x² ≈ (T(i+1, j) - 2*T(i, j) + T(i-1, j)) / dx²

% ∂²T/∂y² ≈ (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / dy²

% Substituting these into Laplace's equation gives the finite difference equation for an interior % node (i, j):

% (T(i+1, j) - 2*T(i, j) + T(i-1, j)) / dx² + (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / dy² = 0



% Rearranging to solve for T(i, j):

% T(i, j) = (dy² * (T(i+1, j) + T(i-1, j)) + dx² * (T(i, j+1) + T(i, j-1))) / (2 * (dx² + dy²))

% Top Surface (j = ny): T(i, ny) = T_top for all i

% Bottom Surface (j = 1): T(i, 1) = T_bottom for all i

% T(0, j) = T(1, j) (Left Edge (i = 1, insulated))

% T(nx+1, j) = T(nx, j) (Right Edge (i = nx, insulated))



% Initialize temperature matrix

T = zeros(ny, nx); 



% Apply boundary conditions

T(1, :) = T_bottom;      % Bottom surface

T(ny, :) = T_top;       % Top surface



% Gauss-Seidel Method for a 2D heat transfer problem

max_iter = 1000;
tolerance = 1e-6;
error = 1;
iter = 0; 

while error > tolerance && iter < max_iter

    T_old = T;

    error_max = 0;

    % Iterate through interior nodes

    for j = 2:ny-1

        for i = 2:nx-1

            T(j, i) = (dy^2 * (T(j, i+1) + T(j, i-1)) + dx^2 * (T(j+1, i) + T(j-1, i))) / (2 * (dx^2 + dy^2));

            error_max = max(error_max, abs(T(j, i) - T_old(j, i)));

        end

    end

    % Apply insulated boundary conditions (after updating interior nodes)

    for j = 1:ny

        T(j, 1) = T(j, 2);      % Left insulated

        T(j, nx) = T(j, nx-1);   % Right insulated

    end

    error = error_max;

    iter = iter + 1;

end



fprintf('Solution converged in %d iterations with a maximum error of %e\n', iter, error);



% Visualization

[X, Y] = meshgrid(linspace(0, lx, nx), linspace(0, ly, ny));figure;

surf(X, Y, T);

xlabel('Length (m)');

ylabel('Width (m)');

zlabel('Temperature (°C)');

title('2D Steady-State Heat Conduction in a Rectangular Plate');

colorbar;

view(3);