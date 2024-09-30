%% Introduction
disp('Group No.: 04');
disp('Course Code: EECE 306');
disp('Course: EECE 20, Sec B');



%% Bus Data
bus_data = [1   1.06   0    0      0     0   0; 
            2   1      0    20     20    40  30; 
            3   1      0   -45    -15    0   0; 
            4   1      0   -40    -5     0   0; 
            5   1      0   -60    -10    0   0];
disp(bus_data);

%% Line Data
transmission_line_data = [1  2  0.02  0.06  0.03; 
                          1  3  0.08  0.24  0.025; 
                          2  3  0.06  0.25  0.02; 
                          2  4  0.06  0.18  0.02; 
                          2  5  0.04  0.12  0.015; 
                          3  4  0.01  0.03  0.01; 
                          4  5  0.08  0.24  0.025];
disp(transmission_line_data);

%% Base Voltage
V_base = [10.6; 10; 10; 10; 10];
disp('Base Voltages (kV):');
disp(V_base);

%% Y-bus Matrix Calculation
Y_bus = zeros(5);
for i = 1:size(transmission_line_data, 1)
    f = transmission_line_data(i, 1); 
    t = transmission_line_data(i, 2);
    Z = transmission_line_data(i, 3) + 1i * transmission_line_data(i, 4);
    Y = 1 / Z;
    Y_shunt = 1i * transmission_line_data(i, 5);
    Y_bus(f, f) = Y_bus(f, f) + Y + Y_shunt;
    Y_bus(t, t) = Y_bus(t, t) + Y + Y_shunt;
    Y_bus(f, t) = Y_bus(f, t) - Y;
    Y_bus(t, f) = Y_bus(t, f) - Y;
end
disp('Y-bus Matrix:');
disp(Y_bus);

%% Gauss-Seidel Iteration using while loop
S_base = 100;  % Base power (MVA)
V = bus_data(:, 2);  % Initial guess for voltages (Voltage magnitude in per unit)
P = bus_data(:, 4);  % Active power demand (MW)
Q = bus_data(:, 5);  % Reactive power demand (MVAR)
tol = 1e-5;  % Tolerance for convergence
max_iter = 100;  % Maximum number of iterations
iteration = 0;  % Initialize iteration count
error = inf;  % Initial error value

while error > tol && iteration < max_iter
    V_prev = V;  % Store previous voltage values for convergence check
    iteration = iteration + 1;  % Increment iteration count
    
    for i = 2:5  % Loop over PQ buses (start from bus 2)
        sum_YV = 0;  % Initialize the sum of Y(i,j)*V(j)
        for j = 1:5
            if i ~= j
                sum_YV = sum_YV + Y_bus(i, j) * V(j);  % Accumulate for all 
                                                       %buses except the current bus i
            end
        end
        
        % Update the voltage at bus i using the Gauss-Seidel formula
        V(i) = (P(i) - 1i * Q(i)) / S_base / conj(V(i)) - sum_YV;
        V(i) = V(i) / Y_bus(i, i);  % Normalize by the diagonal element of the Y-bus matrix
        
        % Update with the new voltage magnitude
        V(i) = abs(V(i)) * exp(1i * angle(V(i)));  % Retain voltage magnitude and phase
    end
    
    % Calculate the error as the maximum voltage change
    error = max(abs(V - V_prev));
    
    % Check for convergence after each iteration
    disp(['Iteration ', num2str(iteration), ': Max voltage error = ', num2str(error)]);
end

% Display results
if error <= tol
    disp(['Converged in ', num2str(iteration), ' iterations']);
else
    disp('Did not converge within the maximum number of iterations');
end

%% Display Final Results
V_angle = angle(V) * 180 / pi;  % Convert angles from radians to degrees
fprintf('%-5s %-10s %-10s\n', 'Bus', 'Voltage (p.u.)', 'Angle (degrees)');
fprintf('-----------------------------\n');
for i = 1:5
    fprintf('%-6d %-6.4f %10.3f \n', i, abs(V(i)), V_angle(i));
end

%% Convert to Bus Voltages in kV
V_bus = abs(V) .* V_base;
disp('Bus Voltages in kV:');
fprintf('%-5s %-10s %-10s\n', 'Bus', 'Voltage (kV)', 'Angle (degrees)');
fprintf('-----------------------------\n');
for i = 1:5
    fprintf('%-6d %-.4f %9.3f\n', i, V_bus(i), V_angle(i));
end

%% Power Calculations for Each Bus
disp('Power Calculations:');
fprintf('%-5s %-10s %-10s\n', 'Bus', 'P (MW)', 'Q (MVAR)');  % Table headers
fprintf('-------------------------\n');  % Separator

for i = 1:5
    S_bus = V(i) * conj(Y_bus(i,:) * V);  % Calculate complex power for each bus
    fprintf('%-5d %-10.4f %-10.4f\n', i, real(S_bus) * S_base, imag(S_bus) * S_base);  % Output values
end
