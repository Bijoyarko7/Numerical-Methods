
%{
% Setup parameters
p_start = 0.2;
p_step = 0.02;
p_end = 10;
n_iter = length(p_start:p_step:p_end);
r = 2^(8/9) - 2^(1/9);

% Pre-allocate results
results_u1 = zeros(n_iter, 3);
results_u2 = zeros(n_iter, 3);

% Setup optimization options
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

% First household optimization
counter = 1;
for p = 0.2:0.02:10
    % Initial guess
    x0 = [1, 1];
    
    % Budget constraint: x1 + p*x2 <= 2 + p*r
    A = [1 p];
    b = 2 + p*r;
    
    % Bounds: non-negativity
    lb = [0, 0];
    ub = [];
    
    % Optimize
    [x, fval] = fmincon(@(x) -utility1(x(1), x(2)), x0, A, b, [], [], lb, ub, [], options);
    
    % Store results
    results_u1(counter,:) = [p, x(1), x(2)];
    counter = counter + 1;
end

% Second household optimization
counter = 1;
for p = 0.2:0.02:10
    % Initial guess
    x0 = [1, 1];
    
    % Budget constraint: x1 + p*x2 <= r + p*2
    A = [1 p];
    b = r + p*2;
    
    % Bounds: non-negativity
    lb = [0, 0];
    ub = [];
    
    % Optimize
    [x, fval] = fmincon(@(x) -utility2(x(1), x(2)), x0, A, b, [], [], lb, ub, [], options);
    
    % Store results
    results_u2(counter,:) = [p, x(1), x(2)];
    counter = counter + 1;
end

% Create tables and join results
colnames1 = {'Price', 'cons1_g1', 'cons1_g2'};
colnames2 = {'Price', 'cons2_g1', 'cons2_g2'};
opt_u1 = array2table(results_u1, 'VariableNames', colnames1);
opt_u2 = array2table(results_u2, 'VariableNames', colnames2);
results = join(opt_u1, opt_u2, 'Keys', 'Price');

% Calculate excess demand for good 2
results.resDemand_g2 = results.cons1_g2 + results.cons2_g2 - (r + 2);

% Plot excess demand
figure;
plot(results.resDemand_g2, results.Price, 'LineWidth', 1.5)
xlabel('Excess demand good 2')
ylabel('Price')
title('Residual Demand')
grid on

% Utility functions
function u = utility1(x1, x2)
    if x2 <= 0
        u = -inf;
    else
        u = x1 - (1/8)*x2^(-8);
    end
end

function u = utility2(x1, x2)
    if x1 <= 0
        u = -inf;
    else
        u = x2 - (1/8)*x1^(-8);
    end
end
%}


%Q2
%{
counter = 1;
n_iter = length(0:0.02:1);
results_soc = zeros(n_iter, 3);
r = 2^(8/9) - 2^(1/9);

% Set optimization options
options = optimoptions('fminunc',...
    'Display', 'off',...
    'MaxIterations', 1000,...
    'MaxFunctionEvaluations', 3000,...
    'OptimalityTolerance', 1e-6,...
    'StepTolerance', 1e-6);

for a = 0:0.02:1
    % Initial guess
    x0 = [1, 1];
    
    % Define objective function with error handling
    objFun = @(x) objectiveWithChecks(x, a, r);
    
    try
        % Optimize
        [x, fval] = fminunc(objFun, x0, options);
        
        % Store results
        results_soc(counter,:) = [a, x(1), x(2)];
    catch
        % Handle optimization failure
        results_soc(counter,:) = [a, NaN, NaN];
    end
    counter = counter + 1;
end
%}

options = optimoptions('fmincon',...
    'Display', 'off',...
    'Algorithm', 'interior-point',...
    'MaxIterations', 1000);

for a = 0:0.1:1
    % Initial guess
    x0 = [1, 1];
    
    % Constraints
    lb = [0, 0];  % Lower bounds
    ub = [r+2, r+2];  % Upper bounds
    
    try
        [x, fval] = fmincon(@(x) objectiveWithChecks(x, a, r), x0, [], [], [], [], lb, ub, [], options);
        results_soc(counter,:) = [a, x(1), x(2)];
    catch
        results_soc(counter,:) = [a, NaN, NaN];
    end
    counter = counter + 1;
end



% Objective function with checks
function f = objectiveWithChecks(x, a, r)
    x1_1 = x(1);
    x2_1 = x(2);
    x1_2 = r + 2 - x(1);
    x2_2 = r + 2 - x(2);
    
    % Check for invalid values
    if x2_1 <= 0 || x1_2 <= 0
        f = inf;
        return;
    end
    
    % utilities
    u1 = utility1(x1_1, x2_1);
    u2 = utility2(x1_2, x2_2);
    
    % Social welfare function
    f = -(a*u1 + (1-a)*u2);
end

% Utility functions
function u = utility1(x1, x2)
    if x2 <= 0
        u = -inf;
    else
        u = x1 - (1/8)*x2^(-8);
    end
end

function u = utility2(x1, x2)
    if x1 <= 0
        u = -inf;
    else
        u = x2 - (1/8)*x1^(-8);
    end
end

% Plot
figure;
plot(results_soc(:,2), results_soc(:,3), 'b.-', 'LineWidth', 2)
hold on
plot(2, r, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10)  % Endowment point
xlabel('Good 1')
ylabel('Good 2')
title('Pareto Optima and Initial Endowment')
grid on
legend('Pareto Optima', 'Initial Endowment')