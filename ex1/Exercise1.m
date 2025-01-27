%Numerical_Methods
%Problem Set 1
%Author: Bijoy Ratan Ghosh

%Problem_1

%{

h_values = zeros(16,1);    
er_one = zeros(16,1); %error for one sided finite difference       
er_two = zeros(16,1); %error for two sided finite difference
er_st = zeros(16,1);  %error for four point stencil
er_cs = zeros(16,1);  %error for complex step differentiation

for n = 1:16
    h_values(n) = 10^(-n);
    er_one(n) = abs((f1(1+h_values(n))-f1(1))/h_values(n) - exp(-1/4)*(3.5*sin(8)-24*cos(8)));
    er_two(n) = abs((f1(1+h_values(n))-f1(1-h_values(n)))/(2*h_values(n)) - exp(-1/4)*(3.5*sin(8)-24*cos(8)));
    er_st(n) = abs((-f1(1+2*h_values(n))+8*f1(1+h_values(n))-8*f1(1-h_values(n))+f1(1-2*h_values(n)))/(12*h_values(n)) - exp(-1/4)*(3.5*sin(8)-24*cos(8)));
    er_cs(n) = abs(imag(f1(1+h_values(n)*1i))/h_values(n) - exp(-1/4)*(3.5*sin(8)-24*cos(8)));
end


% Create a figure with a good size
figure('Position', [100, 100, 800, 600])
x=1
% Plot all errors on log-log scale
loglog(h_values, er_one, 'b-o', ...     % Blue line with circles for one-sided
       h_values, er_two, 'r-s', ...     % Red line with squares for two-sided
       h_values, er_st, 'g-^', ...      % Green line with triangles for stencil
       h_values, er_cs, 'm-d', ...      % Magenta line with diamonds for complex step
       'LineWidth', 2, 'MarkerSize', 8)

% Add labels and title
xlabel('Step size (h)', 'FontSize', 12)
ylabel('Absolute Error', 'FontSize', 12)
title('Comparison of Derivative Approximation Methods', 'FontSize', 14)

% Add a grid to make it easier to read values
grid on

% Add legend with descriptive names
legend('One-sided difference', ...
       'Two-sided difference', ...
       'Four-point stencil', ...
       'Complex step', ...
       'Location', 'northwest')

% Make the plot easier to read with white background
set(gca, 'Color', 'w')


fprintf("\n nice job!")



figure('Position', [100, 100, 800, 600])
% Regular plot (linear scale) using the same markers and colors
plot(h_values, er_one, 'b-o', ... % Blue line with circles for one-sided
     h_values, er_two, 'r-s', ... % Red line with squares for two-sided
     h_values, er_st, 'g-^', ... % Green line with triangles for stencil
     h_values, er_cs, 'm-d', ... % Magenta line with diamonds for complex step
     'LineWidth', 2, 'MarkerSize', 8)

% Essential labels and title
xlabel('Step size (h)')
ylabel('Absolute Error')
title('Comparison of Methods (Linear Scale)')

% Legend in the same position as before
legend('One-sided difference', ...
       'Two-sided difference', ...
       'Four-point stencil', ...
       'Complex step', ...
       'Location', 'northwest')

% Keep the grid for readability
grid on

%}

%Problem_2
er_hermite = zeros(7,1);
er_mc = zeros(4,1);
iter = transpose([2 3 5 9 11 15 31]);

for count=1:7
    n= iter(count,1);
    er_hermite(count) = abs(2*exp(0.5) - (2/sqrt(pi))*gausshermi(@f_0,0,0,n));
end

% Define sample sizes
sizes = [10 100 1000 10000];

for i = 1:length(sizes)
    % Generate random sample
    sample = randn(sizes(i), 1);
    er_mc(i) = abs(2*exp(0.5) - mean((sample.^2).*exp(sample)));
end


% First plot: Gauss-Hermite
figure('Position', [100 100 400 400])
loglog(iter, er_hermite, 'b-o', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Gauss-Hermite Quadrature Error')
xlabel('Number of Quadrature Points')
ylabel('Absolute Error')

% Second plot: Monte Carlo
figure('Position', [500 100 400 400])
loglog(sizes, er_mc, 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Monte Carlo Integration Error')
xlabel('Sample Size')
ylabel('Absolute Error')


% Function for hermite
function y = f_0(x)
 y = (x.^2).*exp(sqrt(2)*x);
end