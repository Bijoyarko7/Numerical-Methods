%Problem_1

%{
h_values = zeros(16,1);    
er_one = zeros(16,1);       
er_two = zeros(16,1);
er_st = zeros(16,1);

for n = 1:16
    h_values(n) = 10^(-n);
    er_one(n) = abs((f1(1+h_values(n))-f1(1))/h_values(n) - exp(-1/4)*(3.5*sin(8)+24*cos(8)));
    er_two(n) = abs((f1(1+h_values(n))-f1(1-h_values(n)))/2*h_values(n) - exp(-1/4)*(3.5*sin(8)+24*cos(8)));
    er_st(n) = abs((-f1(1+2*h_values(n))+8*f1(1+h_values(n))-8*f1(1-h_values(n))+f1(1-2*h_values(n)))/12*h_values(n) - exp(-1/4)*(3.5*sin(8)+24*cos(8)));
end

% Function
function y = f1(x)
    if x == 0
        error('x cannot be zero due to negative power in sine term')
    end
    y = x^4 * exp(-(x^2)/4) * sin(8*x^(-3));
end
%}

%Problem_2
er_hermite = zeros(7,1);
iter = transpose([2 3 5 9 11 15 31]);

for i=1:7
    n= iter(i,1)
    er_hermite(i) = abs(2*exp(0.5) - gausshermi(@f_0,0,0,n));
end

% Function for hermite
function y = f_0(x)
 y = (x.^2).*exp(sqrt(2)*x);
end


