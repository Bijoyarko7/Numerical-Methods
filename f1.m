function y = f1(x)
if x == 0
    y = 0;
else
    y = x^4 * exp(-(x^2)/4) * sin(8*x^(-3));
end