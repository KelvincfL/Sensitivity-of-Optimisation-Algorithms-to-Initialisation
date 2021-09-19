theta = pi/6;
h = 0.01;
grad_f = @(x) 4*x^3-3*x-cos(theta);
X = zeros(101,1);   % to store the final outcome
for i = 1:101
    x_a = (i-51)/50;        % set initial value
    for j = 1:1000
        x_b = x_a - h*grad_f(x_a);      % alg iteration
        x_a = x_b;
    end
    X(i) = x_b;
end
T=table([-1:0.02:1].',X);
disp(T)
table2latex(T,'q2table.tex')
