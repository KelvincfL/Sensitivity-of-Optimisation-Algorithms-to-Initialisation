theta = pi/4;
grad_f = @(x) 4*x^3-3*x-cos(theta);
T = 10;
mu_estimator = zeros(11,1);     % store values of mu_h
var_estimator = zeros(11,1);      % store values of variance of X
for k = 0:10
    h = 0.1*2^(-k);
    iter_required = T/h;        % number of iterations per sample
    N = 2^(20-k);           % number of samples used
    temp_results=zeros(N,1);        % stores outcome of samples
    for n = 1:N
        x_a = 2*rand -1;         % random number in [-1,1]
        for i = 1:iter_required
            x_b = x_a - h*grad_f(x_a);      % alg iteration
            x_a = x_b;
        end
        temp_results(n) = x_a;
    end
    mu_estimator(k+1) = mean(temp_results);     % find the expectation using averages
    var_estimator(k+1) = var(temp_results);
end
table([0:10].',mu_estimator,var_estimator)
%table2latex(tbl,'q4.tex')
