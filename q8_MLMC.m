format long
T = 10;
L = 10;         %fixed truncation level
results=zeros(2^6,1);       
for k = 1:2^6
    theta = k*pi/2^7;
    grad_f = @(x) 4*x^3-3*x-cos(theta); 
    mu_estimator = 0;       %stores value of mu
% below is used to calculate Y_0    
    h_0 = 0.1*2^0;
    iter_required = T/h_0;
    N = 2^10;       % change to 5 for Q8 part i)
    Y = zeros(1,N);
    for i = 1:N
        x_a = 2*rand -1;         % random number in [-1,1]
        for j = 1:iter_required
            x_b = x_a - h_0*grad_f(x_a);      % alg iteration
            x_a = x_b;
        end
        Y(i) = x_a;
    end
    mu_estimator = mu_estimator + mean(Y);
% below is used to calculate Y_l, l not 0
    for l = 1:10
        N = 2^(L-l);        % change to 5 for Q8 part i)
        Y = zeros(1,N);        
        h_1 = 0.1*2^(-l);
        h_2 = 0.1*2^(-l+1);
        iter_h_1 = T/h_1;        % number of iterations for h_l
        iter_h_2 = T/h_2;        % number of iterations for h_(l-1)
        for i = 1:N
            x_0 = 2*rand -1;         
         
            x_a = x_0;      %copy of initial point for X_(h_l)
            for j = 1:iter_h_1
                x_b = x_a - h_1*grad_f(x_a);      
                x_a = x_b;
            end
         
            x_c = x_0;      %copy of initial point for X_(h_(l-1))
            for j = 1:iter_h_2
                x_d = x_c - h_2*grad_f(x_c);      
                x_c = x_d;
            end
            Y(i) = x_a - x_c;
        end
        mu_estimator = mu_estimator + mean(Y);     % sum up the mean of Y_l
    end
    results(k) = mu_estimator;
end
T=table([1:2^6].',results);
disp(T)
