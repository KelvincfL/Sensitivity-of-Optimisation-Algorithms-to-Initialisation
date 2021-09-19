format long
A_2 = 1;        % setting constants
A_3 = 1;
A_4 = 1;
C = 10^7;
T = 10;
N(1) = (C*sqrt(A_2))/(A_3*T*(10*sqrt(A_2)+(1+sqrt(2))*sqrt(6*A_4)));
l=1;
while floor(N(l)) > 0       % calculates level sizes and number of levels
    l = l+1;
    N(l) = (C*sqrt(A_4))/((A_3*T*2^((3*l-1)/2))*(50*sqrt(3*A_2) ...
        + (30+15*sqrt(2))*sqrt(A_4)));    
end
N = floor(N);
L = numel(N)-2;     % highest level
results_p1 = zeros(1,2^6);      %stores value of p1
results_p2 = zeros(1,2^6);      %stores value of p2

for k = 1:2^6
    theta = k*pi/2^7;
    grad_f = @(x) 4*x^3-3*x-cos(theta);
    m1 = cos((theta+2*pi)/3);
    m2 = cos(theta/3);
    mu_estimator = 0;       %stores value of mu
% below is used to calculate Y_0    
    h_0 = 0.1*2^0;
    iter_required = T/h_0;
    n_0 = N(1);       
    Y = zeros(1,n_0);
    for i = 1:n_0
        x_a = 2*rand -1;         % random number in [-1,1]
        for j = 1:iter_required
            x_b = x_a - h_0*grad_f(x_a);      % alg iteration
            x_a = x_b;
        end
        Y(i) = x_a;
    end
    mu_estimator = mu_estimator + mean(Y);
% below is used to calculate Y_l, l not 0
    for l = 1:L
        n_l = N(l+1);        
        Y = zeros(1,n_l);        
        h_1 = 0.1*2^(-l);
        h_2 = 0.1*2^(-l+1);
        iter_h_1 = T/h_1;        % number of iterations for h_l
        iter_h_2 = T/h_2;        % number of iterations for h_(l-1)
        for i = 1:n_l
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
    results_p1(k) = (m2-mu_estimator)/(m2-m1);      % calculates p1
    results_p2(k) = (mu_estimator-m1)/(m2-m1);      % calculates p2
end
%below plots the graph
X = linspace(1,64,64);
plot(X,results_p1,"DisplayName","p_1")
hold on
plot(X,results_p2,"DisplayName","p_2")
hold off
legend
