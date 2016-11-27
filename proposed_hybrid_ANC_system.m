%% MIAS project 2015/2016 ("A New Feedforward Hybrid ANC System", Xiao, Wang, 2011)

close all, clear all, clc

play = 0; % play error signal after noise cancellation
verbose_identification = 0; % show residual error and S_hat(z) estimation accuracy

%% Parameters
N = 40000; % number of iterations
M_p = 41; % P(z) order
M = 21; % S(z) order
M_hat = 31; % S_hat(z) order
L = 51; % H(z) order (adaptive control filter)
frequencies = [0.03*pi; 0.06*pi; 0.09*pi]; % sinusoidal noise component
a_r = [2.0; 1.0; 0.5]; % DFCs 
b_r = [-1.0; -0.5; 0.1]; % DFCs
sigma_p = 0.1; % variance of uncorrelated noise affecting e(n)
sigma_w = 1.0; % variance of reference broadband noise component
mu_h = 0.0007; % broadband stepsize
mu_c_min = 0.00001; 
alpha = 0.9998;
mu_c = zeros(N,1); % sinusoidal noise canceller stepsize
mu_c_max = 0.05; % initial value
mu_n = 0.01; % narrowband stepsize

% fast convergence and low-order adaptive filter
% L = 21; mu_h = 0.0018; %(0.1/L = 0.048)
% firework noise 
% mu_h = 0.010;

%% Filters
P = [fir1(M_p-1, 0.4)]'; % primary path P(z)
S = [fir1(M-1, 0.4)]'; % secondary path S(z)
S_hat = zeros(M_hat,1); % secondary path estimation
h = zeros(L,1); % adaptive control filter
a_c_hat = zeros(length(frequencies),1); % adaptive coefficients (Sinusoidal Noise Canceller)
b_c_hat = zeros(length(frequencies),1);
a_hat = zeros(length(frequencies),1); % adaptive coefficients (Narrowband ANC Subsystem)
b_hat = zeros(length(frequencies),1);


%% Signals (see Fig. 1)
p_0 = zeros(N,1); % primary noise
e_0 = zeros(N,1);
e_0_prime = zeros(N,1);
e = zeros(N,1);
y_w = zeros(N,1);
v_p = randn(N,1)*sqrt(sigma_p); % uncorrelated noise acting on the error microphone
x = zeros(N,1); % reference noise after SNC action (only wideband component remains)
x_hat = zeros(N,1);
y = zeros(N,1); % = y_w + y_f
y_f = zeros(N,1);


% compute reference noise x_r(n)
x_r = randn(N,1)*sqrt(sigma_w); % broadband component
q = length(frequencies); % number of frequencies
time = 1:N;
for ii = 1:q
    x_a_i = [cos(frequencies(ii)*time)]';
    x_b_i = [sin(frequencies(ii)*time)]';
    x_r = x_r + a_r(ii)*x_a_i + b_r(ii)*x_b_i;
end
x_r = x_r/std(x_r);
% normalise so that the signal power is 1

x_a = zeros(N, q); % reference signals for narrowband ANC
x_b = zeros(N, q); 
x_a_hat = zeros(N, q); % reference signals convoluted with S^(z)
x_b_hat = zeros(N, q);
for ii = 1:q
    x_a(:, ii) = cos(frequencies(ii)*time);
    x_b(:, ii) = sin(frequencies(ii)*time);
end



%% Identification of S_hat(z)
mu_s = 0.001;
omega = 1;
N_s = 10000;
x_s = randn(N_s,1)*omega; % input noise for identification
d_s = zeros(N_s,1);
y_s = zeros(N_s,1);
e_s = zeros(N_s,1);
disp('Identifying S_hat(z)...')
for n = M_hat+1:N_s
    % compute d_s(n)
    X_s = x_s(n:-1:n-M+1);
    d_s(n) = S'*X_s;
    % compute y_s(n)
    X_s = x_s(n:-1:n-M_hat+1);
    y_s(n) = S_hat'*X_s;
    % compute e(n) and update weights
    e_s(n) = d_s(n) - y_s(n);
    S_hat = S_hat + mu_s*e_s(n)*X_s;
end
disp('Identification complete.')

if verbose_identification == 1
    figure
    plot(e_s), axis([1 N_s -1.1 1.1])
    title('e(n)')
    figure, freqz(S,1)
    figure, freqz(S_hat,1)
    pause
    close all
end


%% Algorithm (conventional broadband ANC)
first_sample = max([M_p, M, L, M_hat]);
mu_c(first_sample) = mu_c_max;
for n = first_sample:N
    %% sinusoidal noise canceller (SNC) subsystem
    x(n) = x_r(n);
    for ii = 1:q
        x(n) = x(n) - a_c_hat(ii)*x_a(n, ii) - b_c_hat(ii)*x_b(n, ii);
    end
    
    % update coefficients and stepsize
    for ii = 1:q
        a_c_hat(ii) = a_c_hat(ii) + mu_c(n)*x(n)*x_a(n, ii);
        b_c_hat(ii) = b_c_hat(ii) + mu_c(n)*x(n)*x_b(n, ii);
    end
    mu_c(n+1) = alpha*mu_c(n) + (1-alpha)*mu_c_min; % adaptive stepsize
    
    
    %% broadband and narrowband ANC subsystems
    % compute primary noise p_0(n)
    X_r = x_r(n:-1:n-M_p+1);
    p_0(n) = P'*X_r;
    
    % compute y_w(n)
    X_r = x_r(n:-1:n-L+1);
    y_w(n) = h'*X_r;
    
    % compute y_f(n)
    for ii = 1:q
        y_f(n) = y_f(n) + a_hat(ii)*x_a(n, ii) + b_hat(ii)*x_b(n, ii);
    end
    
    % compute y(n)
    y(n) = y_f(n) + y_w(n);
    
    % compute e_0(n)
    e_0(n) = p_0(n) - y(n);
    
    % compute e_0_prime(n)
    E_0 = e_0(n:-1:n-M+1);
    e_0_prime(n) = S'*E_0;
    
    e(n) = e_0_prime(n) + v_p(n); % add (uncorrelated) measurement noise
    
    % compute broadband reference signal
    X = x(n:-1:n-M_hat+1);
    x_hat(n) = S_hat'*X;
    
    % convolve narrowband reference signals
    for ii = 1:q
        X_a = x_a(n:-1:n-M_hat+1, ii);
        x_a_hat(n, ii) = S_hat'*X_a;
        X_b = x_b(n:-1:n-M_hat+1, ii);
        x_b_hat(n, ii) = S_hat'*X_b;
    end
    
    % update broadband filter H(z)
    X_hat = x_hat(n:-1:n-L+1);
    h = h + mu_h*e(n)*X_hat;
    
    % update narrowband coefficients
    for ii = 1:q
        a_hat(ii) = a_hat(ii) + mu_n*e(n)*x_a_hat(n, ii);
        b_hat(ii) = b_hat(ii) + mu_n*e(n)*x_b_hat(n, ii);
    end
    
    if e(n) > 5 % diverging
       disp(['out at ', num2str(n)]), return
    end
    %}
end

%{d
figure
res_pow = e(first_sample*2:end).^2;
res_pow(res_pow < sigma_p/100) = sigma_p/100;
subplot(211), plot(20*log10(res_pow)), title('Residual noise power')
subplot(212), plot(e(first_sample*2:end)), title('Residual noise')%, axis([1, N, -1 1])
%subplot(313)
%[spectrum_e,omega] = periodogram(e(first_sample*2:end));
%plot(omega, 20*log10(spectrum_e)), title('Residual noise spectrum')
%}

%{d
figure
res_pow = e(first_sample:end).^2;
res_pow(res_pow < 0.001) = 0.001;
plot(20*log10(res_pow)), xlabel('time'), ylabel('dB'), title('Residual noise power')
%}
if play == 1
    sound(e(first_sample:end), 20000)
end

%{
% check if SNC works with the given alpha
close all,  figure
[spectrum_x, omega] = periodogram(x(first_sample*2:end));
plot(omega, spectrum_x)
%}