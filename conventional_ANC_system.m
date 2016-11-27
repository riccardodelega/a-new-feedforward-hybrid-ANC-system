%% MIAS project 2015/2016 ("A New Feedforward Hybrid ANC System", Xiao, Wang, 2011)

close all, clear all, clc

simulation = input('Input 0 for broadband noise only, 1 for broadband and narrowband noise: ');
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
mu_h = 0.0007; % stepsize of FxLMS (0.1/L = 0.02)


% fast convergence and low-order adaptive filter
% L = 21; mu_h = 0.0018; %(0.1/L = 0.048)
% firework noise when broadband only
% mu_h = 0.01999;
% firework noise when mixed
mu_h = 0.0027;


%% Filters
P = [fir1(M_p-1, 0.4)]'; % primary path P(z)
S = [fir1(M-1, 0.4)]'; % secondary path S(z)
S_hat = zeros(M_hat,1); % secondary path estimation S^(z)
h = zeros(L,1); % adaptive control filter H(z)


%% Signals (see Fig. 1)
p_0 = zeros(N,1); % primary noise
e_0 = zeros(N,1);
e_0_prime = zeros(N,1);
e = zeros(N,1);
y = zeros(N,1);
x_r_hat = zeros(N,1);
v_p = randn(N,1)*sqrt(sigma_p); % uncorrelated noise acting on the error microphone

% compute reference noise x_r(n)
x_r = randn(N,1)*sqrt(sigma_w); % broadband component
if simulation == 1
    q = length(frequencies); % number of frequencies
    time = 1:N;
    for ii = 1:q
        x_a_i = [cos(frequencies(ii)*time)]';
        x_b_i = [sin(frequencies(ii)*time)]';
        x_r = x_r + a_r(ii)*x_a_i + b_r(ii)*x_b_i;
    end
    x_r = x_r/std(x_r);
    % normalise so that the signal power in both simulation is 1
    % bandpower(x_r) == 1 in both cases
    % sqrt(var(x)?y(n)^2/var(y)), which simplifies to std(x)/std(y)
end

%% Identification of S_hat(z)
mu_s = 0.001;
omega = 1; % variance of excitation signal
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
for n = first_sample:N
    % compute primary noise p_0(n)
    X_r = x_r(n:-1:n-M_p+1);
    p_0(n) = P'*X_r;
    
    % compute y(n)
    X_r = x_r(n:-1:n-L+1);
    y(n) = h'*X_r;
    
    e_0(n) = p_0(n) - y(n);
    
    % compute e_0_prime(n)
    E_0 = e_0(n:-1:n-M+1);
    e_0_prime(n) = S'*E_0;
    
    e(n) = e_0_prime(n) + v_p(n); % add (uncorrelated) measurement noise
    
    % compute reference signal
    X_r = x_r(n:-1:n-M_hat+1);
    x_r_hat(n) = S_hat'*X_r;
    
    % update filter
    X_r_hat = x_r_hat(n:-1:n-L+1);
    h = h + mu_h*e(n)*X_r_hat;
    
    if e(n) > 5 % diverging
       disp(['out at ', num2str(n)]), return
    end
end

%{d
figure
res_pow = e(first_sample:end).^2;
res_pow(res_pow < sigma_p/100) = sigma_p/100;
subplot(211), plot(20*log10(res_pow))
xlabel('time'), ylabel('dB'), title('Residual noise power')
subplot(212), plot(e(first_sample:end))
xlabel('time'), ylabel('linear'), title('Residual noise')
%}
%{d
figure
res_pow = e(first_sample:end).^2;
res_pow(res_pow < sigma_p/100) = sigma_p/100;
plot(20*log10(res_pow)), xlabel('time'), ylabel('dB'), title('Residual noise power')
%}
%{
figure
[~,F,T,P] = spectrogram(e);
imagesc(T,F,10*log(abs(P))), colorbar
%}

if play == 1
    sound(e(first_sample*2:end), 20000)
end
