% 4. For each of the above signals determine the number of FS coeffients
%    required to capture capture 99% of the power of the signal. Use this
%    criteria to determine the number of smapling points that you should
%    use in each case.

clear;
To = 10^-3; 
fo = 1/To; % fundamental frequency

N = 200;
dt = To/N;
t=(0:N-1)*dt; % time vector

% Define Signals -- Comment out the ones that are not needed
x = t; % Sawtooth Signal

%x = sin(2*pi*fo*t + pi / 2); % Half-Wave Recitfier
%x(x<0) = 0;

%x = cos(2*pi*fo*(t + (1/(2*To)).*t.^2)); % Chirp Signal

% Calculating FS Coefficients X_1 to X_N-1
X(1:N) = 0; % Set FS Coefficients to 0
for k=0:N-1
    X(1+k) = (1/N)*x*(exp(j*2*pi*k*fo*t))';  
end
%[X(1), X(2), X(3), ..., X(100)]

% Power of original signal
P_t = sum(x.^2*dt) / To;
% P_1 = To^2 / 3; % For testing

P_x = 0; % Define power of signal using Fourier Coefficients
numterms = 0;

for k=0:N-1 % Do it for the total number of terms
    P_x = P_x + (abs(X(k+1)))^2; % Parseval's Theorem

    % Check if power from Fourier Coefficients is enough to capture 99%
    % of power of original signal
    if P_x > 0.99*P_t
        numterms = k+1;
        break
    end
end

% numterms is the number of FS coefficients required
