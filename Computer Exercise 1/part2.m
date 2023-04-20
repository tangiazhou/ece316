% 2. Find the FS coefficients (plot magnitute and phase) for a half-wave rectified sine wave. Choose 
%    the origin so that the signal is even, and To = 10^(-3). Again
%    consider different choices for the number of points in one period.

To = 10^-3;  % Fundamental period of the signal. Default: 1 ms in this example.
fo = 1/To;   % Fundamental frequency

N = 100;           % Number of samples in [0,To]
dt = To/N;        % Sampling period. The sampling frequency is 1/dt. It should be greater or equal to the Nyquist sampling frequency.
t = (0:N-1)*dt;  

% Half-wave rectifier function
x = sin(2*pi*fo*t + pi / 2);
x(x<0) = 0;

figure
plot(t, x);
title("One Period of Half-Wave Rectifier")
xlabel("Time (s)")
ylabel("Amplitude")

Np=3;
tp = (0:Np*N-1)*dt;       % create a timeline over Np periods.

figure
plot(tp,periodicExtension(x,3));  % Plot 3 periods
title("Three Periods of Half-Wave Rectifier")
xlabel("Time (s)")
ylabel("Amplitude")

% Calculating FS Coefficients X_1 to X_N-1
X(1:N) = 0; % Set FS Coefficients to 0
for k=0:N-1
    X(1+k) = (1/N)*x*(exp(j*2*pi*k*fo*t))';  
end

df = fo; % Frequency increment between successive harmonics
f = (0:N-1)*df; % Frequency vectors    

% Plotting Magnitude of the Fourier Coefficients
figure
stem(f,abs(X)) % The FS coefficients are complex numbers. Here we plot the magnitude
title("Magnitude of the Fourier Coefficients")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
ymax = max(abs(X));
ylim([0 ymax])

% Plotting Phase of the Fourier Coefficients
figure
plot(f, angle(X))
title("Phase of the Fourier Coefficients")
xlabel("Frequency (Hz)")
ylabel("Phase")

function y = periodicExtension(x,n)
% Return the periodic extension of the vector x.
% The number of periods is n
L = length(x);
y(1:n*L) = 0;
for k=1:n
    y((k-1)*L+1: k*L) = x;
end
end
