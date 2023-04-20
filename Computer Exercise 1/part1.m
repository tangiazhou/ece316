% 1. Find the FS coefficients for a sawtooth wave. Use the period To =
%    10^(-3). Choose a suitable number of points and experiment with different
%    numbers of points, e.g. N=100, N=200. Now, determine the original signal from the FS coefficients
%    and see how the accuracy of the reconstructed signal depends on the
%    number of points which is a function of the sampling frequency.
%    Plot the magnitude and phase of the coefficients.
clear;
To = 10^-3; 
fo = 1/To; % fundamental frequency

N = 100;
dt = To/N;
t=(0:N-1)*dt; % Time vector of one period of the signal

% Define sawtooth signal 
x = t;
figure
plot(t, x)
title("One Period of the Sawtooth Signal")
xlabel("Time (s)")
ylabel("Amplitude")

% Extend to 10 periods
Np=10;
tp = (0:Np*N-1)*dt;       % create a timeline over Np periods.
%x = periodicExtension(x, Np);
figure
plot(tp, periodicExtension(x, Np));  % Do a stem plot
title("Periodic Extension of Sawtooth")
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
title("Magnitude of Fourier Coefficients")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
ymax = max(abs(X));
ylim([0 ymax])

% Plotting Phase of the Fourier Coefficients
figure
plot(f, angle(X))
title("Phase of Fourier Coefficients")
xlabel("Frequency (Hz)")
ylabel("Phase")

% Plotting Reconstruction of Signal with Fourier Coefficients
% number of points which is a function of the sampling frequency.

% Inverse Fourier Transform of Fourier Coefficients
x = N*ifft(X);

figure
plot(t, real(x))
title("Reconstruction of Signal Using Fourier Coefficients")
xlabel("Time(s)")
ylabel("Phase")

%FUNCTIONS
function y = periodicExtension(x,n)
% Return the periodic extension of the vector x.
% The number of periods is n
L = length(x);
y(1:n*L) = 0;
for k=1:n
    y((k-1)*L+1: k*L) = x;
end
end
