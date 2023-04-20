% 3. Define the time signal x(t) on the interval [0,To], where To = 10^(-3)
%    as x(t) = cos(2*pi*fo(t + (1/(2*To))*t^2)) where fo=1/To. This is
%    known as a chirp signal. Find the Fourier series coefficients using a
%    suitable number of points in the sampled signal, i.e a suitable
%    sampling frequency.

To = 10^-3;  % Fundamental period of the signal. Default: 1 ms in this example.
fo = 1/To;   % Fundamental frequency

N = 100;         
dt = To/N; % Time intervals in t vector   
t=(0:N-1)*dt; 

% Chirp signal
x = cos(2*pi*fo*(t + (1/(2*To)).*t.^2));
figure
plot(t, x)
title("Chirp Signal")
xlabel("Time (s)")
ylabel("Amplitude")

Np=10;
tp = (0:Np*N-1)*dt;       % create a timeline over Np periods.

figure
plot(tp,periodicExtension(x,Np));  % Plot 3 periods
title("Periodic Extension of Chirp Signal")
xlabel("Time (s)")
ylabel("Amplitude")

X(1:N) = 0;   
for k=0:N-1
    X(1+k) = (1/N)*x*(exp(j*2*pi*k*fo*t))';  
end

df = fo; % Frequency increment between successive harmonics
f = (0:N-1)*df; % Frequency vectors    

% Plotting Magnitude of the Fourier Coefficients
figure
stem(f,abs(X)) % The FS coefficients are complex numbers. Here we plot the magnitude
title("Fourier Coefficients of Chirp Signal")
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
