% Computer Exercise 2 Part A

% Part A
% Create three different alarm signals using the above program.
% a) Signal starting with low frequency and ending on high frequency
% b) Signal starting with high frequency and ending with low frequency
% c) Signal starting with low frequency then high frequency, then back to
%    low frequency.
% d) Signal Starting with high freq. then low frequency, then back to high
%    frequency.
% Pick reasonably good parameters for the frequencies and durations.

clear all
close all

% a) Signal starting with low frequency and ending on high frequency
To = 1;             % Pulse Period. This is the period of one of the pulses, i.e. chirps.
Tp = 5;             % Length of signal to play, in sec. This is the total length to play.
Np = floor(Tp/To);  % Number of periods
fo = 1/To;          % Fundamental Frequency
F1 = 100;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 2000;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

Fs = 2.5*max(F1,F2);  % Sampling rate. A few times greater than Nyquist rate
Nt = To*Fs; % Total number of samples at this sampling rate.
N = 2^(ceil(log2(Nt))); % Round upward to a power of 2 for efficient FFT implementation
dt = To/N;  % Sampling period at the rounded up rate.
Fs = 1/dt;
t=(0:N-1)*dt; % time vector

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

figure 
plot(t,x)     % show a plot of the signal with continuous time
title("Chirp Signal Low Frequency to High Frequency")
xlabel("Time (s)")
ylabel("Amplitude")

player = audioplayer(x, Fs);
play(player)

% b) Signal starting with high frequency and ending with low frequency
F1 = 2000;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 100;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

figure 
plot(t,x)     % show a plot of the signal with continuous time
title("Chirp Signal High Frequency to Low Frequency")
xlabel("Time (s)")
ylabel("Amplitude")

player = audioplayer(x, Fs);
play(player)

% c) Signal starting with low frequency then high frequency, then back to
%    low frequency.
To = 1;             % Pulse Period. This is the period of one of the pulses, i.e. chirps.
Tp = 5;             % Length of signal to play, in sec. This is the total length to play.
Np = floor(Tp/To);  % Number of periods
fo = 1/To;          % Fundamental Frequency
F1 = 100;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 2000;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

Fs = 2.5*max(F1,F2);  % Sampling rate. A few times greater than Nyquist rate
Nt = To*Fs; % Total number of samples at this sampling rate.
N = 2^(ceil(log2(Nt))); % Round upward to a power of 2 for efficient FFT implementation
dt = To/N;  % Sampling period at the rounded up rate.
Fs = 1/dt;
t=(0:N-1)*dt; % time vector

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

player = audioplayer(x, Fs);
play(player)

F1 = 2000;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 100;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

player = audioplayer(x, Fs);
play(player)

% d) Signal Starting with high freq. then low frequency, then back to high
%    frequency.
To = 1;             % Pulse Period. This is the period of one of the pulses, i.e. chirps.
Tp = 5;             % Length of signal to play, in sec. This is the total length to play.
Np = floor(Tp/To);  % Number of periods
fo = 1/To;          % Fundamental Frequency
F1 = 2000;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 100;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

Fs = 2.5*max(F1,F2);  % Sampling rate. A few times greater than Nyquist rate
Nt = To*Fs; % Total number of samples at this sampling rate.
N = 2^(ceil(log2(Nt))); % Round upward to a power of 2 for efficient FFT implementation
dt = To/N;  % Sampling period at the rounded up rate.
Fs = 1/dt;
t=(0:N-1)*dt; % time vector

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

player = audioplayer(x, Fs);
play(player)

F1 = 100;          % Starting Frequency. Start frequency for the chirp pulse.
F2 = 2000;           % Ending Frequency - End frequency for the chip pulse.
B = F2 - F1;        % Frequency span

x = cos(2*pi*F1*(t + (B/(2*To*F1))*t.^2)); % our example signal

player = audioplayer(x, Fs);
play(player)
