close all;

% Part 1
% Consider 16QAM with a raised cosine pulse shape with 𝛼 = 0.5. 
% You should truncate the pulse so that p(t) = 0 for |t| > 4T. 
% Now take the following text, 
% "I am doing the computer exercises 5 for ece316 – Fall 2020." 
% Translate the above text into a binary sequence using the ASCII code. 
% Include all the characters including spaces except the quotes. After you translate into a 
% binary sequence use a serial to parallel converter to produce two streams. Even bits go 
% to the In-phase channel (cos), and odd bits go to the quadrature channel (sin). 
% Then produce a 16QAM signal using the above pulse shape. Use a sampling period 
% equal to T/64 to submit an array with the signal. The carrier frequency of the signal 
% should be fc = 8/T. 
% Submit your code along with a file containing an array with the QAM signal. 

text = 'I am doing the computer exercises 5 for ece316 – Fall 2020.'; 

% Convert the string to ASCII binary 
% 8 bits * 59 elements
binary_msg = dec2bin(unicode2native(text, 'US-ASCII'), 8);

% Combine everything into one vector
binary = reshape(binary_msg.',1,[]) - '0';

% Length of the binary vector
msg_len = length(binary);

% Split the binary vector into even and odd streams
binary_even = binary(2:2:msg_len); % in phase 
binary_odd = binary(1:2:msg_len); % quadrature

% Mapping the bits to pulse amplitudes
% 00: -3
% 01: -1
% 10: 1
% 11: 3
amplitude_i = 2*(2*binary_even(1:2:msg_len/2)-1.5 + binary_even(2:2:msg_len/2));
amplitude_q = 2*(2*binary_odd(1:2:msg_len/2)-1.5 + binary_odd(2:2:msg_len/2));

% Time vector shifted by 4 as there are no negative indices 
t = (0:1/64:8) - 4;
a = 0.5;

% Using the cosine rolloff formula and plugging in the original time vector
% which is supposed to be t = (0:T/64:8*T) - 4T, we can factor out the T which
% cancels out with the Ts in the denominator of the formula to get the 
% formula for the pulse below:
pulse = (sin(pi*t).*cos(a*pi*t)) ./ ((pi*t).*(1-(2*a*t).^2));

% There are possible discontinuities in the pulse formula, so adjust them
pulse(isnan(pulse))=1;
pulse(isinf(pulse))=0;

% Plot of the truncated cosine rolloff pulse
figure
plot(t, pulse)
title("Plot of Truncated Cosine Rolloff Pulse");
xlabel("Time (In Terms of T)")
ylabel("Amplitude")

% Total number of samples in the entire transmitted signal 
% There are 118 symbols being transmitted and each symbol is 8T, so the
% total length of the transmitted signal is 125T
% Since there are 64 samples in T, we have 125*64 + 1
num_samples = (125*64) + 1; 

% Vectors containing indices for all the samples 
sample_i = zeros(1, num_samples);
sample_q = zeros(1, num_samples);

% Number of samples in a pulse (8*64)+1
pulse_samples = length(t);

% Adding all 118 symbols together to create the transmitted signal
for i = 1:118
    index_1 = (i-1)*64+1; % Beginning sample number for symbol
    index_2 = index_1 + pulse_samples -1; % Ending sample number for symbol
    sample_i(index_1:index_2) = sample_i(index_1:index_2) + amplitude_i(i) * pulse;
    sample_q(index_1:index_2) = sample_q(index_1:index_2) + amplitude_q(i) * pulse;
end

% QAM carrier
% cos(2*pi*fc*t) where fc = 8/T, and t is from (0:T/64:125*T)
% The T can be taken out and simplied to the below:
cos_carrier = cos(16*pi*(0:1/64:125));
sin_carrier = sin(16*pi*(0:1/64:125));

% Modulated QAM signal
modulated_signal = sample_i.*cos_carrier + sample_q.*sin_carrier;
figure
plot(modulated_signal)
title("Plot of Transmitted Signal");
xlabel("Samples")
ylabel("Amplitude")

% Part 2
% Take the above QAM signal and demodulate it using a QAM demodulator. Draw the block 
% diagram for the QAM demodulator. First Use a local oscillator signal that is synchronized 
% with the carrier of the received signal. Then use a local Oscillator that has a phase error 
% equal 𝜋/8, i.e. the local oscillator is cos(2𝜋𝑓_ct + 𝜋/8)

% Note that you should implement the 
% low pass filters using convolution in the time domain with the appropriate impulse 
% response, e.g. some sinc() pulse. After demodulating the signals plot them on the plane. 
% Then perform a transformation on the output signals in order to compensate for the effect 
% of carrier phase error. 

% Demodulating using a synchronized oscillator
%channel1 = pulse_sample.*cos(16*pi*(0:1/64:125));
%channel2 = pulse_sample.*sin(16*pi*(0:1/64:125));

% Demodulating with an oscillator with a phase error
channel1 = modulated_signal.*cos(16*pi*(0:1/64:125) + pi/8);
channel2 = modulated_signal.*sin(16*pi*(0:1/64:125) + pi/8);

% Low pass filter
N = num_samples; % Number of samples in the signal

fft_1 = fft(channel1);
fft_2 = fft(channel2);

filter(1:N) = 1;
Nb = floor(((1+a)*N)/128);
filter(Nb+2:N-Nb) = 0; % Set the middle area to be 0 to extract the baseband signal

% Channel outputs with phase errors
channel1out = ifft(filter.*fft_1);
channel2out = ifft(filter.*fft_2);

figure
plot(channel1out)
title("Channel 1 Output with Phase Error");
xlabel("Samples")
ylabel("Amplitude")

figure
plot(channel2out)
title("Channel 2 Output with Phase Error");
xlabel("Samples")
ylabel("Amplitude")

% Compensating for phase error
channel1correct = (channel1out + sample_q*sin(pi/8)) / cos(pi/8);
channel2correct = (channel2out - sample_i*sin(pi/8)) / cos(pi/8);

figure
plot(channel1correct)
title("Channel 1 Output without Phase Error");
xlabel("Samples")
ylabel("Amplitude")

figure
plot(sample_i)
title("Original In Phase Message Signal");
xlabel("Samples")
ylabel("Amplitude")

figure
plot(channel2correct)
title("Channel 2 Output without Phase Error");
xlabel("Samples")
ylabel("Amplitude")

figure
plot(sample_q)
title("Original Quadrature Message Signal");
xlabel("Samples")
ylabel("Amplitude")