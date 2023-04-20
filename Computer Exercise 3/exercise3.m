% Computer Exercise 3

clear;

% Note: Fantasie in C was created for a project last semester so it sounds like 8 bit music
infoWAV= audioinfo("Fantasie in C.wav");
[SigTime,Fs] = audioread("Fantasie in C.wav"); 
Duration = infoWAV.Duration;
N = infoWAV.TotalSamples;

dt = 1/Fs; % Time intervals
%nBits = infoWAV.BitsPerSample;

df = 1/(N*dt); % Frequency steps
t = (0:N-1)*dt; % Time vector
f = (0:N-1)*df; % Frequency Vector

% 1. Extract the two channels. Find the average power of the signal in the two channels.
%    P1 and P2. Listen to each channel and note any differences.
% Stereo Channels
Channel1 = SigTime(:,1);
Channel2 = SigTime(:,2);

P1 = 10*log10(Channel1'*Channel1/N); % Power of channel 1
P2 = 10*log10(Channel2'*Channel2/N); % Power of channel 2

soundsc(Channel1, Fs);
soundsc(Channel2, Fs);

% 2. Form the sum of the two channels. Determine the average power. Listen to it.
%    Below if we refer to an operate on "the signal" you should use the sum
%    signal.
ChannelSum = Channel1 + Channel2;

% Average power of the summation of the two channels
Pp = 10*log10((Channel1+Channel2)'*(Channel1+Channel2)/N);

soundsc(ChannelSum, Fs);

% 3. Form the difference of the two channels. Determine the average power.
%    Note that when two signals are correlated then the power of the sum is not equal to
%    the sum of the powers.
ChannelDiff = Channel1 - Channel2;

% Average power of the difference
Pm = 10*log10((Channel1-Channel2)'*(Channel1-Channel2)/N);

soundsc(ChannelDiff, Fs);

% 4. Perform an FFT on the signal and plot the absolute value. Note the
%    spectral content. What is the bandwidth of the signal using a 99%
%    power criteria?

% If N is odd, subtract 1
if rem(N,2) == 1
    N1 = N-1;
else
    N1 = N;
end

SigFreq = fft(ChannelSum(1:N1)); 
figure;
plot(f(1:length(SigFreq))/10^3,20*log10(abs(SigFreq))); % Plot signal in frequency domain - in dB
title("abs(FFT) in dB");

p_total = sum(abs(SigFreq).^2); % Total power
p = 0;
BW = 0;

for i=1:0.5*N1 % Go to half of the spectrum
    p = p + 2*abs(SigFreq(i))^2; % Add the power at that frequency 
    if p/p_total >= 0.99 % Check if p reaches 99% of the total power
        BW = i*df;
        break;
    end
end

% 5. Filter The signal to a bandwidth of 10 KHz, 5 KHz, 3HKz, and 1 KHz.
%    Listen to it in each case and comment of the changes.

% 10kHz
filter10 = 0;

if(filter10)
    Filter(1:N) = 0; % Set vector to 0
    Bf = 10; % Bandwidth to 10kHz
    Nb = floor(Bf*10^3/df); % Number of bins in terms of the FFT
    
    Filter(1:N)=1; 
    Filter(Nb+2:N-Nb)=0;
    
    % Plot the LPF
    figure
    plot(f/10^3,Filter);
    title("Filter Transfer Function: Low Pass 10kHz");
    
    output10 = Filter'.*SigFreq;
    figure
    plot(f(1:length(SigFreq))/10^3, abs(output10));
    title("Output of LPF 10kHz");
    soundsc(ifft(output10), Fs);
    % Since the cutoff frequency is relatively large, the sound is very similar
    % to the original signal
end

% 5kHz
filter5 = 0;

if(filter5)
    Filter(1:N) = 0; % Set vector to 0
    Bf = 5; % Bandwidth to 10kHz
    Nb = floor(Bf*10^3/df); % Number of bins in terms of the FFT
    
    Filter(1:N)=1; 
    Filter(Nb+2:N-Nb)=0;
    figure;
    plot(f/10^3,Filter);
    title("Filter Transfer Function: Low Pass 5kHz");
    
    output5 = Filter'.*SigFreq;
    plot(f(1:length(SigFreq))/10^3, abs(output5));
    title("Output of LPF 5kHz");
    soundsc(ifft(output5), Fs);
    % Although there isn't much of a difference, you can hear that not all of
    % the high pitches are being played
end

% 3kHz
filter3 = 0;

if(filter3)
    Filter(1:N) = 0; % Set vector to 0
    Bf = 3; % Bandwidth to 10kHz
    Nb = floor(Bf*10^3/df); % Number of bins in terms of the FFT
    
    Filter(1:N)=1; 
    Filter(Nb+2:N-Nb)=0;
    figure;
    plot(f/10^3,Filter);
    title("Filter Transfer Function: Low Pass 3kHz");
    
    output3 = Filter'.*SigFreq;
    plot(f(1:length(SigFreq))/10^3, abs(output3));
    title("Output of LPF 3kHz");
    soundsc(ifft(output5), Fs);
    % There is more of a difference compared with the original sound, with more
    % high frequencies removed
end

% 1kHz
filter1 = 0;

if(filter1)
    Filter(1:N) = 0; % Set vector to 0
    Bf = 1; % Bandwidth to 1kHz
    Nb = floor(Bf*10^3/df); % Number of bins in terms of the FFT
    
    Filter(1:N)=1; 
    Filter(Nb+2:N-Nb)=0;
    figure;
    plot(f/10^3,Filter);
    title("Filter Transfer Function: Low Pass 1kHz");
    
    output1 = Filter'.*SigFreq;
    plot(f(1:length(SigFreq))/10^3, abs(output1));
    title("Output of LPF 1kHz");
    soundsc(ifft(output1), Fs);
    % You can clearly tell that the output is much deeper and lower in
    % frequency
end

% 6. Now modulate the signal with a very small carrier frequency, e.g.
%    100 Hz and listen to the signal. This is what happens if ther is a
%    frequency error in the demodulation of an audio signal.
Modulation = 0;

if(Modulation)
    fc = 100;
    modulated_signal = ChannelSum.*cos(2*pi*fc*t');
    plot(t, modulated_signal);
    title("Modulated Signal with Frequency Offset");
    xlabel("Time (s)")
    ylabel("Signal Value")
    soundsc(modulated_signal, Fs)
end

% 7. Perform the Hilbert transform of the sum signal and listen to it and
%    the original signal. Is there a different. You should also look at a
%    small segment of the signal (e.g. over 100 samples) and compare the
%    original signal with the Hilbert Transformed signal.
HT = 0;

if(HT)
    HTFilter(1)=0;
    HTFilter(2:N/2) = -1i;
    HTFilter(N/2+1:N) = 1i;
    SigFreq = SigFreq'.*HTFilter;
    SigTime = ifft(SigFreq);

    % Play the Hilbert Transform signal
    soundsc(real(SigTime), Fs);

    % Play the original signal
    soundsc(ChannelSum, Fs);

    % Comparison of the original signal and the Hilbert transform using 100
    % points
    figure;
    plot(1:100, real(SigTime(1:100)), 1:100, ChannelSum(1:100))

end

% 8. Play with the graphic equalizer. Enhance the low frequency components,
%    Enhance the high frequency components, etc.
Eq = 0;

if(Eq)
    % Enhancing the low frequency components
    Equalizer(1:N)=1;
    Bf = 1; % Set bandwidth to 1kHz
    Nb = floor(Bf*10^3/df); % Number of bins in terms of the FFT

    % Equalizer Points
    EP = floor([ .1 1 ]*Nb);

    % Equalizer_Gains 
    EG = [ 10 1 ];

    Equalizer(1:EP(1))=EG(1);
    Leq = length(EP);
    for k=2: Leq
        Equalizer(EP(k-1)+1: EP(k)) = EG(k);
    end  

    figure
    xlim([0 N/2]);
    ylim([0 max(EG)]);
    figure
    plot(f(1:Nb)/10^3,Equalizer(1:Nb));
    title("Graphic Equalizer Enhancing Lower Frequency Components");

    output_low = Equalizer'.*SigFreq;
    figure
    plot(f(1:length(SigFreq))/10^3, abs(output_low));
    title("Output of Graphic Equalizer Enhancing Lower Frequency Components");

    % Enhancing the high frequency components
    Equalizer(1:N)=1;

    % Equalizer Points
    EP = floor([ .9 1 ]*Nb);

    % Equalizer_Gains 
    EG = [ 1 10 ];

    Equalizer(1:EP(1))=EG(1);
    Leq = length(EP);
    for k=2: Leq
        Equalizer(EP(k-1)+1: EP(k)) = EG(k);
    end  

    figure
    xlim([0 N/2]);
    ylim([0 max(EG)]);
    plot(f(1:Nb)/10^3,Equalizer(1:Nb));
    title("Graphic Equalizer Enhancing Higher Frequency Components");

    output_high = Equalizer'.*SigFreq;
    figure
    plot(f(1:length(SigFreq))/10^3, abs(output_high));
    title("Output of Graphic Equalizer Enhancing Higher Frequency Components");
end

