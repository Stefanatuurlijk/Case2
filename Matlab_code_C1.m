%% MATLAB CODE 1.1:FFT
T=1;                    % some epoch length
fs=100;                 % a high enough sampling freq.
Ts=1/fs;                % sampling period
t=0:Ts:T;
x=cos(2*pi*5*t)+cos(2*pi*10*t)/2;
%x=x.*hann(length(x))'; % windowed
N=length(x);
M=2^11;                 % take M ? N
f=(0:M/2)*fs/M;         % gridded [0, fnyq]
xhat=fft(x,M)/fs;       % = Ts*fft(x,M)
xhat=xhat(1+(0:M/2));   % only need first half
plot(f,abs(xhat));      % plot of |xˆ(f )|
% MATLAB CODE 1.2: FFT
xa=abs(fft(x,M))*2/N; % amplitudes
xa=xa(1+(0:M/2)); % first half
xa(1)=xa(1)/2; % freq 0 is special
plot(f,xa)
xlabel('Frequency (Hz)')
ylabel('Unit [x]')
%% MATLAB CODE 1.3: Periodogram
[S,fm]=periodogram(x,[],[],fs); % S=periodogram
plot(fm,S) % fm=frequency
xlabel('Frequency (Hz)')
ylabel('[x]^2/[f]')
powx=sum(x.^2)/N % power of x
pows=sum(S)*(fm(2)-fm(1)) % power via periodogram
% We can also code the periodogram ourselves:
% M=2^11; % as before
% f=(0:M/2)*fs/M; % as before
% P=abs(fft(x,M)).^2/fs/N; % classic periodogram
% P=2*P(1+(0:M/2)); % one sided version:
% P(1)=P(1)/2; % halve 1st entry
% P(end)=P(end)/(2-rem(M,2)); % halve if M is even
% plot(f,P);
% xlabel('frequency (Hz)')
% ylabel('unit [x]')
% powp=sum(P)*(f(2)-f(1)) % same power again

%% MATLAB CODE 1.4: Welch
fs=100; % sampling freq.
Ts=1/fs; % sampling period
fr=0.03; % initial frequency resolution
T=floor(2/fr); % corresponding epoch length
[x,t]=signalA(Ts,T); % see Appendix ?
% default Welch
[Sw,fw]=pwelch(x,[],[],[],fs);
% Overrule defaults
fr=1; % set frequency resolution
M=ceil(2/fr/Ts); % then M is the number of samples
[Swx,fwx]=pwelch(x,M,[],[],fs); % .. of Hamming window

%% MATLAB CODE 1.5:Spectogram
T=20;
fs=100;
Ts=1/fs;
t=0:Ts:T;
x=cos(2*pi*(5+15*stepfun(t,10)).*t);
spectrogram(x,[],[],[],fs) % default
% with a window of 50 or 200 samples:
spectrogram(x,50,[],[],fs) % window of 0.5 seconds
spectrogram(x,200,[],[],fs) % window of 2.0 seconds
xlabel('Frequency (Hz)')
ylabel('Time [sec]')
colorbar

%% MATLAB CODE 1.6: Hilbert
t=0:.1:100;
x=t.^2.*(50-t).*exp(-t/8).*cos(t);
z=hilbert(x); % z = x + ix˜
A=abs(z); % amplitude
phi=phase(z); % phase
plot(t,x,t,A);
xlabel('Time (sec)')
ylabel('Amplitude (A)')
plot(t,phi,t,t);
xlabel('Time (sec)')
ylabel('Phase (phi)')


