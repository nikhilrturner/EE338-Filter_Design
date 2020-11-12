f_samp = 260e3;

%Band Edge speifications
fs1 = 63.2e3;
fp1 = 59.2e3;
fp2 = 87.2e3;
fs2 = 83.2e3;

wp1 = fp1*2*pi/f_samp;
wp2 = fp2*2*pi/f_samp;
ws1 = fs1*2*pi/f_samp;
ws2 = fs2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*0.04*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 40;

%Ideal bandstop impulse response of length "n"
w_c1 = (wp1 + ws1)/2;
w_c2 = (wp2 + ws2)/2;
bs_ideal =  ideal_lp(pi,n) -ideal_lp(w_c2,n) + ideal_lp(w_c1,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))


yline(0.15,'-r','Magnitude = 0.15');
yline(0.85,'-r','Magnitude = 0.85');
yline(1.15,'-r','Magnitude = 1.15');
xline(59.2e3,'-g','f = 59.2kHz');
xline(63.2e3,'-g','f = 63.2kHz');
xline(83.2e3,'-g','f = 83.2kHz');
xline(87.2e3,'-g','f = 87.2.2kHz');

xlim([30e3,115e3])
ylim([0,1.5])
xlabel('Frequency (in Hz)')
ylabel('Magnitude')
title('Magnitude Plot')
hold on
grid on 