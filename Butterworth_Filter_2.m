%Butterworth Analog LPF parameters
Wc = 1.08;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + 1i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - 1i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + 1i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - 1i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + 1i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - 1i*Wc*sin(pi/2 + pi/16+3*pi/8);

%Band Edge speifications
fp1 = 59.2e3;
fs1 = 63.2e3;
fs2 = 83.2e3;
fp2 = 87.2e3;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260e3;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 260e3);
plot(f,abs(H))


yline(0.15,'-r','Magnitude = 0.15');
yline(0.85,'-r','Magnitude = 0.85');
yline(1.15,'-r','Magnitude = 1.15');
xline(59.2e3,'-g','f = 59.2kHz');
xline(63.2e3,'-g','f = 63.2kHz');
xline(83.2e3,'-g','f = 83.2kHz');
xline(87.2e3,'-g','f = 87.2.2kHz');

xlim([20e3,130e3])
ylim([0,1.5])
xlabel('Frequency (in Hz)')
ylabel('Magnitude')
title('Magnitude Plot')
hold on
grid on 