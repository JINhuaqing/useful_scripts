% the signal
t = linspace(0, 1, 1000);
x = 1* sin(2*pi*10*t) + sin(2*pi*20*t);
%plot(t, x)

fmin = 1; %Hz
fmax = 40;
fvec = (linspace(fmin,fmax, 100)).'; % freq range in Hz
fs=1000;
[q, fsamples] = get_spectral(x, fs, 1000);
plot(fsamples, q);
pow2db(q) % pow to db 

