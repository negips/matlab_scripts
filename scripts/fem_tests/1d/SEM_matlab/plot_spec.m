function [freq,ff,h]=plot_spec(f,dt)

nfft = 2^nextpow2(length(f));
ff = fft(f(:),nfft);
ff = abs(ff(2:nfft/2+1));
freq = [1:nfft/2]' /(nfft*dt) ;
h=loglog( freq, ff );
xlabel('Frequency')
