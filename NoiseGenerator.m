function [Inoise] = NoiseGenerator(T,dt)
filterfrequency = 100; %Hz, frequencies to filter
% make the time and frequency axis
Tms = T*1000;
t = 0:dt:Tms;
dt_ins = dt/1000;
df = 1/(T+dt_ins);% freq resolution df = 1/(T+dt_ins)
fidx = 1:length(t)/2;% it has to be N/2 pts, where N = length(t) 
faxis = (fidx-1)*df;

%make the phases
Rr = randn(size(fidx));% ~N(0,1) over [-1,1]
distribphases = exp(1i*pi*Rr);% normal distributed phases on the unit circle
%make the amplitudes - filtered
%filterf = sqrt(1./(1+faxis.^2/filterfrequency^2));
filterf = sqrt(1./((2*pi*filterfrequency)^2+(2*pi*faxis).^2)); % see the PSD of an OU process, 
fourierA = distribphases.*filterf; % representation in fourier domain

% make it conj-symmetric so the ifft is real
fourierB = fliplr(conj(fourierA)); 
nss = [0,fourierA,fourierB];
signal = ifft(nss);
if ~isreal(signal)
    disp('trouble');
end
Inoise = signal;
scaling = std(Inoise);
Inoise = Inoise./scaling;
%trsf = fft(Inoise);
%PS = trsf.*conj(trsf);
%theoPS = 1./((2*pi*filterfrequency)^2+(2*pi*faxis).^2)*(1/scaling)^2;
%loglog(faxis,PS(1:length(faxis)),'-r',faxis,theoPS,'-k');

return% of function