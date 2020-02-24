function ad = SteerVectorGenerate(ulaPos, sgnDoa)

% the ULA lies on the y-axis

c = 344;
f = 1000;
d = c/f/2; %half-wavelength
k = 2 * pi * f/c;


ad = exp(1j * k * ulaPos * [cos(sgnDoa), sin(sgnDoa), 0]');


end