function [rx,filtrado] = Fourier(var)

% var = var(~isnan(var));
x = fft(var);
[a1, k1] = max(abs(x(2:((length(x)/2)-1))'));
t_k = length(x)/k1;
n1i = k1 + 1;
n1f = length(x) - k1 + 1;
rx = zeros(length(x),1);
rx(n1f) = x(n1f);
rx(n1i) = x(n1i);
rx = real(ifft(rx));
filtrado = rx + nanmean(var); 

end