function [fft_f, fft_res] = compute_fft(data, fs, fft_n)
%COMPUTE_FFT 计算傅里叶变换的结果及其对应频率
%   如果指定采样率则输出频率，否则输出角频率
%RETURN 傅里叶变换对应频率和结果

narginchk(1, 4);

if nargin < 3
    fft_n = length(data);
end
    
fft_res = fft(data, fft_n);
fft_f = linspace(-pi, pi, fft_n + 1);
fft_f = fft_f(1 : fft_n);
if nargin >= 2
    fft_f = fft_f / (2 * pi) * fs;
end

fft_res = fftshift(fft_res);

end

