function [up_freq, down_freq, up_h, down_h, up_idx, down_idx] = get_up_down(fft_res, fs, tag, threshold, baseband_freq)
%GET_UP_DOWN 获取上下边带的信息
%   fft_res 数据傅里叶变换的结果
%   fs 采样率
%   tag 标签频移的参考频率
%   threshold 标签频移频率的误差范围
%   baseband_freq 基带频率
%RETURN 上下边带的信息

narginchk(4, 5);

if nargin < 5
    baseband_freq = 0;
end
fft_n = length(fft_res);
fft_f = linspace(-fs / 2, fs / 2, fft_n + 1);
fft_f = fft_f(1 : fft_n);

up_idx_ref = floor(fft_n / 2) + (tag + baseband_freq) / fs * fft_n;
down_idx_ref = floor(fft_n / 2) - (tag - baseband_freq) / fs * fft_n + 2;
up_idx_range = up_idx_ref - threshold / fs * fft_n + 1 : up_idx_ref + threshold / fs * fft_n;
down_idx_range = down_idx_ref - threshold / fs * fft_n + 1 : down_idx_ref + threshold / fs * fft_n;
up_idx_range = floor(up_idx_range);
down_idx_range = ceil(down_idx_range);

if (any(up_idx_range <= 0 | down_idx_range < 0 | up_idx_range > fft_n | down_idx_range > fft_n))
    disp("Baseband frequncy error!")
    baseband_freq = 0;
    up_idx_ref = floor(fft_n / 2) + (tag + baseband_freq) / fs * fft_n;
    down_idx_ref = floor(fft_n / 2) - (tag - baseband_freq) / fs * fft_n + 2;
    up_idx_range = up_idx_ref - threshold / fs * fft_n + 1 : up_idx_ref + threshold / fs * fft_n;
    down_idx_range = down_idx_ref - threshold / fs * fft_n + 1 : down_idx_ref + threshold / fs * fft_n;
    up_idx_range = floor(up_idx_range);
    down_idx_range = ceil(down_idx_range);
end

[~, up_idx] = max(abs(fft_res(up_idx_range)));
up_idx = up_idx + up_idx_range(1) - 1;
up_freq = fft_f(up_idx);
up_h = fft_res(up_idx);

[~, down_idx] = max(abs(fft_res(down_idx_range)));
down_idx = down_idx + down_idx_range(1) - 1;
down_freq = fft_f(down_idx);
down_h = fft_res(down_idx);

end

