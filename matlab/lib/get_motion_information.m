function [freq, amp, minimax, period, snr] = get_motion_information(dis, window_t, freq_begin, freq_end)
%GET_MOTION_INFORMATION 根据轨迹获得运动信息
%   dis 轨迹
%   window_t 时间窗口长度
%   freq_begin 频率范围起点
%   freq_end 频率范围终点
%   snr 感知频率信噪比
%RETURN 运动频率、幅度、极差、周期

narginchk(4, 4);

[fft_f, fft_res] = compute_fft(dis, 1 / window_t, 1e6);
fft_f = fft_f(length(fft_f) / 2 + 1 : end);
fft_res = fft_res(length(fft_res) / 2 + 1 : end);
[~, locs] = findpeaks(abs(fft_res), 'SortStr', "descend");

fft_f = fft_f(locs);
fft_f = fft_f(fft_f >= freq_begin & fft_f <= freq_end);
if isempty(fft_f)
    freq = 0;
    amp = 0;
    minimax = 0;
    period = 0;
    return
end

freq = fft_f(1);
loc = locs(1);
fft_res_signal = abs(fft_res(loc));
fft_res_noise = mean(abs(fft_res));
snr = 10 * log10(fft_res_signal / fft_res_noise);

if 0.75 / freq / window_t > length(dis) / 2
    amp = 0;
    minimax = 0;
    period = 0;
    return
end
[pks, locs1] = findpeaks(dis, "MinPeakDistance", 0.75 / freq / window_t);
max_pk = mean(pks(3 : end - 2));
[pks, locs2] = findpeaks(-dis, "MinPeakDistance", 0.75 / freq / window_t);
min_pk = mean(pks(3 : end - 2));
locs1 = locs1(3 : end - 2);
locs2 = locs2(3 : end - 2);
if length(locs1) < length(locs2)
    locs2 = locs2(1 : length(locs1));
else
    locs1 = locs1(1 : length(locs2));
end
amp = max_pk + min_pk;
minimax = max(dis) - min(dis);
period = 2 * abs(mean(locs2 - locs1)) * window_t;

end

