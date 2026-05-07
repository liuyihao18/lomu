function conjugate_de_chirps = conjugate_demodulate(data_chirps, BW, T, fs)
%CONJUGATE_DEMODULATE 共轭解调模块
%   chirps 待解调的信号（行向量）
%   BW 标准Chirp的带宽
%   T 标准Chirp的持续时间
%   fs 采样率

narginchk(4, 4);

%% 生成 Down-Chirp
k = BW / T;
t = 0 : 1 / fs : T - 1 / fs;
up_chirp = chirp(t, -BW / 2, T, BW / 2, "complex");
down_chirp = conj(up_chirp);

%% 标准解调与共轭解调
n_chirp = length(data_chirps) / length(down_chirp);
conjugate_de_chirps = zeros(length(data_chirps), 1);
for i = 1 : n_chirp
    data_chirp = data_chirps((i - 1) * length(down_chirp) + 1 : i * length(down_chirp));
    de_chirp = data_chirp .* down_chirp;
    [fft_f, fft_res] = compute_fft(de_chirp, fs);
    [~, idx] = max(abs(fft_res));
    f0 = BW / 2 + fft_f(idx);
    if f0 > BW / 2
        f0 = f0 - BW;
    end
    delta = floor((BW / 2 - f0) / k * fs);
    if delta < 0 || delta + 1 > length(down_chirp)
        delta = 0;
    end
    conjugate_chirp = [down_chirp(end - delta + 1 : end) down_chirp(1 : end - delta)];
    conjugate_de_chirps((i - 1) * length(down_chirp) + 1 : i * length(down_chirp)) = data_chirp .* conjugate_chirp;
end
end

