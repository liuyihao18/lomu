%[text] # 在线处理数据
clear
clc
addpath('lib');
%%
%[text] ## 参数
% 读取配置文件
env = loadenv('.env');
data_dir = char(env("DATA_DIR"));
config_file = fopen([data_dir, '\config.json']);
if config_file == -1
    error("配置文件不存在！");
else
    config = jsondecode(fscanf(config_file, "%s"));
    fclose(config_file);
end

fs = config.fs;                     % 采样率
center_freq = config.center_freq;   % 载波中心频率
tag = config.tag(1 : 1);            % 标签频移参考频率
tag_num = length(tag);              % 标签数量
threshold = 2e3;                    % 标签频移频率的误差范围
fft_n = 1e5;                        % 傅里叶变换长度
baseband_freq_ref = 0;              % 基带频率

every_window_flag = true;           % 是否逐窗口寻找标签频率
%%
%[text] ## 信号
T = config.T;
BW = config.BW;
if isfield(config, 'f0')
    f0 = config.f0;
else
    f0 = -BW / 2;
end
k = BW / T;
t = 0 : 1 / fs : T - 1 / fs;
up_chirp = chirp(t, f0, T, f0 + BW, "complex");
down_chirp = conj(up_chirp);
write_complex_binary(up_chirp, [data_dir '\chirp.cf32']);
%%
%[text] ## 连接TCP
% 创建TCP服务器
server = tcpserver("127.0.0.1", 2000, "Timeout", 1);
% 设置接收参数
epoch = 0; % 轮次
begin_epoch = 10; % 忽略掉的轮次
read_data_length = 2e5; % 接收数据的长度
read_data_time = read_data_length / fs; % 接收数据对应的时间
% 设置历史记录
history_iq_period = 2.5;
history_dis_period = 10;
history_max_iq_cnt = history_iq_period / read_data_time;
history_max_dis_cnt = history_dis_period / read_data_time;
history_cnt = 0;
history_iq = cell(tag_num, 1);
history_time = cell(tag_num, 1);
history_dis = cell(tag_num, 1);
history_iq_abs_mean = zeros(tag_num, 1);
history_dis_mean = zeros(tag_num, 1);
history_freq = zeros(tag_num, 1);
history_amp = zeros(tag_num, 1);
history_gamma = 0.01;
for j = 1 : tag_num
    history_iq{j} = 1;
    history_time{j} = 0;
    history_dis{j} = 0;
end
dump_dis = cell(tag_num, 1);
for j = 1 : tag_num
    dump_dis{j} = 0;
end
motion_freq_ref = zeros(tag_num, 1);
% 设置窗口
f = figure("Name", "监视器");
if tag_num == 1
    % set(f, "Position", [560, 140, 800, 800]);
    set(f, "Position", [800, 100, 800, 800]);
elseif tag_num == 2
    set(f, "Position", [160, 140, 1600, 800]);
elseif tag_num == 3
    set(f, "Position", [160, 90, 1600, 900]);
end
% 连接
while ~server.Connected
    pause(0.01);
end
% 处理数据
while server.Connected
    raw = read(server, 2 * read_data_length, "single");
    epoch = epoch + 1;
    if epoch < begin_epoch
        continue;
    end
    if length(raw) < read_data_length
        continue;
    end
    data0 = raw(1 : 2 : end) + 1j * raw(2 : 2 : end);
    window_size = length(up_chirp);
    % 对齐起点与信号拼接
    if ~exist('chirp_start', "var")
        data = data0;
        pre = data(1 : min(length(data), floor(33 * T * fs)));
        [corr_res, ~] = xcorr(pre, up_chirp);
        corr_res = abs(corr_res(length(pre) + 1 : length(corr_res)));
        [~, chirp_starts] = findpeaks(corr_res, "MinPeakDistance", (T / 2) * fs);
        chirp_start = chirp_starts(floor(length(chirp_starts) / 2)) + 1;
        % chirp_start = 1;
    else
        if window_end < window_begin
            data = [data(1: end) data0];
        else
            data = [data(window_end - window_size + 1 : end) data0];
            chirp_start = 1;
        end
    end
    window_begin = chirp_start;
    window_end = length(data);
    n_window = floor((window_end - window_begin) / window_size);
    window_end = window_begin + n_window * window_size - 1;
    if n_window == 0
        continue;
    end
    window_t = window_size * 1 / fs;
    processed_data = data(window_begin : window_end);
    if isfield(config, "dc_offset") && config.dc_offset
        dc_offset = mean(processed_data);
        processed_data = processed_data - dc_offset;
    end
    % 共轭解调
    processed_data = conjugate_demodulate(processed_data, BW, T, fs);
    up_freq = zeros(tag_num, 1);
    down_freq = zeros(tag_num, 1);
    up_idx = zeros(tag_num, 1);
    down_idx = zeros(tag_num, 1);
    freqs = zeros(1, n_window);
    up_freqs = zeros(tag_num, n_window);
    down_freqs = zeros(tag_num, n_window);
    hs = zeros(1, n_window);
    up_hs = zeros(tag_num, n_window);
    down_hs = zeros(tag_num, n_window);
    [fft_f, fft_res] = compute_fft(processed_data, fs, fft_n);
    [pk, ~] = max(abs(fft_res));
    [~, locs] = findpeaks(abs(fft_res), "MinPeakHeight", 0.3 * pk);
    baseband_freqs = fft_f(locs);
    [~, idx] = min(abs(baseband_freqs - baseband_freq_ref));
    baseband_freq = baseband_freqs(idx);
    for j = 1 : tag_num
        [up_freq(j), down_freq(j), ~, ~, up_idx(j), down_idx(j)] = get_up_down(fft_res, fs, tag(j), threshold, baseband_freq);
    end
    for i = 1 : n_window
        window_data = processed_data((i - 1) * window_size + 1 : i * window_size);
        % 窗函数
        window_data = window_function(window_data);
        [fft_f, fft_res] = compute_fft(window_data, fs, fft_n);
        if every_window_flag
            [pk, ~] = max(abs(fft_res));
            [~, locs] = findpeaks(abs(fft_res), "MinPeakHeight", 0.3 * pk);
            baseband_freqs = fft_f(locs);
            [~, idx] = min(abs(baseband_freqs - baseband_freq_ref));
            idx = locs(idx);
            baseband_freq = fft_f(idx);
            for j = 1 : tag_num
                [~, ~, ~, ~, up_idx(j), down_idx(j)] = get_up_down(fft_res, fs, tag(j), threshold, baseband_freq);
            end
        end
        % 提取信号
        freqs(i) = fft_f(idx);
        hs(i) = fft_res(idx);
        up_freqs(:, i) = fft_f(up_idx);
        down_freqs(:, i) = fft_f(down_idx);
        up_hs(:, i) = fft_res(up_idx);
        down_hs(:, i) = fft_res(down_idx);
    end
    % 修正
    if every_window_flag
        % 持续跟踪消除傅里叶变换跳变算法
        for i = 2 : n_window
            hs(i) = hs(i) .* exp(1j * pi * (freqs(i) - freqs(1)) * window_t);
            up_hs(:, i) = up_hs(:, i) .* exp(1j * pi * (up_freqs(:, i) - up_freqs(:, 1)) * window_t);
            down_hs(:, i) = down_hs(:, i) .* exp(1j * pi * (down_freqs(:, i) - down_freqs(:, 1)) * window_t);
        end
    end
    history_cnt = history_cnt + 1;
    for j = 1 : tag_num
        % 收发端偏移消除
        [up_hs_div, down_hs_div] = syn_transceivers(hs, up_hs(j, :), down_hs(j, :));
        % 节点频率漂移消除 + 持续跟踪消歧义算法
        up_down_sqrt_hs = up_down_sqrt(up_hs_div, down_hs_div);
        up_down_sqrt_hs = up_down_sqrt_hs ./ exp(1j * angle(up_down_sqrt_hs(1)));
        time = 0 : window_t : (length(up_down_sqrt_hs) - 1) * window_t;
        lambda = 3e8 / center_freq * 1000; % 单位（mm）
        dis = unwrap(angle(up_down_sqrt_hs)) ./ (2 * pi) * lambda / 2;
        if history_cnt < history_max_iq_cnt
            history_iq{j} = [history_iq{j}(1 : end - 1) exp(1j * angle(history_iq{j}(end))) * up_down_sqrt_hs];
        else
            history_iq{j} = [history_iq{j}(length(up_down_sqrt_hs) : end - 1) exp(1j * angle(history_iq{j}(end))) * up_down_sqrt_hs];
        end
        if history_cnt < history_max_dis_cnt
            history_dis{j} = [history_dis{j}(1 : end - 1) history_dis{j}(end) + dis];
            history_time{j} = [history_time{j}(1 : end - 1) history_time{j}(end) + time];
        else
            history_dis{j} = [history_dis{j}(length(dis) : end - 1) history_dis{j}(end) + dis];
            history_time{j} = [history_time{j}(length(time) : end - 1) history_time{j}(end) + time];
        end
        history_iq_abs_mean(j) = history_iq_abs_mean(j) * history_gamma + (1 - history_gamma) * mean(abs(history_iq{j}));
        history_dis_weight = (1 : length(history_dis{j})) .^ 2;
        history_dis_weight = history_dis_weight' ./ sum(history_dis_weight);
        history_dis_mean(j) = history_dis_mean(j) * history_gamma + (1 - history_gamma) * history_dis{j} * history_dis_weight;
        dump_dis{j} = [dump_dis{j}(1 : end - 1) dump_dis{j}(end) + dis];
        figure(f)
        subplot(3, tag_num, j)
        scatter(real(history_iq{j}), imag(history_iq{j}), linspace(1, 10, length(history_iq{j})), linspace(0.6, 0.3, length(history_iq{j})), 'filled')
        % xlim(history_iq_abs_mean(j) * [-3 3])
        ylim(history_iq_abs_mean(j) * [-1.25 1.25])
        axis equal
        xlabel("I")
        ylabel("Q")
        title(['Tag ' num2str(j)])
        subplot(3, tag_num, j + tag_num)
        % 各种滤波器
        %curr_dis = my_filter3(history_dis{j}, window_t, diag([1e-2 1e-1 1e-0]), 1);
        %curr_dis = curr_dis(1, :);
        curr_dis = history_dis{j};
        %curr_dis = movmean(history_dis{j}, floor(5 * (1 + epoch/1000)));
        plot(history_time{j}, curr_dis(1, :))
        xlim([history_time{j}(1) history_time{j}(1) + history_dis_period])
        ylim(history_dis_mean(j) + [-150 150])
        xlabel("时间（s）")
        ylabel("距离（mm）")
        title(['Tag ' num2str(j)])
        % 运动信息提取
        [freq, amp, ~, ~] = get_motion_information(curr_dis(floor(length(curr_dis) * 0.618) : end), window_t, 0.01, 2.0);
        if (epoch - begin_epoch + 1) > 3
            history_freq(j) = history_freq(j) * history_gamma + (1 - history_gamma) * freq;
            if ~isnan(amp)
                history_amp(j) = history_amp(j) * history_gamma + (1 - history_gamma) * amp;
            end
        end
    end
    if exist('txt1', "var")
        delete(txt1);
    end
    % if exist('txt2', "var")
    %     delete(txt2);
    % end
    for j = 1 : tag_num
        subplot(3, tag_num, j + 2 * tag_num)
        if history_amp(j) > 3
            txt1(j) = text(0.3, 0.5, ['频率：' num2str(history_freq(j)) 'Hz'], "FontSize", 16);
            % txt2(j) = text(0.3, 0.4, ['振幅：' num2str(history_amp(j)) 'mm'], "FontSize", 16);
        else
            txt1(j) = text(0.3, 0.5, '频率：0Hz', "FontSize", 16);
            % txt2(j) = text(0.3, 0.4, '振幅：0mm', "FontSize", 16);
        end
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40.5}
%---
