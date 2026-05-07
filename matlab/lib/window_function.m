function data = window_function(data)
%WINDOW_FUNCTION 窗函数
%   data 信号

narginchk(1, 1);

[~, col] = size(data);
ham = hamming(length(data));
if col == 1
    data = data .* ham;
else
    data = data .* ham.';
end

