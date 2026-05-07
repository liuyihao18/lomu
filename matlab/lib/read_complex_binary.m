function data = read_complex_binary(filename, count)
%READ_COMPLEX_BINARY 读取复数二进制文件
%   Format is interleaved float IQ
%   I,Q 32-bit float IQIQIQ ...
%RETURN 读取到的复数

narginchk(1, 2);

if nargin < 2
    count = Inf;
end

f = fopen(filename, 'rb');
if f < 0
    data = 0;
else
    t = fread(f, [2, count], 'float');
    fclose(f);
    data = t(1, :) + 1j * t(2, :);
    [r, c] = size(data);
    data = reshape(data, c, r);
end

end

