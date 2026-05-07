function res = write_complex_binary(data, filename)
%WRITE_COMPLEX_BINARY 把复数写进二进制文件
%   Format is interleaved float IQ
%   I,Q 32-bit float IQIQIQ ...
%RETURN 是否成功写入

narginchk(2, 2);

f = fopen(filename, 'wb');
if f < 0
    res = 0;
else
    re = real(data)';
    im = imag(data)';
    re = re(:)';
    im = im(:)';
    y = [re; im];
    y = y(:);
    res = fwrite(f, y, 'float');
    fclose(f);
end

end

