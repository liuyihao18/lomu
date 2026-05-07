function [f_range, loss] = joint_estimate(y1, y2, t, f_begin, f_end, f_n, lambda1, lambda2)
%JOINT_ESTIMATE 联合估计
%   y1      函数值1 -> 相位
%   y2      函数值2 -> 幅度
%   t       时间
%   f_begin 搜索频率起点
%   f_end   搜索频率终点
%   f_n     搜索频率个数
%   lambda1 控制参数1
%   lambda2 控制参数2
%RETUAN 搜索的结果

narginchk(5, 8)
if nargin < 6
    f_n = 1e4;
end
if nargin < 7
    lambda1 = 1;
end
if nargin < 8
    lambda2 = 1;
end

f_range = linspace(f_begin, f_end, f_n + 1);
f_range = f_range(1 : end - 1);
loss = zeros(1, length(f_range));
for i = 1 : length(f_range)
    f = f_range(i);
    f_sin = sin(2 * pi * f * t);
    f_cos = cos(2 * pi * f * t);
    A = [f_cos * f_cos', f_sin * f_cos', sum(f_cos)
         f_sin * f_cos', f_sin * f_sin', sum(f_sin) 
         sum(f_cos),     sum(f_sin),     length(t)];
    b1 = [f_cos * y1'
          f_sin * y1'
          sum(y1)'];
    b2 = [f_cos * y2'
          f_sin * y2'
          sum(y2)'];
    x1 = A^(-1) * b1;
    x2 = A^(-1) * b2;
    loss(i) = ...
              lambda1 * norm(y1 - x1(1) * f_cos - x1(2) * f_sin - x1(3))^2 ...
            + lambda2 * norm(y2 - x2(1) * f_cos - x2(2) * f_sin - x2(3))^2;
end

end

