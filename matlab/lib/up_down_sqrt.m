function [h] = up_down_sqrt(h_up, h_down)
%UP_DOWN_SQRT 对上下边带开平方
%   h_up 上边带
%   h_down 下边带
%RETURN 开平方的结果

narginchk(2, 2);

h_h = h_up .* h_down;
h = sqrt(h_h);
for i = 2 : length(h_h)
    if abs(angle(-h(i) ./ h(i - 1))) < abs(angle(h(i) ./ h(i - 1)))
        h(i) = -h(i);
    end
end

end

