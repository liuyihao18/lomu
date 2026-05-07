function data = remove_outlier(data)
%REMOVE_OUTLIER 移除离群值
%   data 数据

narginchk(1, 1);

idx = find(isoutlier(data));
for i = idx
        cnt = 0;
        data(i) = 0;
        for j = i - 25 : i + 25
            if j > 0 && j <= length(data) && ~any(j == idx)
                cnt = cnt + 1;
                data(i) = data(i) + data(j);
            end
        end
        if cnt > 0
            data(i) = data(i) / cnt;
        end
end
data = data - mean(data);

end

