function [up_hs_div, down_hs_div] = syn_transceivers(hs, up_hs, down_hs)
%SYN_TRANSCEIVERS 同步收发端
%   hs 激励信号
%   h_up 上边带信号
%   h_down 下边带信号

narginchk(3, 3);

up_hs_div = up_hs ./ hs;
down_hs_div = down_hs ./ hs;
end

