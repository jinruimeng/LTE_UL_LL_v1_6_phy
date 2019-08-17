function [output] = myGetAbs(input,break_num)
    %MYNUM2CELL �˴���ʾ�йش˺�����ժҪ
    %   �˴���ʾ��ϸ˵��
    [m,n] = size(input);
    p = n / break_num;
    output = zeros(m * p, break_num);
    for i = 1:m
        for j = 1:p
            for k = 1:break_num
                output((i - 1)* p + j,k) = abs(input(i,(j - 1) * break_num + k));
            end
        end
    end
end

