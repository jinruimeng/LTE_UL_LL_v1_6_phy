function [output] = myNum2Cell(input,break_num)
    %MYNUM2CELL �˴���ʾ�йش˺�����ժҪ
    %   �˴���ʾ��ϸ˵��
    [m,n] = size(input);
    p = n / break_num;
    output=strings(m * p, break_num);
    for i = 1:m
        for j = 1:p
            for k = 1:break_num
                output((i - 1)* p + j,k) = num2str(input(i,(j - 1) * break_num + k));
            end
        end
    end
    output = replace(output,'i','j'); 
    output = cellstr(output);
end
