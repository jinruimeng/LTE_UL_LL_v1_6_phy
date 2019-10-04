function [Real,Imag] = myGetRealAndImag(input,break_num)
    %MYNUM2CELL 此处显示有关此函数的摘要
    %   此处显示详细说明
    [m,n] = size(input);
    p = n / break_num;
    Real = zeros(m * p, break_num);
    Imag = zeros(m * p, break_num);
    for i = 1:m
        for j = 1:p
            for k = 1:break_num
                Real((i - 1)* p + j,k) = real(input(i,(j - 1) * break_num + k));
                Imag((i - 1)* p + j,k) = imag(input(i,(j - 1) * break_num + k));
            end
        end
    end
end