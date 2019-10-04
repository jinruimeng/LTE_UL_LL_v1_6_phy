N = 1024;
M = 1600;
data_num = 2;
for data_index = 0:data_num
    fprintf('%d start!\n', data_index)
    eval(['data',num2str(data_index) ,'2=','data',num2str(data_index),'(:,1:N);']);
    for j = 1:N
        eval(['tmp =  mean(data',num2str(data_index),'2(:,j));']);
        eval(['data', num2str(data_index), '2(:,j) = data',  num2str(data_index), '2(:,j) - tmp;']);
    end
    eval(['C',num2str(data_index), '=cov(data', num2str(data_index),'2);']);
    eval(['power',num2str(data_index),'1 = 0;']);
    for i = 1:N
       eval(['power',num2str(data_index),'1 =  power',num2str(data_index),'1 + C',num2str(data_index),'(i,i);']);
    end
    eval(['power',num2str(data_index),'1 = power',num2str(data_index),'1 /N;']);
    eval(['s',num2str(data_index), '= svd(C', num2str(data_index), ');']);
    eval(['power', num2str(data_index), '2 = sum(s', num2str(data_index), ')/N;']);
    eval(['power',num2str(data_index),'3 = 0;']);
    for i = 1:M
        for j = 1 : N
            eval(['power', num2str(data_index), '3 = power', num2str(data_index), '3 + data', num2str(data_index), '2(i,j)^2;']);
        end
    end
    eval(['power', num2str(data_index), '3 = power', num2str(data_index), '3/M/N;']);
    eval(['npower',num2str(data_index), '= min(s', num2str(data_index), ');']);
    eval(['snr', num2str(data_index), '1 = (power', num2str(data_index), '3/npower', num2str(data_index), ') - 1;']);
    eval(['snr', num2str(data_index), '2 = 10*log10(snr', num2str(data_index), '1);']);
end

% data12 = data1(:,1:1024);
% for j = 1:1024
%     tmp = mean(data12(:,j));
%     data12(:,j) = data12(:,j) - tmp;
% end
% C1 = cov(data12);
% power11 = 0;
% for i = 1:1024
%    power11 =  power11 + C1(i,i);
% end
% power11 = power11 /1024;
% s1 = svd(C1);
% power12 = sum(s1)/1024;
% power13 = 0;
% for i = 1:1600
%     for j = 1 : 1024
%         power13 = power13 + data12(i,j)^2;
%     end
% end
% power13 = power13 / 1600 /1024;
% npower1 = min(s1);
% snr11 = (power13/npower1) - 1;
% snr12 = 10*log10(snr11);