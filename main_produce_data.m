%%���ݴ��� ��ɶ�ȡ������excel
close all;
clear;
clc;

N_subframes = 3;%֡������
col = 2400;%һ֡�е�Ƶ������
break_num = 600;%�ضϳ���
time = 15;%����
time2 = 1;%��һ���ٶ��ϱ��ֵ�����
SNR = 100;
path = 'E:\workspace\keyan\channelDataP.xlsx';
% types = ["PedA", "PedB", "PedBcorr", "VehA", "VehB", "TU", "RA", "HT"];
types = ["VehA", "VehB", "TU", "RA", "HT"];
% types = ["VehA"];

for k = 1:length(types)
    type = char(types(k));
    for i = 1:time
        user_speed = 5 + 5 * (i - 1);
        for j = 1:time2
            sheet = types(k) + num2str(i) + "_" + num2str(j);
            % xlswrite(path,myGetAbs(myGenChannel(N_subframes,type,user_speed,1,SNR), break_num), sheet);
            xlswrite(path,myNum2Cell(myGenChannel(N_subframes,type,user_speed,1,SNR), break_num), sheet);
            str = types(k) + " " + sprintf("�� %d ��,�� %d ��������ˣ�", i,j);
            disp(str);
        end
        if i == 1 && k == 1
            %�ж�Excel�Ƿ��Ѿ��򿪣����Ѵ򿪣����ڴ򿪵�Excel�н��в�����
            %����ʹ�Excel
            try
                Excel=actxGetRunningServer('Excel.Application');
            catch
                Excel = actxserver('Excel.Application'); 
            end
            %set(Excel, 'Visible', 1); %����Excel����Ϊ�ɼ�,Ч���ǵ���Excel����
            %����Excel���������
            Workbooks = Excel.Workbooks;
            %�������ļ����ڣ��򿪸ò����ļ��������½�һ���������������棬�ļ���Ϊ����.Excel
            if exist(path,'file')
                Workbook = invoke(Workbooks,'Open',path);
            else
                Workbook = invoke(Workbooks, 'Add'); 
                Workbook.SaveAs(path);
            end
            Sheets = Excel.ActiveWorkBook.Sheets; 
            %���ع��������ж��ٹ�������
            %Count = Excel.ActiveWorkbook.Sheets.Count;
            Sheets.Item(1).Delete;%ɾ����deletesheetҳ�Ĺ�����
            Excel.ActiveWorkbook.Save;
            Excel.ActiveWorkbook.Close;
            Excel.Quit;
            Excel.delete;
        end
    end
end

sprintf("�����ˣ�");