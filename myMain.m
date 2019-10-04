close all;
clear;
clc;

N_subframes = 800;%֡������
col = 2400;%һ֡�е�Ƶ������
break_num = 1200;%�ضϳ���
time = 1;%����
time2 = 1;%��һ���ٶ��ϱ��ֵ�����
SNR_vec = 10;
path_a = 'E:\workspace\keyan\channelDataP_a.xlsx';
path_b = 'E:\workspace\keyan\channelDataP_b.xlsx';
% types = ["PedA", "PedB", "PedBcorr", "VehA", "VehB", "TU", "RA", "HT"];
types = ["VehA", "VehB", "TU", "RA", "HT"];
% types = ["RA"];

for k = 1:length(types)
    type = char(types(k));
    for i = 1:time
        user_speed = 5 + 5 * (i - 1);
        for j = 1:time2
            sheet = types(k) +"_" +  num2str(i) + "_" + num2str(j);
            [channelData_a,channelData_b] = myGenChannel(N_subframes,type,user_speed,1,SNR_vec);
            [Real_a,Imag_a] = myGetRealAndImag(channelData_a, break_num);
            [Real_b,Imag_b] = myGetRealAndImag(channelData_b, break_num);
            xlswrite(path_a, Real_a, sheet + "_Real");
            xlswrite(path_a, Imag_a, sheet + "_Imag");
            xlswrite(path_b, Real_b, sheet + "_Real");
            xlswrite(path_b, Imag_b, sheet + "_Imag");
            %xlswrite(path_a,myGetAbs(channelData_a, break_num), sheet);
            %xlswrite(path_b,myGetAbs(channelData_b, break_num), sheet);
            %xlswrite(path_a,myNum2Cell(channelData_a, break_num), sheet);
            %xlswrite(path_b,myNum2Cell(channelData_b, break_num), sheet);
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
            %ɾ���ձ�
            if exist(path_a,'file')
                Workbook = invoke(Workbooks,'Open',path_a);
            else
                Workbook = invoke(Workbooks, 'Add'); 
                Workbook.SaveAs(path_a);
            end
            Sheets = Excel.ActiveWorkBook.Sheets; 
            %���ع��������ж��ٹ�������
            %Count = Excel.ActiveWorkbook.Sheets.Count;
            Sheets.Item(1).Delete;%ɾ����deletesheetҳ�Ĺ�����
            Excel.ActiveWorkbook.Save;
            Excel.ActiveWorkbook.Close;
            Excel.Quit;
            Excel.delete;
            
            try
                Excel=actxGetRunningServer('Excel.Application');
            catch
                Excel = actxserver('Excel.Application'); 
            end
            Workbooks2 = Excel.Workbooks;          
            if exist(path_b,'file')
                Workbook2 = invoke(Workbooks2,'Open',path_b);
            else
                Workbook2 = invoke(Workbooks2, 'Add'); 
                Workbook2.SaveAs(path_b);
            end
            Sheets = Excel.ActiveWorkBook.Sheets; 
            Sheets.Item(1).Delete;
            Excel.ActiveWorkbook.Save;
            Excel.ActiveWorkbook.Close;
            Excel.Quit;
            Excel.delete;
        end
    end
end

sprintf("�����ˣ�");


