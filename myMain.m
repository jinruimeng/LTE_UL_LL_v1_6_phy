close all;
clear;
clc;

N_subframes = 800;%帧的数量
col = 2400;%一帧中导频的数量
break_num = 1200;%截断长度
time = 1;%轮数
time2 = 1;%在一个速度上保持的轮数
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
            str = types(k) + " " + sprintf("第 %d 轮,第 %d 次轮完成了！", i,j);
            disp(str);
        end
        if i == 1 && k == 1
            %判断Excel是否已经打开，若已打开，就在打开的Excel中进行操作，
            %否则就打开Excel
            try
                Excel=actxGetRunningServer('Excel.Application');
            catch
                Excel = actxserver('Excel.Application'); 
            end
            %set(Excel, 'Visible', 1); %设置Excel属性为可见,效果是弹出Excel窗口
            %返回Excel工作簿句柄
            Workbooks = Excel.Workbooks;
            %若测试文件存在，打开该测试文件，否则，新建一个工作簿，并保存，文件名为测试.Excel
            %删除空表
            if exist(path_a,'file')
                Workbook = invoke(Workbooks,'Open',path_a);
            else
                Workbook = invoke(Workbooks, 'Add'); 
                Workbook.SaveAs(path_a);
            end
            Sheets = Excel.ActiveWorkBook.Sheets; 
            %返回工作簿中有多少工作表数
            %Count = Excel.ActiveWorkbook.Sheets.Count;
            Sheets.Item(1).Delete;%删除第deletesheet页的工作表
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

sprintf("跑完了！");


