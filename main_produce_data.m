%%数据处理 完成读取数据至excel
close all;
clear;
clc;

N_subframes = 3;%帧的数量
col = 2400;%一帧中导频的数量
break_num = 600;%截断长度
time = 15;%轮数
time2 = 1;%在一个速度上保持的轮数
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
            if exist(path,'file')
                Workbook = invoke(Workbooks,'Open',path);
            else
                Workbook = invoke(Workbooks, 'Add'); 
                Workbook.SaveAs(path);
            end
            Sheets = Excel.ActiveWorkBook.Sheets; 
            %返回工作簿中有多少工作表数
            %Count = Excel.ActiveWorkbook.Sheets.Count;
            Sheets.Item(1).Delete;%删除第deletesheet页的工作表
            Excel.ActiveWorkbook.Save;
            Excel.ActiveWorkbook.Close;
            Excel.Quit;
            Excel.delete;
        end
    end
end

sprintf("跑完了！");