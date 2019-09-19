function myPlot(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 19-Sep-2019 13:22:10 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','聚类协方差','Marker','o','Color',[0 0 0]);
set(plot1(2),'DisplayName','无反馈PCA','Marker','x');
set(plot1(3),'DisplayName','小波变换','Marker','diamond',...
    'Color',[0 0.447058826684952 0.74117648601532]);
set(plot1(4),'DisplayName','协方差空间均匀量化','Marker','square');

% 创建 ylabel
ylabel('密钥一致率');

% 创建 xlabel
xlabel('密钥比特数');

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[24 125]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'FontName','黑体');
% 创建 legend
legend(axes1,'show');

