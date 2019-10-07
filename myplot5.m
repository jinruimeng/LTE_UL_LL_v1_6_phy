function myplot5(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 07-Oct-2019 22:39:54 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','有交互PCA','Marker','v',...
    'Color',[0.600000023841858 0.200000002980232 0]);
set(plot1(2),'DisplayName','小波变换','Marker','^','Color',[1 0 1]);
set(plot1(3),'DisplayName','基于协方差空间均匀量化','Marker','square','Color',[0 0 1]);
set(plot1(4),'DisplayName','基于协方差码本','Marker','*','Color',[1 0 0]);
set(plot1(5),'DisplayName','原始信号','LineStyle','--','Color',[0 0 0]);

% 创建 ylabel
ylabel('最大安全密钥速率');

% 创建 xlabel
xlabel('保留分量个数');

box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XGrid','on','YGrid','on');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');

