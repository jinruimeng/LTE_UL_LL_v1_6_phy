function myPlot3(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 07-Oct-2019 10:32:44 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(YMatrix1,'LineWidth',1);
set(plot1(1),'DisplayName','PCA','Marker','*','Color',[1 0 0]);
set(plot1(2),'DisplayName','小波变换','Marker','^','Color',[1 0 1]);

% 创建 ylabel
ylabel('分量累计贡献能量');

% 创建 xlabel
xlabel('保留分量个数');

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[1 20]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XMinorGrid','on','YGrid','on');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');
