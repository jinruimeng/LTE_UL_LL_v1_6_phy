function myplot4(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 07-Oct-2019 22:38:05 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','有交互PCA','Marker','v',...
    'Color',[0.87058824300766 0.490196079015732 0]);
set(plot1(2),'DisplayName','基于协方差空间均匀量化','Marker','square','Color',[0 0 1]);
set(plot1(3),'DisplayName','基于协方差码本','Marker','*','Color',[1 0 0]);

% 创建 ylabel
ylabel('最大安全密钥速率');

% 创建 xlabel
xlabel('协方差空间量化阶数/聚类中心个数');

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[1 20]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XGrid','on','YGrid','on');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');
