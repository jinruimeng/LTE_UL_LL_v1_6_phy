function myPlot2(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 06-Oct-2019 12:45:26 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 semilogy 的矩阵输入创建多行
semilogy1 = semilogy(YMatrix1,'LineWidth',1,'Parent',axes1);
set(semilogy1(1),'DisplayName','无交互PCA','Marker','v','Color',[0 0 0]);
set(semilogy1(2),'DisplayName','小波变换','Marker','^','Color',[1 0 1]);
set(semilogy1(3),'DisplayName','基于协方差空间均匀量化','Marker','square',...
    'Color',[0 0 1]);
set(semilogy1(4),'DisplayName','基于协方差码本','Marker','*','Color',[1 0 0]);

% 创建 ylabel
ylabel('密钥不一致率');

% 创建 xlabel
xlabel('量化比特数(bit)');

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[1 13]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[0.008 0.303]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XGrid','on','XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'XTickLabel',...
    {'12','16','20','24','28','32','36','40','44','48','52','56','60'},'YGrid',...
    'on','YMinorTick','on','YScale','log');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');
