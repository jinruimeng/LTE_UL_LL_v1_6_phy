function myPlot3(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 07-Oct-2019 10:32:44 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ʹ�� plot �ľ������봴������
plot1 = plot(YMatrix1,'LineWidth',1);
set(plot1(1),'DisplayName','PCA','Marker','*','Color',[1 0 0]);
set(plot1(2),'DisplayName','С���任','Marker','^','Color',[1 0 1]);

% ���� ylabel
ylabel('�����ۼƹ�������');

% ���� xlabel
xlabel('������������');

% ȡ�������е�ע���Ա����������� X ��Χ
% xlim(axes1,[1 20]);
box(axes1,'on');
% ������������������
set(axes1,'XMinorGrid','on','YGrid','on');
% ���� legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');
