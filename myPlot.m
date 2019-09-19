function myPlot(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x ���ݵ�����
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 19-Sep-2019 13:22:10 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ʹ�� plot �ľ������봴������
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','����Э����','Marker','o','Color',[0 0 0]);
set(plot1(2),'DisplayName','�޷���PCA','Marker','x');
set(plot1(3),'DisplayName','С���任','Marker','diamond',...
    'Color',[0 0.447058826684952 0.74117648601532]);
set(plot1(4),'DisplayName','Э����ռ��������','Marker','square');

% ���� ylabel
ylabel('��Կһ����');

% ���� xlabel
xlabel('��Կ������');

% ȡ�������е�ע���Ա����������� X ��Χ
% xlim(axes1,[24 125]);
box(axes1,'on');
% ������������������
set(axes1,'FontName','����');
% ���� legend
legend(axes1,'show');

