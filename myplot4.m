function myplot4(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 07-Oct-2019 22:38:05 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ʹ�� plot �ľ������봴������
plot1 = plot(YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'DisplayName','�н���PCA','Marker','v',...
    'Color',[0.87058824300766 0.490196079015732 0]);
set(plot1(2),'DisplayName','����Э����ռ��������','Marker','square','Color',[0 0 1]);
set(plot1(3),'DisplayName','����Э�����뱾','Marker','*','Color',[1 0 0]);

% ���� ylabel
ylabel('���ȫ��Կ����');

% ���� xlabel
xlabel('Э����ռ���������/�������ĸ���');

% ȡ�������е�ע���Ա����������� X ��Χ
% xlim(axes1,[1 20]);
box(axes1,'on');
% ������������������
set(axes1,'XGrid','on','YGrid','on');
% ���� legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');
