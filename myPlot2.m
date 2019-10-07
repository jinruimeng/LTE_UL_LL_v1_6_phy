function myPlot2(YMatrix1)
%CREATEFIGURE1(YMatrix1)
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 06-Oct-2019 12:45:26 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ʹ�� semilogy �ľ������봴������
semilogy1 = semilogy(YMatrix1,'LineWidth',1,'Parent',axes1);
set(semilogy1(1),'DisplayName','�޽���PCA','Marker','v','Color',[0 0 0]);
set(semilogy1(2),'DisplayName','С���任','Marker','^','Color',[1 0 1]);
set(semilogy1(3),'DisplayName','����Э����ռ��������','Marker','square',...
    'Color',[0 0 1]);
set(semilogy1(4),'DisplayName','����Э�����뱾','Marker','*','Color',[1 0 0]);

% ���� ylabel
ylabel('��Կ��һ����');

% ���� xlabel
xlabel('����������(bit)');

% ȡ�������е�ע���Ա����������� X ��Χ
% xlim(axes1,[1 13]);
% ȡ�������е�ע���Ա����������� Y ��Χ
% ylim(axes1,[0.008 0.303]);
box(axes1,'on');
% ������������������
set(axes1,'XGrid','on','XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'XTickLabel',...
    {'12','16','20','24','28','32','36','40','44','48','52','56','60'},'YGrid',...
    'on','YMinorTick','on','YScale','log');
% ���� legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');
