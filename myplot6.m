function myplot6(zdata1)
%CREATEFIGURE1(zdata1)
%  ZDATA1:  surface zdata

%  �� MATLAB �� 07-Oct-2019 22:50:16 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ���� surf
surf(zdata1,'Parent',axes1);

% ���� ylabel
ylabel('�������');

% ���� xlabel
xlabel('�������');

grid(axes1,'on');
% ���� colorbar
colorbar('peer',axes1);

