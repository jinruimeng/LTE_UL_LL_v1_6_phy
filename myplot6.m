function myplot6(zdata1)
%CREATEFIGURE1(zdata1)
%  ZDATA1:  surface zdata

%  由 MATLAB 于 07-Oct-2019 22:50:16 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 surf
surf(zdata1,'Parent',axes1);

% 创建 ylabel
ylabel('分量序号');

% 创建 xlabel
xlabel('分量序号');

grid(axes1,'on');
% 创建 colorbar
colorbar('peer',axes1);

