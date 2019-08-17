function createfigure(X1, YMatrix1, X2, YMatrix2, X3, YMatrix3, X4, YMatrix4, X5, YMatrix5, X6, YMatrix6, X7, YMatrix7)
%CREATEFIGURE(X1, YMATRIX1, X2, YMATRIX2, X3, YMATRIX3, X4, YMATRIX4, X5, YMATRIX5, X6, YMATRIX6, X7, YMATRIX7)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data
%  X3:  vector of x data
%  YMATRIX3:  matrix of y data
%  X4:  vector of x data
%  YMATRIX4:  matrix of y data
%  X5:  vector of x data
%  YMATRIX5:  matrix of y data
%  X6:  vector of x data
%  YMATRIX6:  matrix of y data
%  X7:  vector of x data
%  YMATRIX7:  matrix of y data

%  Auto-generated by MATLAB on 16-Oct-2015 13:15:56

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YMinorTick','on','YScale','log');
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 40]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.0001 1]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'on');

% Create ylabel
ylabel('BER');

% Create xlabel
xlabel('SNR [dB]');

% Create title
title('UE BER');

% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(X1,YMatrix1,'MarkerSize',5,'Marker','o');
set(semilogy1(1),'DisplayName','UE 1, coded');
set(semilogy1(2),'DisplayName','UE 1, uncoded');

% Create multiple lines using matrix input to semilogy
semilogy2 = semilogy(X2,YMatrix2);
set(semilogy2(1),'Color',[0 0.447 0.741]);
set(semilogy2(2),'Color',[0.85 0.325 0.098]);

% Create multiple lines using matrix input to semilogy
semilogy3 = semilogy(X3,YMatrix3);
set(semilogy3(1),'Color',[0 0.447 0.741]);
set(semilogy3(2),'Color',[0.85 0.325 0.098]);

% Create multiple lines using matrix input to semilogy
semilogy4 = semilogy(X4,YMatrix4);
set(semilogy4(1),'Color',[0 0.447 0.741]);
set(semilogy4(2),'Color',[0.85 0.325 0.098]);

% Create multiple lines using matrix input to semilogy
semilogy5 = semilogy(X5,YMatrix5);
set(semilogy5(1),'Color',[0 0.447 0.741]);
set(semilogy5(2),'Color',[0.85 0.325 0.098]);

% Create multiple lines using matrix input to semilogy
semilogy6 = semilogy(X6,YMatrix6);
set(semilogy6(1),'Color',[0 0.447 0.741]);
set(semilogy6(2),'Color',[0.85 0.325 0.098]);

% Create multiple lines using matrix input to semilogy
semilogy7 = semilogy(X7,YMatrix7);
set(semilogy7(1),'Color',[0 0.447 0.741]);
set(semilogy7(2),'Color',[0.85 0.325 0.098]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','best','FontSize',9);

