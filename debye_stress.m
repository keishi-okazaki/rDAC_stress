function [fitresult, gof] = debye_stress(XX, YY, z) 
% XX: デバイリングの切り開いた角度の位置、YY: d値のピーク位置、z: 装置の傾き, rDACはz度固定にしている。現verではzを振ると解が安定しない。
[xData, yData] = prepareCurveData(XX, YY);

ft = fittype( 'D0*(1 - T/2*cos(z/180*pi)*sin(2*x/180*pi)+U/6*(1-3*(cos(z/180*pi))^2*(cos(x/180*pi))^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf z];
opts.StartPoint = [0.084 0.02 0.02 z];
opts.Upper = [Inf Inf Inf z];
% フィッティングがうまくいかない場合は初期値の変更の必要あり。

[fitresult, gof] = fit( xData, yData, ft, opts );
