function [fitresult, gof] = pVoigtFit(X, Y)

[xData, yData] = prepareCurveData( X, Y );

ft = fittype('A*(n*(2/pi)*(w./(4*(x-xc).^2 + w^2)) + (1 - n)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-(4*log(2)/w^2).*(x-xc).^2))', 'independent', 'x', 'dependent', 'y' );
% ft = fittype('A*(n*(2/pi)*(w./(4*(x-xc).^2 + w^2)) + (1 - n)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-(4*log(2)/w^2).*(x-xc).^2))', 'independent', 'x', 'dependent', 'y' );
% 普通の式？Type 1 pseudo-Voigt function 
% (https://www.lightstone.co.jp/origin/flist13.html)
% (https://www.originlab.com/doc/origin-help/psdvoigt1-fitfunc)

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0]; % A n w xc
opts.StartPoint = [1 0.2 0.02 mean(X)];
% xcの初期値だけ各データ平均のd値を使うようにしている。
% フィッティングがうまくいかない場合は初期値の変更の必要あり。[0.13 0.47 0.012 0.245]
opts.Upper = [Inf 1 Inf Inf];

% モデルをデータに近似。
[fitresult, gof] = fit( xData, yData, ft, opts );

% データの近似をプロット、メインスクリプトに移動させた。
% figure;
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', '新規近似 1', 'Location', 'NorthEast' );
% % ラベル Axes
% xlabel x
% ylabel yyy
% grid on
% hold on


