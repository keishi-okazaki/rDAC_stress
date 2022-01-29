function [fitresult, gof] = pVoigtFit(X, Y)

[xData, yData] = prepareCurveData( X, Y );

ft = fittype('A*(n*(2/pi)*(w./(4*(x-xc).^2 + w^2)) + (1 - n)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-(4*log(2)/w^2).*(x-xc).^2))', 'independent', 'x', 'dependent', 'y' );
% ft = fittype('A*(n*(2/pi)*(w./(4*(x-xc).^2 + w^2)) + (1 - n)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-(4*log(2)/w^2).*(x-xc).^2))', 'independent', 'x', 'dependent', 'y' );
% ���ʂ̎��HType 1 pseudo-Voigt function 
% (https://www.lightstone.co.jp/origin/flist13.html)
% (https://www.originlab.com/doc/origin-help/psdvoigt1-fitfunc)

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0]; % A n w xc
opts.StartPoint = [1 0.2 0.02 mean(X)];
% xc�̏����l�����e�f�[�^���ς�d�l���g���悤�ɂ��Ă���B
% �t�B�b�e�B���O�����܂������Ȃ��ꍇ�͏����l�̕ύX�̕K�v����B[0.13 0.47 0.012 0.245]
opts.Upper = [Inf 1 Inf Inf];

% ���f�����f�[�^�ɋߎ��B
[fitresult, gof] = fit( xData, yData, ft, opts );

% �f�[�^�̋ߎ����v���b�g�A���C���X�N���v�g�Ɉړ��������B
% figure;
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', '�V�K�ߎ� 1', 'Location', 'NorthEast' );
% % ���x�� Axes
% xlabel x
% ylabel yyy
% grid on
% hold on


