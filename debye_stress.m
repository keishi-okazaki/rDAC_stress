function [fitresult, gof] = debye_stress(XX, YY, z) 
% XX: �f�o�C�����O�̐؂�J�����p�x�̈ʒu�AYY: d�l�̃s�[�N�ʒu�Az: ���u�̌X��, rDAC��z�x�Œ�ɂ��Ă���B��ver�ł�z��U��Ɖ������肵�Ȃ��B
[xData, yData] = prepareCurveData(XX, YY);

ft = fittype( 'D0*(1 - T/2*cos(z/180*pi)*sin(2*x/180*pi)+U/6*(1-3*(cos(z/180*pi))^2*(cos(x/180*pi))^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf z];
opts.StartPoint = [0.084 0.02 0.02 z];
opts.Upper = [Inf Inf Inf z];
% �t�B�b�e�B���O�����܂������Ȃ��ꍇ�͏����l�̕ύX�̕K�v����B

[fitresult, gof] = fit( xData, yData, ft, opts );
