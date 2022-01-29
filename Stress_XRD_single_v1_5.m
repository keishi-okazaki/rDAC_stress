%% stress analysis from XRD Debye-rings for rDAC
% Keishi Okazaki, JAMSTEC, 6/21/2021
% 6/29/2021, ver 1.1, added a function on iterative calculations
% MATLAB & Curve Fitting Toolbox are required to execute this script
% 7/3/2021, ver 1.2, some improvements on the detrend procedure
% but now Signal Processing Toolbox is also required for using 'islocalmin' and 'findpeak' options
% 11/2/2021, ver. 1.3, slight modifications for debye_stress fitting
% 11/28/2021, ver. 1.4, slight modifictions for imporing data processes from csv file.
% 12/15/2021, ver. 1.5, data selection using R2 for debye_stress fitting

% Psuedo-Voigt function
% peak = A*(mu*(2/pi)*(w./(4*(x-xc).^2 + w^2)) + (1 - mu)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-(4*log(2)/w^2).*(x-xc).^2))

% Elastic distorsion of the Debye ring under the anisotropic stress condition
% d = d0*(1 - T/2*cos(z/180*pi)*sin(2*x/180*pi)+U/6*(1-3*(cos(z/180*pi))^2*(cos(x/180*pi))^2))
% T = tau/G, U = sigmaU/G, z = the tilt angle of the rDAC, currently 30 degree

clear
close all

%% �Ȃ񂩃p�����[�^�Ƃ�
tilt_angle = 30; % tilt angle of rDAC, in degree. 30 degree for 2020AB and 2021AB.
angle_rot = 0; % rotation of the Debye ring, in degree. ����0�x���f�t�H���g�BIPAnalizer����f�t�H���g�ŏ����o����-90����́H
detrend_opt = 'islocalmin'; % �f�[�^�̂Ńg�����h�����邩�ǂ����B'on', 'off', 'islocalmin', 'findpeak'�̂ǂꂩ
average_data = 'on'; % �f�[�^�̈ړ����ς���邩�ǂ����B'on' or 'off'
average_para = 3; % ���̈ړ����ς��Ƃ邩�B3���ƑO��1�_�̕��ς��Ƃ�B�s�[�N�̌`����Ώ̂��ƃs�[�N�������\��������̂Œ��ӁB
deg_use = [0 23 45 68 90 113 248 270 293 315 338]; % �C���|�[�g�����f�[�^�̂����g���p�x���w��B�����_�ȉ��͎l�̌ܓ�����Ă���B
%deg_use = [0 23 45 68 90 113 135 225 248 270 293 315 338]; % �C���|�[�g�����f�[�^�̂����g���p�x���w��B�����_�ȉ��͎l�̌ܓ�����Ă���B
threshold_R2 = 0.8; % �s�[�N�t�B�b�e�B���O��R2�����l�ȏ�̌��ʂ��g���ĉ��͌v�Z����B

% detrend_opt�̐ݒ�p�����[�^
ave_detrend = 5; % detrend'on'��'islocalmin'��'findpeak'�p�B �f�t�H���g��5���炢�������H�ŏ��ƍŌ��5�_���g���ăf�[�^��detrend
ave_detrend2 = 2; % 'findpeak'��'islocalmin'�p�B�f�t�H���g��2���炢�������H�ɏ��l+/-2�_���g���ăf�[�^��detrend
hashi = 25; % �ɏ��l��ROI�̒[����hashi�f�[�^�|�C���g�ȓ���T���B
mPP = 0.4; % �ɏ��l�̒Ⴓ�̒l��threshold, �傫���قǌ�����,1���炢�ł���������

% detrend options:
% 'off': �f�g�����h���Ȃ��B�o�b�N�O���E���h���Ȃ��Ȃ炱�����ł����͂��B���Ԃ�Ȃ��Ȃ�ĂȂ����ǁB
% 'on': �ŏ��ƍŌ��ave_detrend�œ��͂������̓_���g���ăf�[�^��detrend�Bmatlab��detrend funciton���͌����ڂ̓}�V�B
% �אڂ���s�[�N�Ƃ��܂�������ROI��ݒ�ł����ꍇ�͂���ł����B
% 'islocalmin','findpeak': 
% �ׂ̃s�[�N���ǂ����Ă������Ă��܂��ꍇ�̓s�[�N�ƃs�[�N�̊Ԃ̋ɏ��l���Ȃ�Ƃ������Ă��̈ʒu���o�b�N�O���E���h�̒l�Ƃ��� (�悤�ɓw�͂���)�B
% Signal processing toolbox���K�v�B�ɏ��l��������Ȃ������ꍇ��'on'�Ɠ������Ƃ�����B

summarytable =[]; % ��s�������Ă���
% summarytable.Properties.VariableNames = {'filename','Peak_d_value','tauG','errortauG','sigmaUG','errorsigmaUG'};

%% �v�Z���[�`��
renzoku = 0; % ���̒l��1�ɂȂ�܂ŌJ��Ԃ��B
d = 2; % �f�[�^�ė��p�����邩�ǂ����B2��ڈȍ~�̌v�Z�Ŏg�p�B

while renzoku < 1
    close all
    
    %�t�@�C�����J���A2��ڂ̌v�Z�œ���XRD�f�[�^���g�p����ۂ̓X�L�b�v�\�B
    if d == 2
        filename = uigetfile('*.xlsx;*.csv'); % IPAnalizer����ꊇ�����o�������t�@�C���Bcsv or xlsx
        [num,txt,raw] = xlsread(filename); % num: ���������Atxt: �e�L�X�g�f�[�^�����Araw: �S���̃f�[�^
        
        d0 = extractAfter(txt,' - '); % ' ? '���p�x�f�[�^�̈ʒu��������̂Ɏg�p���Ă���̂Ńt�@�C�����Ȃǂɂ͒���
        d0_char = d0(1,1:3:end); % �p�x�̃e�L�X�g�f�[�^
        d0_char2 = strrep(d0_char,'whole','999'); % "whole" ���_�~�[�p�x�ɒu���B
        diff_angle = str2double(d0_char) + angle_rot; % �e�L�X�g�f�[�^�œǂݍ��񂾊p�x�𐔒l�f�[�^�ɕϊ�
%         dataN = length(diff_angle); %�C�ӂ̊p�x�Ő؂蕪���ꂽ�f�[�^�Z�b�g�̐��B���˂���\�t�g�g�p�Bd�l�ŏo�͂��Ă��Ȃ��ƌv�Z�ł��Ȃ��B
        
        x_d_all = num(:,1:3:end); % �e�f�[�^��d�l�ix���j
        y_all = num(:,2:3:end); % �e�f�[�^��intensity�iy���j
        use_angle = ismember(diff_angle,deg_use);
        diff_angle = diff_angle(use_angle);
        dataN = length(diff_angle); %�C�ӂ̊p�x�Ő؂蕪���ꂽ�f�[�^�Z�b�g�̐��B���˂���\�t�g�g�p�Bd�l�ŏo�͂��Ă��Ȃ��ƌv�Z�ł��Ȃ��B
        
        x_d = x_d_all(:,use_angle); 
        y = y_all(:,use_angle);
        
    else
        % ����XRD���ĉ�͂���ꍇ�̓t�@�C�����J���������X�L�b�v����B
    end
    
    % ROI�I��
    f1 = figure;
    for i = 1:dataN
        figure(f1);
        plot(x_d(:,i),y(:,i))
        title('save the ROI as "ROI" in the workspace using the brush tool') % �S�p�x�Ńv���t�@�C�����S������悤�ɂ���B�ł��傫�������ǂ��Ȃ��B
        ylabel('Intensity');
        xlabel('d spacing [10^{-10} m]');
        hold on
    end
    
    pause % ROI��ۑ����enter�����čĊJ
    
    % ROI�I��pfigure window���J���Ă���������B
    try
        close(f1)
    catch
    end
    
    % ROI����f�[�^�؂�o��
    x_min_ROI = min(ROI(:,1));
    x_max_ROI = max(ROI(:,1));
    
    x_min = find(x_d(:,1)==x_min_ROI);
    x_max = find(x_d(:,1)==x_max_ROI);
    
    % �t�B�b�e�B���O�������f�[�^�����܂��؂�o���Ă��邩�m�F
    % for i = 1:dataN
    %     plot(x_d([x_min:x_max],i),y([x_min:x_max],i)+10*(i-1))
    %     hold on
    % end
    %
    % pause
    
    %% �؂蕪���f�[�^��������Ƀt�B�b�e�B���O
    f2 = figure('Position',[0 100 1250 750]);
    %     ave_detrend = 5; % detrend'on'��'findpeak'�p�B �f�t�H���g��5���炢�������H�ŏ��ƍŌ��5�_���g���ăf�[�^��detrend
    %     ave_detrend2 = 2; % 'findpeak'�p�B�f�t�H���g��2���炢�������H�ɏ��l+/-2�_���g���ăf�[�^��detrend
    %     hashi = 30; % ROI�̒[����hashi�|�C���g�ȓ���T���B
    %     mPP = 0.4; % �ɏ��l�̒Ⴓ�̒l��threshold, �傫���قǌ�����1���炢�ł���������
    
    for j = 1:dataN % dataN
        
        switch average_data  % �f�[�^�̈ړ�����
            case 'on'
                xx = movmean(x_d([x_min:x_max],j),average_para);
                yy = movmean(y([x_min:x_max],j),average_para);
                
            case 'off'
                xx = x_d([x_min:x_max],j);
                yy = y([x_min:x_max],j);
                
            otherwise
                error('choose an averaging option')
        end
        
        switch detrend_opt
            case 'on'
                % �ŏ��ƍŌ��ave_detrend�œ��͂������̓_���g���ăf�[�^��detrend�Bmatlab��dtrend funciton���͌����ڂ̓}�V�B
                % �אڂ���s�[�N�Ƃ��܂�������ROI��ݒ�ł����ꍇ�͂���ł����B
                ave_detrend = 5; % �f�t�H���g��5���炢�������H�ŏ��ƍŌ��5�_���g���ăf�[�^��detrend
                p = polyfit(xx([1:ave_detrend end-ave_detrend:end]),yy([1:ave_detrend end-ave_detrend:end]),1);
                yyy = yy - polyval(p,xx);
                
            case 'off' % �f�g�����h���Ȃ��B�o�b�N�O���E���h���Y��Ȃ炱�����ł����͂��B
                yyy = yy;
                
            case {'findpeak','islocalmin'}
                switch detrend_opt
                    case 'findpeak'
                        [pks,locs,width,hight] = findpeaks(-yy,'MinPeakProminence',mPP);
                        
                    case 'islocalmin'
                        [TF,PK] = islocalmin(yy,'MinProminence',mPP,'FlatSelection','center');
                        locs = find(TF==1);
                        hight = PK(locs);
                end
                
                if length(locs) == 2 &&  locs(1) <= hashi && (length(xx) - locs(2)) <= hashi % ���[�̍ŏ��l���^�悭���������ꍇ
                    p = polyfit(xx([locs(1)-ave_detrend2:locs(1)+ave_detrend2 locs(2)-ave_detrend2:locs(2)+ave_detrend2]),...
                        yy([locs(1)-ave_detrend2:locs(1)+ave_detrend2 locs(2)-ave_detrend2:locs(2)+ave_detrend2]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 1 % �ɏ��l�̗L���Ȃǂ̊m�F�p
                elseif length(locs) == 1 && locs <= hashi % left�s�[�N���������������ꍇ
                    p = polyfit(xx([locs(1)-2:locs(1)+ave_detrend2 end-ave_detrend:end]),...
                        yy([locs(1)-ave_detrend2:locs(1)+ave_detrend2 end-ave_detrend:end]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 2
                elseif length(locs) == 1 && (length(xx) - locs) <= hashi % right�s�[�N���������������ꍇ
                    p = polyfit(xx([1:ave_detrend locs(1)-ave_detrend2:locs(1)+ave_detrend2]),...
                        yy([1:ave_detrend locs(1)-ave_detrend2:locs(1)+ave_detrend2]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 3
                elseif length(locs) > 2
                    left_part = find(locs(1) <= hashi);
                    [left_h, left_pos] = max(hight(left_part));
                    right_part = find(locs >= (length(xx) - hashi));
                    [right_h, right_pos] = max(hight(right_part));
                    
                    if isempty(left_part(left_pos)) && isnumeric(right_part(right_pos))% left�ŏ��l��������Ȃ������ꍇ
                        p = polyfit(xx([1:ave_detrend locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),...
                            yy([1:ave_detrend locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),1);
%                         kek = 41
                    elseif isempty(right_part(right_pos)) && isnumeric(left_part(left_pos))% right�ŏ��l��������Ȃ������ꍇ
                        p = polyfit(xx([locs(left_part(left_pos))-2:locs(left_part(left_pos))+ave_detrend2 end-ave_detrend:end]),...
                            yy([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 end-ave_detrend:end]),1);
%                         kek = 42
                    else % �ŏ��l���������������ꍇ
                        p = polyfit(xx([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),...
                            yy([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),1);
%                         kek = 43
                    end
                    
                    yyy = yy - polyval(p,xx);
                    
                else % �ŏ��l��������Ȃ������ꍇ
                    p = polyfit(xx([1:ave_detrend end-ave_detrend:end]),yy([1:ave_detrend end-ave_detrend:end]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 5
                    
                end
                
            otherwise
                error('Choose a detrend optioin')
        end
        
        
        % psuedo-Voigt�֐��Ŋe�؂蕪���p�x�̃s�[�N���t�B�b�e�B���O
        [fitresult_peak, gof] = pVoigtFit(xx,yyy);
        fitpeak(j,:) = coeffvalues(fitresult_peak); %A n w xc
        % a�������An�����[�����c�֐��ƃK�E�X�֐��̔�Ax���s�[�N���Axc���s�[�N�ʒu
        fitpeak95(j,:) = reshape(confint(fitresult_peak),[1,8]);% 95%�M����ԁAa1 min,a1 max, b1 min, b1 max...
        R2_peak(j) = gof.rsquare; %R2�l���ۑ����Ă����B
        
        % �f�[�^�̋ߎ����v���b�g
        figure(f2);
        subplot(2,round(dataN/2),j)
        plot(fitresult_peak, x_d([x_min:x_max]), yyy);
        xlabel('d spaceing [10^{-10} m]')
        ylabel('Intensity')
        title(num2str(diff_angle(j)))
        grid on
        hold on
        set(gca,'FontName','Helvetica','FontSize',14);
        legend('off');
        
        yy = []; yyy = [];
    end
    
    %% ���͌v�Z�̂��߂�fitting
    rmdata = isoutlier(fitpeak(:,4)); % �O��l���폜�Av1.3�Œǉ�
    fitpeak_calc = fitpeak((not(rmdata) & [R2_peak>=threshold_R2]'),4); % �O��l���폜
    diff_angle_calc = diff_angle((not(rmdata) & [R2_peak>=threshold_R2]')); % �O��l���폜
    
    [fitresult_stress, gofstress] = debye_stress(diff_angle_calc,fitpeak_calc,tilt_angle); % D0, T = tau/G, U = sigmaU/G
    fitstress = coeffvalues(fitresult_stress); %D0 T U tilt_angle
    fitstress95(1,:) = reshape(confint(fitresult_stress),[1,8]);% 95%�M����ԁAa1 min,a1 max, b1 min, b1 max...
    R2_stress = gofstress.rsquare;
    
    %% plot ���쐬
    f3 = figure;
    plot(diff_angle_calc,fitpeak_calc,...
        'MarkerFaceColor',[1 1 1],'MarkerSize',10,...
        'Marker','o',...
        'LineWidth',2,...
        'LineStyle','none');
    hold on
    
    diff_angle_fitplot = linspace(0,360,361) + angle_rot;
    fitted_line = fitstress(1).*(1 - fitstress(2)./2.*cos(fitstress(4)./180*pi).*...
        sin(2.*diff_angle_fitplot./180*pi)+fitstress(3)./6.*(1-3.*(cos(fitstress(4)./180*pi)).^2.*...
        (cos(diff_angle_fitplot./180*pi)).^2));
    
    plot(diff_angle_fitplot,fitted_line,'-r')
    
    xlim([min(diff_angle_fitplot) max(diff_angle_fitplot)])
    ylabel('d spacing [10^{-10} m]');
    xlabel('Angle [deg]');
    
    box(gca,'on');
    set(gca,'FontName','Helvetica','FontSize',14);
    
    %% save results
    savefilename = uiputfile(strcat(extractBefore(filename,'.'),'.mat'));
    savefigname = strcat(extractBefore(savefilename,'.'),'.fig');
    savefig(f2,strcat(extractBefore(savefilename,'.'),'_angle.fig'));
    savefig(f3,savefigname)
    save(strcat(savefilename))
    
    % savefigname = extractBefore(filename,'.');
    % savefig(gcf,strcat('result_',savefigname))
    % save(strcat('result_',savefigname))
    
    %% Export tau/G and sigmaU/G with errors as .csv file
    try
        load summarytable.mat
        disp('fitting results are added and saved in summarytable.xlsx')
    catch
        disp('fitting results are newly saved in summarytable.xlsx')
    end
    
    summarytable_new = table(string(filename), mean(fitpeak_calc),fitstress(2),fitstress(2) - fitstress95(1,3),fitstress(3),fitstress(3) - fitstress95(1,5));
    summarytable_new.Properties.VariableNames = {'filename','Peak_d_value','tauG','errortauG','sigmaUG','errorsigmaUG'}
    summarytable = [summarytable;summarytable_new];
    
    save('summarytable.mat','summarytable')
    writetable(summarytable,'summarytable.xlsx')
    
    %% ������x��蒼�� or �ʉ�͂����邩�ǂ�������
    c = menu('More calculations?','Yes','No');
    renzoku = c - 1;
    
    if c == 1
        d = menu('Use same XRD data?','Yes','No, load new XRD data');
    else
    end
    
end
%% END
