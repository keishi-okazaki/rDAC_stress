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

%% なんかパラメータとか
tilt_angle = 30; % tilt angle of rDAC, in degree. 30 degree for 2020AB and 2021AB.
angle_rot = 0; % rotation of the Debye ring, in degree. 下が0度がデフォルト。IPAnalizerからデフォルトで書き出すと-90を入力？
detrend_opt = 'islocalmin'; % データのでトレンドをするかどうか。'on', 'off', 'islocalmin', 'findpeak'のどれか
average_data = 'on'; % データの移動平均を取るかどうか。'on' or 'off'
average_para = 3; % 何個の移動平均をとるか。3だと前後1点の平均をとる。ピークの形が非対称だとピークがずれる可能性もあるので注意。
deg_use = [0 23 45 68 90 113 248 270 293 315 338]; % インポートしたデータのうち使う角度を指定。小数点以下は四捨五入されている。
%deg_use = [0 23 45 68 90 113 135 225 248 270 293 315 338]; % インポートしたデータのうち使う角度を指定。小数点以下は四捨五入されている。
threshold_R2 = 0.8; % ピークフィッティングのR2が一定値以上の結果を使って応力計算する。

% detrend_optの設定パラメータ
ave_detrend = 5; % detrend'on'と'islocalmin'と'findpeak'用。 デフォルトは5くらいがいい？最初と最後の5点を使ってデータをdetrend
ave_detrend2 = 2; % 'findpeak'と'islocalmin'用。デフォルトは2くらいがいい？極小値+/-2点を使ってデータをdetrend
hashi = 25; % 極小値をROIの端からhashiデータポイント以内を探す。
mPP = 0.4; % 極小値の低さの値のthreshold, 大きいほど厳しい,1くらいでもいいかも

% detrend options:
% 'off': デトレンドしない。バックグラウンドがないならこっちでいいはず。たぶんないなんてないけど。
% 'on': 最初と最後のave_detrendで入力した数の点を使ってデータをdetrend。matlabのdetrend funcitonよりは見た目はマシ。
% 隣接するピークとうまく分けてROIを設定できた場合はこれでいい。
% 'islocalmin','findpeak': 
% 隣のピークがどうしても入ってしまう場合はピークとピークの間の極小値をなんとか見つけてその位置をバックグラウンドの値とする (ように努力する)。
% Signal processing toolboxが必要。極小値が見つからなかった場合は'on'と同じことをする。

summarytable =[]; % 空行列を作っておく
% summarytable.Properties.VariableNames = {'filename','Peak_d_value','tauG','errortauG','sigmaUG','errorsigmaUG'};

%% 計算ルーチン
renzoku = 0; % この値が1になるまで繰り返す。
d = 2; % データ再利用をするかどうか。2回目以降の計算で使用。

while renzoku < 1
    close all
    
    %ファイルを開く、2回目の計算で同じXRDデータを使用する際はスキップ可能。
    if d == 2
        filename = uigetfile('*.xlsx;*.csv'); % IPAnalizerから一括書き出ししたファイル。csv or xlsx
        [num,txt,raw] = xlsread(filename); % num: 数字だけ、txt: テキストデータだけ、raw: 全部のデータ
        
        d0 = extractAfter(txt,' - '); % ' ? 'を角度データの位置を見つけるのに使用しているのでファイル名などには注意
        d0_char = d0(1,1:3:end); % 角度のテキストデータ
        d0_char2 = strrep(d0_char,'whole','999'); % "whole" をダミー角度に置換。
        diff_angle = str2double(d0_char) + angle_rot; % テキストデータで読み込んだ角度を数値データに変換
%         dataN = length(diff_angle); %任意の角度で切り分けれたデータセットの数。瀬戸さんソフト使用。d値で出力していないと計算できない。
        
        x_d_all = num(:,1:3:end); % 各データのd値（x軸）
        y_all = num(:,2:3:end); % 各データのintensity（y軸）
        use_angle = ismember(diff_angle,deg_use);
        diff_angle = diff_angle(use_angle);
        dataN = length(diff_angle); %任意の角度で切り分けれたデータセットの数。瀬戸さんソフト使用。d値で出力していないと計算できない。
        
        x_d = x_d_all(:,use_angle); 
        y = y_all(:,use_angle);
        
    else
        % 同じXRDを再解析する場合はファイルを開く部分をスキップする。
    end
    
    % ROI選択
    f1 = figure;
    for i = 1:dataN
        figure(f1);
        plot(x_d(:,i),y(:,i))
        title('save the ROI as "ROI" in the workspace using the brush tool') % 全角度でプロファイルが全部入るようにする。でも大きすぎも良くない。
        ylabel('Intensity');
        xlabel('d spacing [10^{-10} m]');
        hold on
    end
    
    pause % ROIを保存後にenter押して再開
    
    % ROI選択用figure windowが開いていたら消す。
    try
        close(f1)
    catch
    end
    
    % ROIからデータ切り出し
    x_min_ROI = min(ROI(:,1));
    x_max_ROI = max(ROI(:,1));
    
    x_min = find(x_d(:,1)==x_min_ROI);
    x_max = find(x_d(:,1)==x_max_ROI);
    
    % フィッティングしたいデータがうまく切り出せているか確認
    % for i = 1:dataN
    %     plot(x_d([x_min:x_max],i),y([x_min:x_max],i)+10*(i-1))
    %     hold on
    % end
    %
    % pause
    
    %% 切り分けデータを処理後にフィッティング
    f2 = figure('Position',[0 100 1250 750]);
    %     ave_detrend = 5; % detrend'on'と'findpeak'用。 デフォルトは5くらいがいい？最初と最後の5点を使ってデータをdetrend
    %     ave_detrend2 = 2; % 'findpeak'用。デフォルトは2くらいがいい？極小値+/-2点を使ってデータをdetrend
    %     hashi = 30; % ROIの端からhashiポイント以内を探す。
    %     mPP = 0.4; % 極小値の低さの値のthreshold, 大きいほど厳しい1くらいでもいいかも
    
    for j = 1:dataN % dataN
        
        switch average_data  % データの移動平均
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
                % 最初と最後のave_detrendで入力した数の点を使ってデータをdetrend。matlabのdtrend funcitonよりは見た目はマシ。
                % 隣接するピークとうまく分けてROIを設定できた場合はこれでいい。
                ave_detrend = 5; % デフォルトは5くらいがいい？最初と最後の5点を使ってデータをdetrend
                p = polyfit(xx([1:ave_detrend end-ave_detrend:end]),yy([1:ave_detrend end-ave_detrend:end]),1);
                yyy = yy - polyval(p,xx);
                
            case 'off' % デトレンドしない。バックグラウンドが綺麗ならこっちでいいはず。
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
                
                if length(locs) == 2 &&  locs(1) <= hashi && (length(xx) - locs(2)) <= hashi % 両端の最小値が運よく見つかった場合
                    p = polyfit(xx([locs(1)-ave_detrend2:locs(1)+ave_detrend2 locs(2)-ave_detrend2:locs(2)+ave_detrend2]),...
                        yy([locs(1)-ave_detrend2:locs(1)+ave_detrend2 locs(2)-ave_detrend2:locs(2)+ave_detrend2]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 1 % 極小値の有無などの確認用
                elseif length(locs) == 1 && locs <= hashi % leftピークだけが見つかった場合
                    p = polyfit(xx([locs(1)-2:locs(1)+ave_detrend2 end-ave_detrend:end]),...
                        yy([locs(1)-ave_detrend2:locs(1)+ave_detrend2 end-ave_detrend:end]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 2
                elseif length(locs) == 1 && (length(xx) - locs) <= hashi % rightピークだけが見つかった場合
                    p = polyfit(xx([1:ave_detrend locs(1)-ave_detrend2:locs(1)+ave_detrend2]),...
                        yy([1:ave_detrend locs(1)-ave_detrend2:locs(1)+ave_detrend2]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 3
                elseif length(locs) > 2
                    left_part = find(locs(1) <= hashi);
                    [left_h, left_pos] = max(hight(left_part));
                    right_part = find(locs >= (length(xx) - hashi));
                    [right_h, right_pos] = max(hight(right_part));
                    
                    if isempty(left_part(left_pos)) && isnumeric(right_part(right_pos))% left最小値が見つからなかった場合
                        p = polyfit(xx([1:ave_detrend locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),...
                            yy([1:ave_detrend locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),1);
%                         kek = 41
                    elseif isempty(right_part(right_pos)) && isnumeric(left_part(left_pos))% right最小値が見つからなかった場合
                        p = polyfit(xx([locs(left_part(left_pos))-2:locs(left_part(left_pos))+ave_detrend2 end-ave_detrend:end]),...
                            yy([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 end-ave_detrend:end]),1);
%                         kek = 42
                    else % 最小値が両方見つかった場合
                        p = polyfit(xx([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),...
                            yy([locs(left_part(left_pos))-ave_detrend2:locs(left_part(left_pos))+ave_detrend2 locs(right_part(right_pos))-ave_detrend2:locs(right_part(right_pos))+ave_detrend2]),1);
%                         kek = 43
                    end
                    
                    yyy = yy - polyval(p,xx);
                    
                else % 最小値が見つからなかった場合
                    p = polyfit(xx([1:ave_detrend end-ave_detrend:end]),yy([1:ave_detrend end-ave_detrend:end]),1);
                    yyy = yy - polyval(p,xx);
%                     kek = 5
                    
                end
                
            otherwise
                error('Choose a detrend optioin')
        end
        
        
        % psuedo-Voigt関数で各切り分け角度のピークをフィッティング
        [fitresult_peak, gof] = pVoigtFit(xx,yyy);
        fitpeak(j,:) = coeffvalues(fitresult_peak); %A n w xc
        % aが高さ、nがローレンツ関数とガウス関数の比、xがピーク幅、xcがピーク位置
        fitpeak95(j,:) = reshape(confint(fitresult_peak),[1,8]);% 95%信頼区間、a1 min,a1 max, b1 min, b1 max...
        R2_peak(j) = gof.rsquare; %R2値も保存しておく。
        
        % データの近似をプロット
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
    
    %% 応力計算のためのfitting
    rmdata = isoutlier(fitpeak(:,4)); % 外れ値を削除、v1.3で追加
    fitpeak_calc = fitpeak((not(rmdata) & [R2_peak>=threshold_R2]'),4); % 外れ値を削除
    diff_angle_calc = diff_angle((not(rmdata) & [R2_peak>=threshold_R2]')); % 外れ値を削除
    
    [fitresult_stress, gofstress] = debye_stress(diff_angle_calc,fitpeak_calc,tilt_angle); % D0, T = tau/G, U = sigmaU/G
    fitstress = coeffvalues(fitresult_stress); %D0 T U tilt_angle
    fitstress95(1,:) = reshape(confint(fitresult_stress),[1,8]);% 95%信頼区間、a1 min,a1 max, b1 min, b1 max...
    R2_stress = gofstress.rsquare;
    
    %% plot を作成
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
    
    %% もう一度やり直す or 別解析をするかどうか聞く
    c = menu('More calculations?','Yes','No');
    renzoku = c - 1;
    
    if c == 1
        d = menu('Use same XRD data?','Yes','No, load new XRD data');
    else
    end
    
end
%% END
