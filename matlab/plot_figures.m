function plot_figures()
% PLOT_FIGURES  讀 data/curve_*.csv，輸出報告用的圖：figs/*.pdf（向量格式）與 figs/*.png（300 dpi）。
%
% 這支程式「不需要你撰寫」，只要執行並看圖。
%
% 圖的格式：單欄寬度 3.5 吋；Helvetica 8 pt；每個演算法在所有圖裡用同一種顏色 + 線型 + 標記
% （黑白列印也能分辨：線型與標記不同）；座標軸標籤含單位；圖例放在圖外（不遮擋曲線）；
% 陰影帶 = 平均值 ± 1 個標準誤（SE，standard error of the mean；試驗數在各 CSV 的 n_trials）。
%   Fig 1  EVM 對車速（含 ICI 極限），兩組 (fc, Δf) 各一個子圖
%   Fig 2  SINR 對快照數 L（含 SINR_opt 虛線），訓練資料 signalFree / withSignal 各一個子圖
%   Fig 3  SINR 對指向偏差 δθ（L = 100），signalFree / withSignal 各一個子圖
%   Fig 3b 同 Fig 3，L = 12
%   Fig 4  DOA 空間譜（Capon 與 MUSIC），標出真實角度
% 需要：MATLAB R2016b 以上（readtable、print），不需要任何工具箱。
% 先在 repo 根目錄執行  node tests/make_curves.js  產生 data/curve_*.csv（repo 內已附一份）。

here = fileparts(mfilename('fullpath')); root = fileparts(here);
outDir = fullfile(root, 'figs'); if ~exist(outDir, 'dir'), mkdir(outDir); end
rd = @(f) readtable(fullfile(root, 'data', f), 'CommentStyle', '#');

% one style per algorithm, the same in every figure: {color, line style, marker}
ST = struct( ...
    'FOURIER',   {{[0 0 0],        '-',  'o'}}, ...
    'MMSE_M',    {{[0.15 0.25 0.8], '--', 's'}}, ...
    'MMSE_P',    {{[0.8 0.15 0.15], '-.', '^'}}, ...
    'SMI',       {{[0 0.5 0.1],     ':',  'd'}}, ...
    'DL',        {{[0.85 0.5 0],    '-',  'v'}}, ...
    'BEAMSPACE', {{[0.5 0 0.55],    '--', 'x'}}, ...
    'SINR_opt',  {{[0.35 0.35 0.35], ':', 'none'}}, ...
    'ICI',       {{[0.35 0.35 0.35], ':', 'none'}});
key = @(name) strrep(name, '-', '_');
setup = @() set(groot, 'defaultAxesFontName', 'Helvetica', 'defaultAxesFontSize', 8, 'defaultTextFontName', 'Helvetica', 'defaultTextFontSize', 8, 'defaultLineLineWidth', 0.8);
setup();

% ---- Fig 1: EVM vs speed ----------------------------------------------------------------------------------
T = rd('curve_evm_vs_velocity.csv'); algs = {'FOURIER', 'MMSE-M', 'MMSE-P', 'SMI', 'DL', 'BEAMSPACE'};
fig = newfig(3.5, 4.6);
pairs = [5 15; 28 120];
for p = 1:2
    ax = subplot(2, 1, p); hold(ax, 'on'); box(ax, 'on'); h = [];
    for a = 1:numel(algs)
        sel = T.fc_GHz == pairs(p, 1) & T.delta_f_kHz == pairs(p, 2) & strcmp(T.algorithm, algs{a});
        h(end + 1) = band(ax, T.velocity_kmh(sel), T.evm_mean_pct(sel), T.evm_se_pct(sel), ST.(key(algs{a}))); %#ok<AGROW>
    end
    sel = T.fc_GHz == pairs(p, 1) & T.delta_f_kHz == pairs(p, 2) & strcmp(T.algorithm, 'SMI');
    h(end + 1) = plot(ax, T.velocity_kmh(sel), T.ici_floor_evm_pct(sel), ST.ICI{2}, 'Color', ST.ICI{1}, 'LineWidth', 1.2); %#ok<AGROW>
    xlabel(ax, 'Speed (km/h)'); ylabel(ax, 'EVM (%)'); title(ax, sprintf('f_c = %g GHz, \\Deltaf = %g kHz', pairs(p, 1), pairs(p, 2)), 'FontWeight', 'normal');
    if p == 1, legendBelow(ax, h, [algs, {'ICI limit'}]); end
end
saveFig(fig, outDir, 'fig1_evm_vs_speed');

% ---- Fig 2: SINR vs L ---------------------------------------------------------------------------------------
T = rd('curve_sinr_vs_L.csv'); algs = {'SMI', 'DL', 'BEAMSPACE', 'MMSE-P', 'SINR_opt'};
fig = newfig(3.5, 4.6); modes = {'signalFree', 'withSignal'};
for m = 1:2
    ax = subplot(2, 1, m); hold(ax, 'on'); box(ax, 'on'); h = [];
    for a = 1:numel(algs)
        sel = strcmp(T.train_mode, modes{m}) & strcmp(T.algorithm, algs{a});
        h(end + 1) = band(ax, T.L(sel), T.sinr_mean_dB(sel), T.sinr_se_dB(sel), ST.(key(algs{a}))); %#ok<AGROW>
    end
    xline_(ax, 8); xlabel(ax, 'Snapshots L'); ylabel(ax, 'SINR (dB)'); title(ax, ['training data: ' modes{m}], 'FontWeight', 'normal');
    if m == 1, legendBelow(ax, h, strrep(algs, '_', ' ')); end
end
saveFig(fig, outDir, 'fig2_sinr_vs_L');

% ---- Fig 3 / 3b: SINR vs delta_theta -----------------------------------------------------------------------
T = rd('curve_sinr_vs_delta_theta.csv'); algs = {'SMI', 'DL', 'BEAMSPACE', 'MMSE-M', 'MMSE-P', 'SINR_opt'};
for L = [100 12]
    fig = newfig(3.5, 4.6);
    for m = 1:2
        ax = subplot(2, 1, m); hold(ax, 'on'); box(ax, 'on'); h = [];
        for a = 1:numel(algs)
            sel = strcmp(T.train_mode, modes{m}) & T.L == L & strcmp(T.algorithm, algs{a});
            h(end + 1) = band(ax, T.delta_theta_deg(sel), T.sinr_mean_dB(sel), T.sinr_se_dB(sel), ST.(key(algs{a}))); %#ok<AGROW>
        end
        xlabel(ax, 'Pointing error \delta\theta (deg)'); ylabel(ax, 'SINR (dB)'); title(ax, sprintf('%s, L = %d', modes{m}, L), 'FontWeight', 'normal');
        if m == 1, legendBelow(ax, h, strrep(algs, '_', ' ')); end
    end
    saveFig(fig, outDir, sprintf('fig3_sinr_vs_delta_L%d', L));
end

% ---- Fig 4: DOA spectrum -------------------------------------------------------------------------------------
T = rd('curve_doa_spectrum.csv'); fig = newfig(3.5, 2.6); ax = axes(fig); hold(ax, 'on'); box(ax, 'on');
h1 = band(ax, T.angle_deg, T.capon_mean_dB, T.capon_se_dB, {[0 0.45 0.75], '-', 'none'});
h2 = band(ax, T.angle_deg, T.music_mean_dB, T.music_se_dB, {[0.8 0.15 0.15], '--', 'none'});
ylim(ax, [-50 3]); xlim(ax, [-90 90]);
h3 = plot(ax, [-10 -10], [-50 3], 'k-.', 'LineWidth', 0.8); plot(ax, [40 40], [-50 3], 'k-.', 'LineWidth', 0.8);
text(ax, -9, -4, '\theta_1', 'FontSize', 8); text(ax, 41, -4, '\theta_2', 'FontSize', 8);
xlabel(ax, 'Angle (deg)'); ylabel(ax, 'Spectrum (dB, rel. to peak)');
legendBelow(ax, [h1 h2 h3], {'Capon', 'MUSIC (2 sources)', 'true angles'});
saveFig(fig, outDir, 'fig4_doa_spectrum');
fprintf('圖已存到 %s（pdf = 向量格式，png = 300 dpi）。陰影 = 平均 ± 1 SE。\n', outDir);
end

% ----------------------------------------------------------------------------------------------------------------
function fig = newfig(w, h)
fig = figure('Units', 'inches', 'Position', [1 1 w h], 'Color', 'w', 'PaperUnits', 'inches', 'PaperSize', [w h], 'PaperPosition', [0 0 w h]);
end

function hline = band(ax, x, y, se, st)
% mean curve with a +-1 SE band; st = {color, line style, marker}
x = x(:); y = y(:); se = se(:);
patch(ax, [x; flipud(x)], [y + se; flipud(y - se)], st{1}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
k = unique(round(linspace(1, numel(x), min(numel(x), 7))));
hline = plot(ax, x, y, 'LineStyle', st{2}, 'Color', st{1}, 'LineWidth', 0.9);
if ~strcmp(st{3}, 'none'), set(hline, 'Marker', st{3}, 'MarkerIndices', k, 'MarkerSize', 3.5); end   % markers only at a few points (MarkerIndices: R2016b+)
end

function xline_(ax, x0)
yl = get(ax, 'YLim'); plot(ax, [x0 x0], yl, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); set(ax, 'YLim', yl);
end

function legendBelow(ax, h, labels)
% legend outside, below the axes (does not cover the curves)
try
    lg = legend(ax, h, labels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 6.5, 'Box', 'off');
    try, lg.NumColumns = 3; catch, end
catch
    legend(ax, h, labels, 'Location', 'best', 'FontSize', 6.5);
end
end

function saveFig(fig, outDir, name)
print(fig, fullfile(outDir, [name '.pdf']), '-dpdf', '-painters');
print(fig, fullfile(outDir, [name '.png']), '-dpng', '-r300');
close(fig);
end
