function mc_independent()
% MC_INDEPENDENT  在 MATLAB 裡「依規格」獨立產生通道與訊號，重跑 SMI、DL 的平均 SINR 對快照數 L 的曲線，
%                 再和網頁模擬器的結果（data/curve_sinr_vs_L.csv）逐點比較。
%
% 這支程式「不需要你撰寫」，只要執行並閱讀輸出。
%
% 這個比對檢驗的是「實作是否符合規格」，不檢驗規格本身。
% （規格 = PARAMS.md 寫的統一模型。MATLAB 與網頁各自照同一份規格寫，結果若吻合，代表兩邊對規格的
%   理解和實作一致；但規格本身是否符合真實物理，這個比對無法回答。）
%
% 規格（對應 PARAMS.md §3）：
%   v = 0（靜止），N = 8 根天線，間距 d = 0.5λ，K = 20 dB，SNR = 20 dB（每根天線），SIR = -10 dB，
%   目標在 θ1 = 0°，干擾在 θ2 = 40°，天線無相位誤差，目標符元為 16-QAM（平均功率 1），
%   干擾符元為複高斯 CN(0, P_j)（jamWave = gaussian），雜訊為複高斯 CN(0, σ²)。
%   一次試驗只抽一次通道：LoS 路徑增益 sqrt(K/(K+1))·exp(jφ0)，φ0 在 [0, 2π) 均勻；
%   M = 32 條擴散路徑，角度 θ1 + δ_m（δ_m ~ N(0, σθ²)，σθ = 10°），增益 CN(0, 1/((K+1)M))；
%   h = Σ_i β_i a(θ_i)。SINR = |wᴴh|² / (P_j|wᴴa_J|² + σ²‖w‖²)（用真實通道 h，不是樣本）。
%   SMI-MVDR：w = R̂⁻¹a/(aᴴR̂⁻¹a)；訓練資料 signalFree = 只有干擾 + 雜訊，withSignal = 還含目標。
%   DL：R̂ + γI，γ = γ_rel·σ²，γ_rel = +10 dB。L < N 時 R̂ 秩虧：SMI 用偽逆（特徵值門檻 1e-10·最大特徵值）。
%
% 比對方法：對 (訓練模式 × 演算法 × L) 的每個點算 z = (MATLAB 平均 − 網頁平均) / sqrt(SE² + SE²)。
%   比較個數 m，用 Bonferroni 校正的臨界值（α = 0.01，雙尾）：z_crit = sqrt(2)·erfcinv(α/m)；
%   最大 |z| < z_crit 才算 PASS。（同時比較 m 個統計量，最大 |z| 的期望值會隨 m 增加，
%   所以不能用固定的門檻 3。）輸出 m、z_crit、最大 |z| 與 PASS/FAIL，並畫兩條曲線疊圖到
%   figs/mc_independent_check.pdf。
%
% 需要：MATLAB R2016b 以上，不需要任何工具箱。2000 次 × 10 個 L，通常幾十秒到數分鐘。

here = fileparts(mfilename('fullpath')); root = fileparts(here);
rng(20260601);
N = 8; Kdb = 20; K = 10^(Kdb / 10); M = 32; sigTh = 10 * pi / 180; sigma2 = 10^(-20 / 10); Pj = 10^(10 / 10);
th1 = 0; th2 = 40 * pi / 180; gammaRel = 10^(10 / 10); epsRank = 1e-10; nTrials = 2000; alpha = 0.01;
Ls = [2 4 6 8 10 12 16 24 32 48]; Lmax = max(Ls);
modes = {'signalFree', 'withSignal'}; algos = {'SMI', 'DL'};
steer = @(th) exp(-1i * 2*pi * (0:N-1).' * 0.5 * sin(th));
aJ = steer(th2);
lev = [-3 -1 1 3] / sqrt(10);

sums = zeros(numel(modes), numel(algos), numel(Ls)); sq = sums;
for t = 1:nTrials
    % one channel realisation per trial
    phi0 = 2*pi * rand; beta0 = sqrt(K / (K + 1)) * exp(1i * phi0);
    dl = sigTh * randn(M, 1); g = sqrt(1 / ((K + 1) * M) / 2) * (randn(M, 1) + 1i * randn(M, 1));
    h = beta0 * steer(th1);
    for m = 1:M, h = h + g(m) * steer(th1 + dl(m)); end
    s = lev(randi(4, 1, Lmax)) + 1i * lev(randi(4, 1, Lmax));                          % 16-QAM, unit power (1 x Lmax)
    jw = sqrt(Pj / 2) * (randn(1, Lmax) + 1i * randn(1, Lmax));                         % Gaussian jammer symbols
    nz = sqrt(sigma2 / 2) * (randn(N, Lmax) + 1i * randn(N, Lmax));
    Xsf = aJ * jw + nz; Xfull = Xsf + h * s;
    for im = 1:numel(modes)
        if im == 1, Xall = Xsf; else, Xall = Xfull; end
        for il = 1:numel(Ls)
            L = Ls(il); X = Xall(:, 1:L); R = (X * X') / L; a = steer(th1);
            [V, D] = eig((R + R') / 2); d = real(diag(D)); keep = d > epsRank * max(d);
            if sum(keep) < N, x = V(:, keep) * ((V(:, keep)' * a) ./ d(keep)); else, x = R \ a; end
            wS = x / (a' * x);
            xd = (R + gammaRel * sigma2 * eye(N)) \ a; wD = xd / (a' * xd);
            for ia = 1:2
                if ia == 1, w = wS; else, w = wD; end
                v = 10 * log10(abs(w' * h)^2 / (Pj * abs(w' * aJ)^2 + sigma2 * norm(w)^2));
                sums(im, ia, il) = sums(im, ia, il) + v; sq(im, ia, il) = sq(im, ia, il) + v^2;
            end
        end
    end
end
mu = sums / nTrials; se = sqrt(max(sq / nTrials - mu.^2, 0) / (nTrials - 1));

% ---- compare with the web simulator's curve
T = readtable(fullfile(root, 'data', 'curve_sinr_vs_L.csv'), 'CommentStyle', '#');
zs = []; rows = {};
for im = 1:numel(modes), for ia = 1:numel(algos), for il = 1:numel(Ls)
    sel = strcmp(T.train_mode, modes{im}) & strcmp(T.algorithm, algos{ia}) & T.L == Ls(il);
    mw = T.sinr_mean_dB(sel); sw = T.sinr_se_dB(sel);
    z = (mu(im, ia, il) - mw) / sqrt(se(im, ia, il)^2 + sw^2);
    zs(end + 1) = z; %#ok<AGROW>
    rows(end + 1, :) = {modes{im}, algos{ia}, Ls(il), mu(im, ia, il), se(im, ia, il), mw, sw, z}; %#ok<AGROW>
end, end, end
mCmp = numel(zs); zCrit = sqrt(2) * erfcinv(alpha / mCmp); zMax = max(abs(zs));
fprintf('%-11s %-5s %4s %10s %8s %10s %8s %7s\n', 'train', 'algo', 'L', 'MATLAB dB', 'SE', 'web dB', 'SE', 'z');
for q = 1:size(rows, 1), fprintf('%-11s %-5s %4d %10.3f %8.3f %10.3f %8.3f %7.2f\n', rows{q, :}); end
pass = zMax < zCrit;
fprintf('\n比較個數 m = %d，Bonferroni 臨界值 z_crit = %.3f（alpha = %.2f，雙尾），最大 |z| = %.3f  ->  %s\n', mCmp, zCrit, alpha, zMax, tern(pass, 'PASS', 'FAIL'));
if ~pass, fprintf('FAIL：不要調整門檻。請檢查哪一列的 z 最大，並與 PARAMS.md §3 的規格逐項對照（符元、干擾波形、K、擴散路徑數）。\n'); end

% ---- overlay figure
if ~exist(fullfile(root, 'figs'), 'dir'), mkdir(fullfile(root, 'figs')); end
fig = figure('Units', 'inches', 'Position', [1 1 3.5 4.6], 'Color', 'w');
sty = {'-', 'k', 'o'; '--', [0.8 0.4 0], 's'};
for im = 1:2
    ax = subplot(2, 1, im); hold(ax, 'on'); box(ax, 'on');
    for ia = 1:2
        sel = strcmp(T.train_mode, modes{im}) & strcmp(T.algorithm, algos{ia});
        plot(ax, T.L(sel), T.sinr_mean_dB(sel), sty{ia, 1}, 'Color', sty{ia, 2}, 'LineWidth', 0.8);
        errorbar(ax, Ls, squeeze(mu(im, ia, :)), squeeze(se(im, ia, :)), sty{ia, 3}, 'Color', sty{ia, 2}, 'MarkerSize', 3, 'LineStyle', 'none');
    end
    xlabel(ax, 'Snapshots L'); ylabel(ax, 'Mean SINR (dB)'); title(ax, modes{im}, 'FontWeight', 'normal');
    set(ax, 'FontName', 'Helvetica', 'FontSize', 8);
    if im == 1, legend(ax, {'SMI web', 'SMI MATLAB', 'DL web', 'DL MATLAB'}, 'Location', 'southeast', 'FontSize', 6); end
end
set(fig, 'PaperUnits', 'inches', 'PaperSize', [3.5 4.6], 'PaperPosition', [0 0 3.5 4.6]);
print(fig, fullfile(root, 'figs', 'mc_independent_check.pdf'), '-dpdf', '-painters');
fprintf('圖已存到 figs/mc_independent_check.pdf（線 = 網頁，符號 + 誤差棒 = MATLAB，誤差棒為標準誤）。\n');
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
