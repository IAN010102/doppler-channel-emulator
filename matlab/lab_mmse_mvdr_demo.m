% LAB_MMSE_MVDR_DEMO  講義程式後面可以直接貼上的小段程式：理論協方差的 MMSE 與 MVDR（方法 A、方法 B）。
%
% 怎麼用：
%   (1) 在講義程式裡，已經有這些變數的話，把這個檔案的內容貼在後面執行即可（用的是講義的變數名稱）：
%         array_vec_est1   N×1   使用者 1（目標）的導向向量  exp(-j·2π·n·d/λ·sinθ_user1)，n = 0..N-1
%         array_vec_est2   N×1   使用者 2（干擾）的導向向量（講義若叫別的名字，下面一行改名字即可）
%         Es               純量  符元能量（講義是 1）
%         sigma            純量  雜訊功率 σ（講義的 SNR = 30 dB → sigma = Es · 10^(-30/10)）
%         NR               純量  快照數（講義是 1000）
%         rx               N×NR  接收訊號（每一欄是一次快照）
%   (2) 這些變數不存在時，下面的「若沒有就自己做一份」會用講義的情境（N = 8、d = 0.5λ、θ_user1 = -20°、
%       θ_user2 = 30°、SNR = 30 dB、Es = 1、NR = 1000、干擾功率 = Es）產生，直接按 Run 就能看到結果。
% 角度慣例：講義與網頁模擬器是同一個慣例（以天線法線為基準，導向向量 exp(-j·2π·n·d/λ·sinθ)，
% 正負號相同），所以 θ 不需要轉換。
%
% 注意：下面註解裡「講義的哪個式子」是依照講義的內容描述（MMSE 用理論的 R_r 與 r_rs；MVDR 方法 A 用不含目標的 R_u；
% 方法 B 用含目標的 R_r，實務上用樣本 R̂_r）對應的；式子的編號請你自己對照講義頁面補上。
% 需要：任何 MATLAB 版本；不需要工具箱。

if ~exist('Es', 'var'), Es = 1; end
if ~exist('NR', 'var'), NR = 1000; end
if ~exist('array_vec_est1', 'var') || ~exist('rx', 'var')
    % ---- 若沒有就自己做一份（講義情境）
    N = 8; d = 0.5; th1 = -20; th2 = 30; snr_dB = 30;
    n = (0:N-1).';
    array_vec_est1 = exp(-1i * 2*pi * n * d * sind(th1));      % 使用者 1（目標）：a(θ_user1)
    array_vec_est2 = exp(-1i * 2*pi * n * d * sind(th2));      % 使用者 2（干擾）：a(θ_user2)
    sigma = Es * 10^(-snr_dB / 10);                            % 雜訊功率（SNR = Es/σ = 30 dB）
    rng(1);
    qpsk = @(m) sqrt(Es/2) * (sign(randn(1, m)) + 1i * sign(randn(1, m)));   % 功率 = Es 的 QPSK 符元
    s1 = qpsk(NR); s2 = qpsk(NR);                              % 使用者 1、2 的符元（干擾功率 = Es）
    noise = sqrt(sigma/2) * (randn(N, NR) + 1i * randn(N, NR));
    rx = array_vec_est1 * s1 + array_vec_est2 * s2 + noise;    % 接收訊號 N×NR
end
if ~exist('array_vec_est2', 'var'), error('請先定義干擾的導向向量 array_vec_est2（或把這一行改成講義裡的變數名稱）。'); end
N = numel(array_vec_est1); I = eye(N);

% ---- 理論協方差（講義的 R_r 與 R_u；全是理論值，不是從資料估計）
R_r = Es * (array_vec_est1 * array_vec_est1') + Es * (array_vec_est2 * array_vec_est2') + sigma * I;   % R_r = a1·a1ᴴ·Es + a2·a2ᴴ·Es + σ·I（含目標）
R_u = Es * (array_vec_est2 * array_vec_est2') + sigma * I;                                              % R_u = a2·a2ᴴ·Es + σ·I（不含目標）
r_rs = array_vec_est1 * Es;                                                                             % r_rs = a·Es（接收訊號與期望訊號的互相關）

% ---- MMSE（講義：w = R_r⁻¹ · r_rs）
w_mmse = R_r \ r_rs;

% ---- MVDR 方法 A（用不含目標的 R_u）與方法 B（用含目標的 R_r）
x_A = R_u \ array_vec_est1;  w_A = x_A / (array_vec_est1' * x_A);       % w = R_u⁻¹a / (aᴴR_u⁻¹a)
x_B = R_r \ array_vec_est1;  w_B = x_B / (array_vec_est1' * x_B);       % w = R_r⁻¹a / (aᴴR_r⁻¹a)

% ---- 實務：R_r 用樣本協方差 R̂_r = rx·rxᴴ / NR 代替（方法 B 的實用版）
R_hat = (rx * rx') / NR;
x_Bh = R_hat \ array_vec_est1;  w_B_hat = x_Bh / (array_vec_est1' * x_Bh);

% ---- 檢查：MMSE 與 MVDR 方法 B 方向相同（只差一個純量 Es·(aᴴR_r⁻¹a)）
k = Es * (array_vec_est1' * x_B);                       % 純量 Es · aᴴ R_r⁻¹ a
fprintf('w_MMSE = Es(aᴴR_r⁻¹a)·w_B 的相對誤差：%.2e（應該約 1e-13 以下）\n', norm(w_mmse - k * w_B) / norm(w_mmse));
fprintf('方法 A 與方法 B 的權重差（理論協方差）：%.2e\n', norm(w_A - w_B) / norm(w_B));
fprintf('wᴴa：MMSE = %.4f（沒有 wᴴa = 1 的限制，所以輸出有縮放偏差），MVDR-A = %.4f，MVDR-B = %.4f\n', real(w_mmse' * array_vec_est1), real(w_A' * array_vec_est1), real(w_B' * array_vec_est1));

% ---- 輸出 SINR（用理論的訊號、干擾、雜訊功率；不需要 rx）
sinr = @(w) 10 * log10( Es * abs(w' * array_vec_est1)^2 / ( Es * abs(w' * array_vec_est2)^2 + sigma * norm(w)^2 ) );
fprintf('SINR (dB)：MMSE %.2f，MVDR-A %.2f，MVDR-B %.2f，樣本 R̂_r 的 MVDR-B %.2f\n', sinr(w_mmse), sinr(w_A), sinr(w_B), sinr(w_B_hat));

% ---- 原始 EVM（講義 p.28 的算法：不做增益正規化，直接用 s_hat = wᴴ·rx 與 s 比；需要 s1）
if exist('s1', 'var')
    evm = @(w) sqrt( mean(abs(w' * rx - s1).^2) / mean(abs(s1).^2) );
    fprintf('原始 EVM (%%)：MMSE %.2f，MVDR-A %.2f，MVDR-B %.2f\n', 100 * evm(w_mmse), 100 * evm(w_A), 100 * evm(w_B));
    fprintf('（MMSE 的輸出有縮放偏差 |wᴴa - 1|；MVDR 因為 wᴴa = 1，目標訊號不會被放大或縮小。）\n');
end
