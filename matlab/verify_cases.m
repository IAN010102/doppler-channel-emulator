function nFail = verify_cases()
% VERIFY_CASES  用 MATLAB 自己的程式重算網頁匯出的案例，與網頁的結果逐一比對。
%
% 這支程式「不需要你撰寫」，只要執行並閱讀輸出。
%
% 做了什麼：讀 data/matlab_cases.json（23 個固定案例），對每個案例、每個演算法
%   （FOURIER、MMSE-M、MMSE-P、SMI、DL、BEAMSPACE）：
%     1. 用案例裡的快照 X（或不含訊號的 X_sf）在 MATLAB 重算 R̂；
%     2. 依 JSON 檔頭 definitions / weight_normalisation 寫的公式重算權重 w；
%     3. 與網頁算出的 w 比較；再算 SINR 與 SINR_opt，與網頁比較。
%
% 驗證的是什麼：這驗證的是「線性代數與公式的實作」（R̂、偽逆、六種權重、SINR 公式）網頁端與 MATLAB 一致；
% 它「不」驗證通道產生器（h、a_J、X 本身是網頁產生後匯出的，MATLAB 沒有重新產生）。
% 通道產生器是否符合規格，由 mc_independent.m 另外檢查。
%
% 判定（固定，不放寬）：
%   權重：兩邊各自乘上一個常數使 wᴴ·a_assumed = 1 之後，相對誤差 < 1e-9。
%         （Wiener 類的權重 MMSE-M / MMSE-P 沒有 wᴴa = 1 的規定，所以先正規化才能比。）
%   SINR：絕對誤差 < 1e-9 dB。SINR 與 w 的縮放無關（分子、分母同乘一個常數的平方），
%         所以即使兩邊權重差一個常數倍，SINR 仍應相同。
%   SINR_opt：絕對誤差 < 1e-9 dB。
% 若有案例失敗：不要放寬容許誤差。輸出會列出該案例的條件數 κ(R̂)，並以 κ·2.2e-16 估計
% 「單純因為浮點數捨入，兩邊反矩陣就會差多少」。若誤差大小與這個估計同一個數量級，
% 原因是病態矩陣（例如 L = N 時 κ 可達 1e8），不是公式寫錯。案例都保留在 JSON 裡供檢查。
%
% 需要：MATLAB R2016b 以上（jsondecode），不需要任何工具箱。
%
% 輸出：每個案例每個演算法一行 PASS / FAIL 與最大誤差，最後是總通過數。

here = fileparts(mfilename('fullpath')); root = fileparts(here);
J = jsondecode(fileread(fullfile(root, 'data', 'matlab_cases.json')));
cases = J.cases; cz = @(z) z.re + 1i * z.im;
algos = {'FOURIER', 'MMSE_M', 'MMSE_P', 'SMI', 'DL', 'BEAMSPACE'};   % jsondecode turns 'MMSE-M' into 'MMSE_M'
tolW = 1e-9; tolS = 1e-9;
nTot = 0; nPass = 0; nFail = 0;
fprintf('%-42s %-10s %-6s %-12s %-12s\n', 'case', 'algorithm', 'result', 'max w err', 'SINR err dB');
for k = 1:numel(cases)
    c = cases(k); P = c.params; N = P.N; L = P.L;
    X = cz(c.X); Xsf = cz(c.X_sf); s = cz(c.s); h = cz(c.h); aJ = cz(c.a_J); a = cz(c.a_assumed);
    if strcmp(P.trainMode, 'signalFree'), Xtr = Xsf; else, Xtr = X; end
    Rtr = (Xtr * Xtr') / L;            % R̂ of the data SMI / DL / BEAMSPACE train on (signalFree: without the target)
    Rw  = (X * X') / L;                % R̂ of the Wiener data (always with the target)
    kappa = cond((Rtr + Rtr') / 2);
    for q = 1:numel(algos)
        name = algos{q}; web = c.algorithms.(name); wWeb = cz(web.w);
        switch name
            case 'FOURIER',  w = a / N;
            case 'MMSE_M',   w = pinv_or_solve(Rw, a, P);
            case 'MMSE_P',   r = (X * conj(s)) / L; w = pinv_or_solve(Rw, r, P);
            case 'SMI',      x = pinv_or_solve(Rtr, a, P); w = x / (a' * x);
            case 'DL',       Rd = Rtr + P.gamma_abs * eye(N); x = Rd \ a; w = x / (a' * x);
            case 'BEAMSPACE'
                B = cz(c.beamspace.B); RB = B' * Rtr * B; aB = B' * a; xB = RB \ aB; wB = xB / (aB' * xB); w = B * wB;
        end
        % weights: scale both so that w' * a_assumed = 1, then compare
        wn = w / (a' * w); wWn = wWeb / (a' * wWeb);
        eW = norm(wn - wWn) / norm(wWn);
        sMat = 10 * log10(abs(w' * h)^2 / (P.P_j * abs(w' * aJ)^2 + P.sigma2 * norm(w)^2));
        eS = abs(sMat - web.sinr_dB);
        hr = norm(h)^2 - P.P_j * abs(aJ' * h)^2 / (P.sigma2 + P.P_j * norm(aJ)^2);
        eO = abs(10 * log10(hr / P.sigma2) - web.sinr_opt_dB);
        ok = eW < tolW && eS < tolS && eO < tolS;
        nTot = nTot + 1; if ok, nPass = nPass + 1; else, nFail = nFail + 1; end
        fprintf('%-42s %-10s %-6s %-12.2e %-12.2e\n', c.id, strrep(name, '_', '-'), tern(ok, 'PASS', 'FAIL'), eW, max(eS, eO));
        if ~ok
            fprintf('    -> kappa(R_hat) = %.2e, rounding-only estimate kappa*2.2e-16 = %.2e (if the error is of this size, it is conditioning, not a formula error)\n', kappa, kappa * 2.2e-16);
        end
    end
end
% ---- 講義基準案例 E0：用講義的公式（理論協方差 R_r、R_u）重算 MMSE 與 MVDR 方法 A、B，與網頁的理論版結果比對
ie = find(strcmp({cases.id}, 'E0_lab_baseline'));
if ~isempty(ie)
    c = cases(ie); lab = c.lab; P = c.params; a1 = cz(lab.array_vec_est1); a2 = cz(lab.a2); Es = lab.Es; sigma = lab.sigma; N = numel(a1); I = eye(N);
    h = cz(c.h); aJ = cz(c.a_J);
    Rr = Es * (a1 * a1') + P.P_j * (a2 * a2') + sigma * I;          % R_r = a1·a1ᴴ·Es + a2·a2ᴴ·(干擾功率) + σ·I   （講義的 MMSE 理論協方差；E0 的干擾功率 = Es）
    Ru = P.P_j * (a2 * a2') + sigma * I;                              % R_u = a2·a2ᴴ·(干擾功率) + σ·I               （不含目標）
    wM = Rr \ (a1 * Es);                                              % MMSE：w = R_r⁻¹ r_rs，r_rs = a1·Es
    xB = Rr \ a1; wB = xB / (a1' * xB);                               % MVDR 方法 B（含目標的 R_r）
    xA = Ru \ a1; wA = xA / (a1' * xA);                               % MVDR 方法 A（不含目標的 R_u）
    sinr = @(w) 10 * log10(abs(w' * h)^2 / (P.P_j * abs(w' * aJ)^2 + P.sigma2 * norm(w)^2));
    rel = @(x, y) norm(x - y, 'fro') / norm(y, 'fro');
    tests = {'E0 R_r (理論協方差，含目標)', rel(Rr, cz(lab.R_r_theory)), tolW;
             'E0 R_u (理論協方差，不含目標)', rel(Ru, cz(lab.R_u_theory)), tolW;
             'E0 MMSE 權重 w = R_r^-1 a1 Es', rel(wM, cz(lab.mmse_theory.w)), tolW;
             'E0 MVDR-B 權重', rel(wB, cz(lab.mvdr_B_theory.w)), tolW;
             'E0 MVDR-A 權重', rel(wA, cz(lab.mvdr_A_theory.w)), tolW;
             'E0 SINR MMSE (dB)', abs(sinr(wM) - lab.mmse_theory.sinr_dB), tolS;
             'E0 SINR MVDR-B (dB)', abs(sinr(wB) - lab.mvdr_B_theory.sinr_dB), tolS;
             'E0 SINR MVDR-A (dB)', abs(sinr(wA) - lab.mvdr_A_theory.sinr_dB), tolS;
             'E0 w_MMSE = Es·(a1ᴴR_r⁻¹a1)·w_B', rel(wM, Es * (a1' * xB) * wB), tolW};
    for q = 1:size(tests, 1)
        ok = tests{q, 2} < tests{q, 3}; nTot = nTot + 1; if ok, nPass = nPass + 1; else, nFail = nFail + 1; end
        fprintf('%-42s %-10s %-6s %-12.2e\n', tests{q, 1}, '(lecture)', tern(ok, 'PASS', 'FAIL'), tests{q, 2});
    end
end
fprintf('\n總共 %d 項，通過 %d，未通過 %d。\n', nTot, nPass, nFail);
if nFail == 0, fprintf('全部通過：網頁端的線性代數與公式實作和 MATLAB 一致（通道產生器不在此驗證）。\n'); end
end

function x = pinv_or_solve(R, b, P)
% R \ b, or the pseudo-inverse solution when R_hat is rank deficient and smiSingular = 'pinv' (rank threshold epsRank * largest eigenvalue)
Rh = (R + R') / 2; [V, D] = eig(Rh); d = real(diag(D)); keep = d > P.epsRank * max(d);
if strcmp(P.smiSingular, 'pinv') && sum(keep) < size(R, 1)
    x = V(:, keep) * ((V(:, keep)' * b) ./ d(keep));         % R^+ b = sum over kept eigenpairs of (u' b / lambda) u
else
    x = R \ b;
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
