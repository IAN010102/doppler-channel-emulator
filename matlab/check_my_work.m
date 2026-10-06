function nPassed = check_my_work(mode)
% CHECK_MY_WORK  批改你寫的 ex1 ~ ex5：拿網頁模擬器匯出的案例（data/matlab_cases.json）來對答案。
%
% 怎麼用：在 MATLAB 的指令視窗（目前資料夾要是 matlab/）輸入
%       check_my_work
% 它會逐關呼叫 matlab/exercises/ 裡你寫的函式，對照網頁的結果，告訴你每一關
%   通過 / 未通過（附最大誤差與「最可能的原因」）/ 尚未完成（還沒填 %%% TODO，不算失敗）。
% 最後一行是「通過 X / 5 關」。
%
% 比對用的案例：網頁匯出的「訓練資料不含訊號 (signalFree)、靜止 (v = 0)」案例，共 10 個
% （L = 4, 8, 12, 24, 100 各兩個指向偏差 0°、3°）。
%   第 1 關  你的 a(θ)               與案例裡的 a_assumed 比：相對誤差 < 1e-9
%   第 2 關  你的 R̂                  與案例 X_sf 在網頁算出的 R̂ 比：相對誤差 < 1e-9
%   第 3 關  你的 MVDR 權重          用「SINR」與網頁 SMI 比：絕對誤差 < 1e-6 dB，並且要滿足 w' * a = 1（誤差 < 1e-9）
%            （L < N 時 R̂ 沒有反矩陣，網頁用偽逆；這一關只批改 L >= N 的案例）
%   第 4 關  你的 SINR               與網頁的 SINR 比：絕對誤差 < 1e-6 dB
%   第 5 關  你的對角加載權重        用 SINR 與網頁 DL 比：絕對誤差 < 1e-6 dB，並且 w' * a = 1（所有 L）
% 容許誤差是固定的，不會為了讓你過關而放寬。
%
% 為什麼用 SINR 比第 3、5 關？ 權重 w 乘上任何一個常數，輸出 SINR 都不變（分子分母同時乘上同一個
% 常數的平方），所以 SINR 與「權重縮放方式」無關，是公平的比較標準。
%
% check_my_work('solution') 會改用 matlab/solutions/ 的參考解答跑一遍（selftest.m 用它自我測試，應為 5/5）。
%
% 常見錯誤提示（未通過時會自動比對）：
%   第 1 關：相位符號相反、忘了乘 2π、角度忘了換弧度、天線編號從 1 開始、輸出不是 N×1
%   第 2 關：忘了除以 L、共軛放錯（X.' 或 conj）、矩陣大小不對（X'*X）、把虛部丟掉
%   第 3 關：忘了除以 aᴴR⁻¹a、用 R 而不是 R⁻¹、a' 與 a.' 搞混
%   第 4 關：用了 20*log10、w.' 沒共軛、漏了雜訊項或干擾項、忘了平方、忘了 ‖w‖²、忘了換 dB
%   第 5 關：γ 沒乘 σ²、加在整個矩陣而不是只有對角線、忘了正規化

if nargin < 1, mode = 'mine'; end
here = fileparts(mfilename('fullpath')); root = fileparts(here);
if strcmp(mode, 'solution')
    addpath(fullfile(here, 'solutions')); sfx = '_solution';
else
    addpath(fullfile(here, 'exercises')); sfx = '';
end
f1 = str2func(['ex1_steering' sfx]);   f2 = str2func(['ex2_covariance' sfx]); f3 = str2func(['ex3_mvdr' sfx]);
f4 = str2func(['ex4_sinr' sfx]);       f5 = str2func(['ex5_diagload' sfx]);
names = {'導向向量 a(θ)', '樣本協方差 R̂', 'MVDR 權重', 'SINR (dB)', '對角加載 MVDR'};

J = jsondecode(fileread(fullfile(root, 'data', 'matlab_cases.json')));
cases = J.cases; cz = @(z) z.re + 1i * z.im;
tolRel = 1e-9; tolDb = 1e-6;

% status: 1 = pass, 0 = fail, -1 = not done yet, -2 = skipped (L < N); err = largest error; msg = hint of the first failure
status = zeros(5, 1); maxErr = zeros(5, 1); msg = repmat({''}, 5, 1); nDone = zeros(5, 1); failed = false(5, 1); notDone = false(5, 1);
used = 0;
for k = 1:numel(cases)
    c = cases(k); P = c.params;
    if ~strcmp(P.trainMode, 'signalFree') || P.v_kmh ~= 0, continue; end
    used = used + 1;
    N = P.N; L = P.L; Xsf = cz(c.X_sf); Rweb = cz(c.R_hat_train_web); aWeb = cz(c.a_assumed); h = cz(c.h); aJs = sqrt(P.P_j) * cz(c.a_J); s2 = P.sigma2;
    sSMI = c.algorithms.SMI.sinr_dB; sDL = c.algorithms.DL.sinr_dB; wWeb = cz(c.algorithms.SMI.w);
    th = P.theta1_hat_deg + P.delta_theta_deg;

    % ---- 第 1 關
    out = f1(N, P.d_over_lambda, th);
    [st, er, hm] = grade_vec(out, aWeb, tolRel, @() hint1(out, aWeb, N, P.d_over_lambda, th));
    [status, maxErr, msg, notDone, failed] = record(1, st, er, hm, maxErr, msg, notDone, failed);
    % ---- 第 2 關
    out = f2(Xsf);
    [st, er, hm] = grade_mat(out, Rweb, tolRel, @() hint2(out, Rweb, Xsf, L));
    [status, maxErr, msg, notDone, failed] = record(2, st, er, hm, maxErr, msg, notDone, failed);
    % ---- 第 3 關（L < N: 跳過）
    if L >= N
        out = f3(Rweb, aWeb);
        [st, er, hm] = grade_w(out, sSMI, h, aJs, s2, tolDb, aWeb, @() hint3(out, Rweb, aWeb, wWeb));
        [status, maxErr, msg, notDone, failed] = record(3, st, er, hm, maxErr, msg, notDone, failed);
    end
    % ---- 第 4 關
    out = f4(wWeb, h, aJs, s2);
    if ~isnumeric(out) || ~isscalar(out) || isnan(out), st = -1; er = NaN; hm = '';
    else, er = abs(out - sSMI); if er < tolDb, st = 1; hm = ''; else, st = 0; hm = hint4(out, wWeb, h, aJs, s2, sSMI); end
    end
    [status, maxErr, msg, notDone, failed] = record(4, st, er, hm, maxErr, msg, notDone, failed);
    % ---- 第 5 關
    out = f5(Rweb, aWeb, 10^(P.gamma_rel_dB / 10), s2);
    [st, er, hm] = grade_w(out, sDL, h, aJs, s2, tolDb, aWeb, @() hint5(out, Rweb, aWeb, 10^(P.gamma_rel_dB / 10), s2, h, aJs, sDL));
    [status, maxErr, msg, notDone, failed] = record(5, st, er, hm, maxErr, msg, notDone, failed);
end

fprintf('\n===== 批改結果（%d 個案例）=====\n', used);
nPassed = 0;
for q = 1:5
    if notDone(q) && ~failed(q)
        fprintf('第 %d 關 %-14s：尚未完成（還沒填 %%%%%% TODO，不算失敗）\n', q, names{q});
    elseif failed(q)
        fprintf('第 %d 關 %-14s：未通過（最大誤差 %.2e）\n         最可能的原因：%s\n', q, names{q}, maxErr(q), msg{q});
    else
        fprintf('第 %d 關 %-14s：通過（最大誤差 %.2e）\n', q, names{q}, maxErr(q)); nPassed = nPassed + 1;
    end
end
fprintf('通過 %d / 5 關\n', nPassed);
end

% ---------------------------------------------------------------------------------------------- grading helpers
function [status, maxErr, msg, notDone, failed] = record(q, st, er, hm, maxErr, msg, notDone, failed)
status = 0;
if st == -1, notDone(q) = true;
elseif st == 0, failed(q) = true; if isempty(msg{q}), msg{q} = hm; end
end
if ~isnan(er), maxErr(q) = max(maxErr(q), er); end
end

function [st, er, hm] = grade_vec(out, ref, tol, hintFn)
hm = '';
if ~isnumeric(out) || any(isnan(out(:))), st = -1; er = NaN; return; end
if ~isequal(size(out), size(ref)), st = 0; er = Inf; hm = hintFn(); if isempty(hm), hm = '輸出的大小不對：應該是 N×1 的欄向量。'; end; return; end
er = norm(out - ref) / norm(ref);
if er < tol, st = 1; else, st = 0; hm = hintFn(); if isempty(hm), hm = '和網頁的結果不一致：請對照【對應數學式】逐項檢查。'; end; end
end

function [st, er, hm] = grade_mat(out, ref, tol, hintFn)
hm = '';
if ~isnumeric(out) || any(isnan(out(:))), st = -1; er = NaN; return; end
if ~isequal(size(out), size(ref)), st = 0; er = Inf; hm = hintFn(); if isempty(hm), hm = '矩陣大小不對：R̂ 應該是 N×N。'; end; return; end
er = norm(out - ref, 'fro') / norm(ref, 'fro');
if er < tol, st = 1; else, st = 0; hm = hintFn(); if isempty(hm), hm = '和網頁的結果不一致：請對照【對應數學式】逐項檢查。'; end; end
end

function [st, er, hm] = grade_w(out, sinrRef, h, aJs, s2, tol, a, hintFn)
% 通過條件：(1) 權重的 SINR 與網頁相同（SINR 與縮放無關）且 (2) 權重滿足題目寫的輸出規定 w' * a = 1（誤差 < 1e-9）。
hm = '';
if ~isnumeric(out) || any(isnan(out(:))), st = -1; er = NaN; return; end
if numel(out) ~= numel(h), st = 0; er = Inf; hm = '權重 w 的大小不對：應該是 N×1。'; return; end
er = abs(local_sinr(out(:), h, aJs, s2) - sinrRef); con = abs(a' * out(:) - 1);
if er < tol && con < 1e-9, st = 1; else, st = 0; er = max(er, con); hm = hintFn(); if isempty(hm), hm = '權重算出來的 SINR 和網頁不同：請檢查公式與正規化（w'' * a 是否等於 1）。'; end; end
end

function s = local_sinr(w, h, aJs, s2)       % 批改用的 SINR（不使用你寫的 ex4，讓各關獨立批改）
s = 10 * log10(abs(w' * h)^2 / (abs(w' * aJs)^2 + s2 * norm(w)^2));
end

% ---------------------------------------------------------------------------------------------- hints for common mistakes
function m = hint1(out, ref, N, d, th)
m = ''; n = (0:N-1).'; tr = th * pi / 180;
if size(out, 1) ~= N || size(out, 2) ~= 1, m = '輸出不是 N×1 的欄向量（用 (0:N-1).'' 或 (0:N-1)(:) 做欄向量）。'; return; end
if norm(out - conj(ref)) / norm(ref) < 1e-6, m = '導向向量的相位符號相反（式子裡 exp 的指數前面要有負號）。'; return; end
if norm(out - exp(-1i * n * d * sin(tr))) / norm(ref) < 1e-6, m = '忘了乘 2π（相位是 2π · n · (d/λ) · sin θ）。'; return; end
if norm(out - exp(-1i * 2*pi * n * d * sin(th))) / norm(ref) < 1e-6, m = '角度忘了從「度」換成「弧度」（theta_deg * pi / 180）。'; return; end
if abs(abs(out' * ref) / (norm(out) * norm(ref)) - 1) < 1e-6, m = '天線編號可能從 1 開始了：第一根天線是第 0 根（n = 0:N-1），它的相位要是 0。'; return; end
end

function m = hint2(out, ref, X, L)
m = '';
if size(out, 1) ~= size(ref, 1) || size(out, 2) ~= size(ref, 2), m = '矩陣大小不對：X'' * X 會得到 L×L；要用 X * X''（N×N）。'; return; end
if norm(out - ref * L, 'fro') / norm(ref * L, 'fro') < 1e-6, m = '忘了除以 L（R̂ 是 L 次的「平均」，不是總和）。'; return; end
if norm(out - conj(ref), 'fro') / norm(ref, 'fro') < 1e-6 || norm(out - ref.', 'fro') / norm(ref, 'fro') < 1e-6
    m = '共軛放錯邊了（應該是 X * X''，單引號是共軛轉置；X.'' 或 conj(X) * X.'' 是不同的東西）。'; return;
end
if norm(out - real(ref), 'fro') / norm(ref, 'fro') < 1e-6, m = '虛部被丟掉了（不要用 real 或 abs，保留複數）。'; return; end
end

function m = hint3(out, R, a, wWeb)
m = ''; out = out(:);
sc = a' * out;
if abs(sc - 1) > 1e-6 && norm(out / sc - wWeb) / norm(wWeb) < 1e-5, m = '忘了除以 aᴴR⁻¹a：w 要乘上一個常數，使 w'' * a = 1。'; return; end
wA = R * a; wA = wA / (a' * wA);
if norm(out - wA) / norm(wA) < 1e-6, m = '用了 R 而不是 R⁻¹：要解 R \ a，不是 R * a。'; return; end
wB = R \ conj(a); wB = wB / (a.' * wB);
if norm(out - wB) / max(norm(wB), eps) < 1e-6, m = 'a'' 與 a.'' 搞混了：a'' 是共軛轉置，這裡要用 a''。'; return; end
end

function m = hint4(out, w, h, aJs, s2, ref)
m = ''; s = abs(w' * h)^2; i = abs(w' * aJs)^2; nz = s2 * norm(w)^2;
cand = {20*log10(s / (i + nz)), '用了 20*log10；功率比要用 10*log10（20 是用在振幅比）。';
        10*log10(abs(w.' * h)^2 / (abs(w.' * aJs)^2 + nz)), 'w.'' 沒有共軛；權重要用共軛轉置 w''。';
        10*log10(s / i), '漏了雜訊項 σ² · ‖w‖²。';
        10*log10(s / nz), '漏了干擾項 |wᴴ a_J|²。';
        10*log10(abs(w' * h) / (abs(w' * aJs) + sqrt(nz))), '忘了平方：功率 = 絕對值的平方。';
        10*log10(s / (i + s2)), '雜訊項漏了 ‖w‖²（要乘 norm(w)^2）。';
        s / (i + nz), '忘了換成 dB：要 10*log10( ... )。'};
for q = 1:size(cand, 1)
    if abs(out - cand{q, 1}) < 1e-6 && abs(cand{q, 1} - ref) > 1e-6, m = cand{q, 2}; return; end
end
end

function m = hint5(out, R, a, gr, s2, h, aJs, ref)
m = ''; out = out(:); N = size(R, 1); mine = local_sinr(out, h, aJs, s2);
alt = {R + gr * eye(N),             'γ 沒乘 σ²：γ = γ_rel · σ²（對角線上加的是 γ_rel * sigma2）。';
       R + gr * s2 * ones(N),        '加在整個矩陣上了：只能加在「對角線」上（用 eye(N)）。'};
for q = 1:size(alt, 1)
    x = alt{q, 1} \ a; wq = x / (a' * x);
    if abs(local_sinr(wq, h, aJs, s2) - mine) < 1e-6, m = alt{q, 2}; return; end
end
Rd = R + gr * s2 * eye(N); x = Rd \ a;
if abs(a' * out - 1) > 1e-6 && norm(out / (a' * out) - x / (a' * x)) / norm(x / (a' * x)) < 1e-5, m = '忘了除以 aᴴ(R+γI)⁻¹a：w 要縮放到 w'' * a = 1。'; return; end
end
