function w = ex5_diagload(R, a, gamma_rel, sigma2)
% EX5_DIAGLOAD  第 5 關：對角加載 MVDR (diagonal loading)
%
% 先自己寫 → 卡住了再看下面的「提示」→ 最後才看 matlab/solutions/ex5_diagload_solution.m
%
% 【這一關要做什麼】
%   快照太少時，樣本協方差 R̂ 很不準，甚至沒有反矩陣。
%   對角加載的做法很樸素：在 R̂ 的對角線上「加一點點底噪」，讓它一定可以求反矩陣，
%   也讓權重不會被很小的特徵值牽著走。加多少？定成「雜訊功率 σ² 的 γ_rel 倍」。
%
% 【生活比喻】
%   你憑很少的資料猜一個很精細的結論，容易被偶然的數字騙。
%   給每個估計值都加一點「保守的底線」，就比較不會被騙，代價是結論變得稍微粗一點。
%
% 【輸入輸出】
%   輸入  R          N×N 樣本協方差（網頁給的）
%         a          N×1 導向向量
%         gamma_rel  加載比，「倍數」（不是 dB！10 dB 就是 10^(10/10) = 10 倍）
%         sigma2     每根天線的雜訊功率 σ²
%   輸出  w          N×1 權重向量，滿足 w' * a = 1
%
% 【對應數學式】
%       γ = γ_rel · σ²
%       w = (R + γ I)⁻¹ a / ( aᴴ (R + γ I)⁻¹ a )
%   符號：  I   N×N 單位矩陣（MATLAB 寫 eye(N)）
%           γ   真正加到對角線上的數字（單位和 R 的元素相同）
%
% 【提示】（卡住了再往下看，一次只看一條）
%   提示 1：先算 gamma = gamma_rel * sigma2;
%   提示 2：把 R 換成 Rd = R + gamma * eye(size(R, 1));
%   提示 3：接下來和第 3 關完全一樣：x = Rd \ a;  w = x / (a' * x);
%   提示 4：你甚至可以直接呼叫自己寫的 ex3_mvdr(Rd, a)。
%   提示 5：網頁的預設是 +10 dB，也就是 gamma_rel = 10。

w = nan(size(a));                   % 還沒寫時的暫代值（NaN）

%%% TODO 從這裡開始寫（約 3 行）
% gamma = ...;                      % γ = γ_rel · σ²
% Rd = ...;                         % R + γ I
% w = ...;                          % 對 Rd 做 MVDR
%%% TODO 到這裡為止

if any(isnan(w))
    disp('這一關還沒寫（ex5_diagload）：請在 %%% TODO 區填入程式，存檔後再執行 check_my_work。');
end
end
