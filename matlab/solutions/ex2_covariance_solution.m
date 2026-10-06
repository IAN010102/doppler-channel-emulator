function R = ex2_covariance_solution(X)
% 第 2 關的參考解答：樣本協方差 R̂ = (1/L) X Xᴴ。（請先自己寫，卡住再看提示，最後才看這裡！）
% 說明檔頭見 matlab/exercises/ex2_covariance.m。

L = size(X, 2);                     % X 是 N×L：第 2 個維度（欄數）就是快照數 L
R = (X * X') / L;                   % X' 是「共軛轉置」（不只轉置，還把虛部變號）。
                                    % X * X' 把 L 次「x·xᴴ」一次加好；除以 L 才是「平均」。
                                    % 結果 R 是 N×N，且 R 等於 R'（Hermitian）。
end
