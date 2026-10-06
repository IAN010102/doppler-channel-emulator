function w = ex3_mvdr_solution(R, a)
% 第 3 關的參考解答：MVDR 權重 w = R⁻¹a / (aᴴR⁻¹a)。（請先自己寫，卡住再看提示，最後才看這裡！）
% 說明檔頭見 matlab/exercises/ex3_mvdr.m。

x = R \ a;                          % 解線性方程式 R·x = a，得到 x = R⁻¹a。
                                    % 反斜線比 inv(R)*a 更準確也更快；R 要有反矩陣（快照數 L ≥ 天線數 N）。
w = x / (a' * x);                   % a' * x 是一個數字（aᴴR⁻¹a）；除掉它，
                                    % 就保證 w' * a = 1：目標方向的訊號原封不動通過。
end
