function w = ex5_diagload_solution(R, a, gamma_rel, sigma2)
% 第 5 關的參考解答：對角加載 MVDR。（請先自己寫，卡住再看提示，最後才看這裡！）
% 說明檔頭見 matlab/exercises/ex5_diagload.m。

gamma = gamma_rel * sigma2;         % 真正加到對角線的數字 γ = (倍數) × (雜訊功率)；
                                    % gamma_rel 是「倍數」，不是 dB（10 dB = 10 倍）。
Rd = R + gamma * eye(size(R, 1));   % 只在對角線上加 γ（eye(N) 是 N×N 單位矩陣），其餘元素不動。
x = Rd \ a;                         % 和第 3 關一樣：先解 Rd·x = a；因為加了 γ，Rd 一定有反矩陣。
w = x / (a' * x);                   % 再除以 aᴴRd⁻¹a，讓 w' * a = 1。
end
