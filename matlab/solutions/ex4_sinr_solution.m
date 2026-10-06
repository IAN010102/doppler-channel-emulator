function sinr_dB = ex4_sinr_solution(w, h, aJ, sigma2)
% 第 4 關的參考解答：輸出 SINR (dB)。（請先自己寫，卡住再看提示，最後才看這裡！）
% 說明檔頭見 matlab/exercises/ex4_sinr.m。

s  = abs(w' * h)^2;                 % 權重對目標的輸出是一個複數 w'h；功率 = 絕對值的平方。
i  = abs(w' * aJ)^2;                % 干擾的功率，算法相同（aJ 已經含干擾源的振幅）。
nz = sigma2 * norm(w)^2;            % 每根天線有雜訊功率 σ²；加權後總雜訊 = σ² · ‖w‖²。
sinr_dB = 10 * log10( s / (i + nz) );
                                    % 訊號 / (干擾 + 雜訊) 是「功率比」，換 dB 用 10*log10。
                                    % （振幅比才用 20*log10；這裡全部都已經是功率了。）
end
