function a = ex1_steering_solution(N, d_over_lambda, theta_deg)
% 第 1 關的參考解答：導向向量 a(θ)。（請先自己寫，卡住再看提示，最後才看這裡！）
% 說明檔頭見 matlab/exercises/ex1_steering.m。

n = (0:N-1).';                      % 天線編號 0,1,...,N-1。結尾的 .' 把「列向量」轉成「N×1 欄向量」
                                    % （第 0 根天線當參考點，所以編號從 0 開始，不是從 1 開始）
theta_rad = theta_deg * pi / 180;   % MATLAB 的 sin 吃「弧度」；度 × π/180 = 弧度
a = exp(-1i * 2*pi * n * d_over_lambda * sin(theta_rad));
                                    % 每根天線的相位 = -2π · n · (d/λ) · sin θ。
                                    % 1i 是虛數單位；exp(-1i*φ) 是一個長度為 1、角度為 -φ 的複數，
                                    % 所以 a 的每個元素絕對值都是 1，只有「相位」不同。
end
