function selftest()
% SELFTEST  自我測試：把 matlab/solutions/ 的參考解答「暫時當成你的答案」跑一遍 check_my_work，
% 預期 5 / 5 關全部通過。若不是 5/5，代表案例資料或批改程式有問題（不是你的程式有問題）。
%
% 用法：在 MATLAB 指令視窗（目前資料夾 matlab/）輸入   selftest
n = check_my_work('solution');
if n ~= 5
    error('selftest 失敗：參考解答只通過 %d / 5 關。請把上面的輸出貼給助教或開發者。', n);
end
fprintf('selftest OK：參考解答 5 / 5 關通過，批改程式與案例資料自洽。\n');
end
