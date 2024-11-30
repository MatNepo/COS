% исходные данные
Fm = 10;        % частота модулируемого сигнала
Fn = 360;       % Несущая частота
m = 1;          % коэффициент модуляции 0 < m <= 1
Phin = 0;       % фаза несущей
Phim = 0;       % фаза модулируемого сигнала
I = 256;        % число уровней квантования
Kdiscr = 8;     % отношение частоты дискретизации к частоте несущей
Fd = Fn*Kdiscr; % частота дискретизации
Td = 1/Fd;      % период дискретизации
tend = 2;       % время окончания модуляции сигнала
 
% Генерация сигнала
t = 0:Td:tend;
N = length(t);  % общее количество отсчетов
mod_sig = 0.5*(1+m.*sin(2*pi*Fm.*t+Phim)).*sin(2*pi*Fn.*t+Phin);
 
% Построение графика модулированного сигнала
figure;
plot(t, mod_sig);
title('Модулированный сигнал');
xlabel('t, сек');
ylabel('U, B');
axis([0 tend -1 1]);
 
% Создание шума
noise10 = 0.2*rand(1,N)-0.1;
noise35 = 0.7*rand(1,N)-0.35;
noise65 = 1.3*rand(1,N)-0.65;
 
% Добавление шума к сигналу
mod_sig_10 = mod_sig + noise10;
mod_sig_35 = mod_sig + noise35;
mod_sig_65 = mod_sig + noise65;
 
% Построение графиков зашумленных сигналов
figure;         % Шум 10%
plot(t, mod_sig_10);
title('Зашумленный сигнал 10%');
xlabel('t, сек');
ylabel('Амплитуда');
axis([0 tend -1.1 1.1]);
figure;         % Шум 35%
plot(t, mod_sig_35);
title('Зашумленный сигнал 35%');
xlabel('t, сек');
ylabel('Амплитуда');
axis([0 tend -1.4 1.4]);
figure;         % Шум 65%
plot(t, mod_sig_65);
title('Зашумленный сигнал 65%');
xlabel('t, сек');
ylabel('Амплитуда');
axis([0 tend -1.7 1.7]);
 
% АЦП
% Дискретизация сигнала
y_10 = floor((mod_sig_10+m)*I/(2*m));
y_35 = floor((mod_sig_35+m)*I/(2*m));
y_65 = floor((mod_sig_65+m)*I/(2*m));
 
% Построение графиков дискретизации
figure;         % Сигнал после АЦП 10%
plot(t, y_10);
title('Сигнал после АЦП 10%');
xlabel('t, сек');
ylabel('Квантованные значения сигнала');
axis([0 tend 0 I]);
figure;         % Сигнал после АЦП 35%
plot(t, y_35);
title('Сигнал после АЦП 35%');
xlabel('t, сек');
ylabel('Квантованные значения сигнал');
axis([0 tend 0 I]);
figure;         % Сигнал после АЦП 65%
plot(t, y_65);
title('Сигнал после АЦП 65%');
xlabel('t, сек');
ylabel('Квантованные значения сигнал');
axis([0 tend 0 I]);
 
% первый перенос частоты - получим исходный сигнал, убрав несущую
% домножениена Sin и Cos
for i = 1:N
    sin_out_10(i) = y_10(i)*sin((i-1) * 2*pi*Fn/Fd);
    cos_out_10(i) = y_10(i)*cos((i-1) * 2*pi*Fn/Fd);
    
    sin_out_35(i) = y_35(i)*sin((i-1) * 2*pi*Fn/Fd);
    cos_out_35(i) = y_35(i)*cos((i-1) * 2*pi*Fn/Fd);
    
    sin_out_65(i) = y_65(i)*sin((i-1) * 2*pi*Fn/Fd);
    cos_out_65(i) = y_65(i)*cos((i-1) * 2*pi*Fn/Fd);
end
% пропустим результаты домножения на Sin и Cos через БФ 6 порядка
[b_6, a_6] = butter(6, Fm*Kdiscr/(2*Fd));
sin_out_butt_10 = filter(b_6, a_6, sin_out_10);
cos_out_butt_10 = filter(b_6, a_6, cos_out_10);
 
sin_out_butt_35 = filter(b_6, a_6, sin_out_35);
cos_out_butt_35 = filter(b_6, a_6, cos_out_35);
 
sin_out_butt_65 = filter(b_6, a_6, sin_out_65);
cos_out_butt_65 = filter(b_6, a_6, cos_out_65);
% операция детектирования
detection6_10 = sqrt(sin_out_butt_10.^2 + cos_out_butt_10.^2);
detection6_35 = sqrt(sin_out_butt_35.^2 + cos_out_butt_35.^2);
detection6_65 = sqrt(sin_out_butt_65.^2 + cos_out_butt_65.^2);
 
% Построение графиков первой фильтрации
figure;
 
plot(t, detection6_10, t, detection6_35, t, detection6_65);
title('Результат первого переноса частоты');
xlabel('t, сек');
ylabel('Амплитуда');
legend('Шум 10%', 'Шум 35%', 'Шум 65%');
 
% второй перенос частоты
% домножение на Sin и Cos
for i = 1:N
    sin_out2_10(i) = detection6_10(i)*sin((i-1) * 2*pi*Fm/Fd);
    cos_out2_10(i) = detection6_10(i)*cos((i-1) * 2*pi*Fm/Fd);
    
    sin_out2_35(i) = detection6_35(i)*sin((i-1) * 2*pi*Fm/Fd);
    cos_out2_35(i) = detection6_35(i)*cos((i-1) * 2*pi*Fm/Fd);
    
    sin_out2_65(i) = detection6_65(i)*sin((i-1) * 2*pi*Fm/Fd);
    cos_out2_65(i) = detection6_65(i)*cos((i-1) * 2*pi*Fm/Fd);
end
% пропустим результаты домножения на Sin и Cos через БФ 6 порядка
[b_6, a_6] = butter(6, Fm*Kdiscr/(2*8*Fd));
sin_out2_butt_10 = filter(b_6, a_6, sin_out2_10);
cos_out2_butt_10 = filter(b_6, a_6, cos_out2_10);
 
sin_out2_butt_35 = filter(b_6, a_6, sin_out2_35);
cos_out2_butt_35 = filter(b_6, a_6, cos_out2_35);
 
sin_out2_butt_65 = filter(b_6, a_6, sin_out2_65);
cos_out2_butt_65 = filter(b_6, a_6, cos_out2_65);
% операция детектирования
detection_26_10 = sqrt(sin_out2_butt_10.^2 + cos_out2_butt_10.^2);
detection_26_35 = sqrt(sin_out2_butt_35.^2 + cos_out2_butt_35.^2);
detection_26_65 = sqrt(sin_out2_butt_65.^2 + cos_out2_butt_65.^2);
 
% Построение графиков второй фильтрации
mid = 11.3; % Пороговое значение, взятое из lab2 для БФ 6ого порядка
figure;
phi = 0:Td:tend;
plot(phi, detection_26_10, t, detection_26_35, t, detection_26_65);
yline(mid,'-r','LineWidth',2);
title('Результат второго переноса частоты');
xlabel('t, сек');
ylabel('Уровень сигнала');
legend('Шум 10%', 'Шум 35%', 'Шум 65%', 'Пороговое значение, взятое из lab2 для БФ 6ого порядка');
 
% Определение наличия сигнала
lbound = mid;
sig_out6_10(1:N) = 1;
sig_out6_35(1:N) = 1;
sig_out6_65(1:N) = 1;
for k = 1:N
    if detection_26_10(k) <= lbound
        sig_out6_10(k) = 0;
    end
    if detection_26_35(k) <= lbound
        sig_out6_35(k) = 0;
    end
    if detection_26_65(k) <= lbound
        sig_out6_65(k) = 0;
    end
end
 
%  Построение графиков определения наличия сигнала
figure;
p = plot(t, sig_out6_10, "-b", t, sig_out6_35, "-r", t, sig_out6_65,  "-g");
p(1).LineWidth = 2;p(2).LineWidth = 2;p(3).LineWidth = 2;
title('Определение наличия сигнала');
xlabel('t, сек');
ylabel('Наличие сигнала (0 - нет, 1 - есть)');
legend('Шум 10%', 'Шум 35%', 'Шум 65%', 'Location', 'southeast');
axis([0.294 0.315 -0.1 1.1]);
 
% определение задержек определителя
k = 1;
while k <= length(sig_out6_10) && sig_out6_10(k) == 0
    k = k + 1;
end
delay_10=t(k);  % задержка определителя 6-го порядка, шум 10%
k = 1;
while(sig_out6_35(k) == 0)
    k = k + 1;
end
delay_35=t(k);  % задержка определителя 6-го порядка, шум 35%
k = 1;
while(sig_out6_65(k) == 0)
    k = k + 1;
end
delay_65=t(k);  % задержка определителя 6-го порядка, шум 65%
