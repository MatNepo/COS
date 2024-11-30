% Исходные данные
Fm = 10;     % частота модулируемого сигнала
Fn = 360;    % частота несущей
m = 1;       % коэффициент модуляции 0 < m <= 1
Phin = 0;    % фаза модулирующего сигнала
Phim = 0;    % фаза модулируемого сигнала
I = 256;     % число уровней квантования
Kdiscr = 8;  % отношение частоты дискретизации к частоте несущей
Fd = Fn * Kdiscr;   % частота дискретизации
Td = 1 / Fd;   % период дискретизации
tend = 2.0;      % время окончания модуляции сигнала
 
% Генерация временного интервала
t = 0:Td:tend;
N = length(t);  % общее количество отсчетов

% Генерация сигнала
mod_sig = zeros(1, N);
for i = 1:N
    time = t(i);
    if time <= 6 / Fm
        % До 6-го периода частота и амплитуда исходные
        mod_sig(i) = 0.5 * (1 + m * sin(2 * pi * Fm * time + Phim)) * sin(2 * pi * Fn * time + Phin);
    elseif time > 6 / Fm && time <= 9 / Fm
        % 7–9 периоды: скачок частоты модуляции в 2 раза, амплитуда в 4 раза меньше
        mod_sig(i) = 0.125 * (1 + m * sin(2 * pi * (2 * Fm) * time + Phim)) * sin(2 * pi * Fn * time + Phin);
    else
        % После 9-го периода: параметры возвращаются к исходным
        mod_sig(i) = 0.5 * (1 + m * sin(2 * pi * Fm * time + Phim)) * sin(2 * pi * Fn * time + Phin);
    end
end
 
% Построение графика модулированного сигнала
figure;
plot(t, mod_sig);
title('Модулированный сигнал');
xlabel('t, сек');
ylabel('U, B');
axis([0 tend -1 1]);
 
% АЦП (Дискретизация сигнала)
y = floor((mod_sig + m) * I / (2 * m));
 
% Построение графика дискретизированного сигнала
figure;
plot(t, y);
title('Сигнал после АЦП');
xlabel('t, сек');
ylabel('Квантованные значения сигнала');
axis([0 tend 0 I]);

% Первый перенос частоты
sin_out = y .* sin((0:N-1) * 2 * pi * Fn / Fd);
cos_out = y .* cos((0:N-1) * 2 * pi * Fn / Fd);

% Пропускаем через фильтр Буттерворта
Fc1 = Fm * Kdiscr / (2 * Fd);
Fc1 = min(Fc1, 1); % Ограничиваем до допустимого диапазона [0, 1]
[b, a] = butter(4, Fc1);

sin_out_butt = filter(b, a, sin_out);
cos_out_butt = filter(b, a, cos_out);

% Операция детектирования
detection4 = sqrt(sin_out_butt.^2 + cos_out_butt.^2);

% Построение графика первой фильтрации
figure;
plot(t, detection4);
title('Результат первого переноса частоты');
xlabel('t, сек');
ylabel('Амплитуда');

% Второй перенос частоты
sin_out2 = detection4 .* sin((0:N-1) * 2 * pi * Fm / Fd);
cos_out2 = detection4 .* cos((0:N-1) * 2 * pi * Fm / Fd);

% Пропускаем через фильтр Буттерворта
Fc2 = Fm * Kdiscr / (2 * 8 * Fd);
Fc2 = min(Fc2, 1); % Ограничиваем до допустимого диапазона [0, 1]
[b, a] = butter(4, Fc2);

sin_out2_butt = filter(b, a, sin_out2);
cos_out2_butt = filter(b, a, cos_out2);

% Операция детектирования
detection4_2 = sqrt(sin_out2_butt.^2 + cos_out2_butt.^2);

% Построение графика второй фильтрации
mid = max(detection4_2) * 0.5; % Динамическое пороговое значение
figure;
plot(t, detection4_2);
yline(mid, '-r', 'LineWidth', 2);
title('Результат второго переноса частоты');
xlabel('t, сек');
ylabel('Уровень сигнала');
legend('Фильтр 4 порядка', 'Пороговое значение');

% Определение наличия сигнала
sig_out = detection4_2 > mid;

% Построение графика определения наличия сигнала
figure;
p = plot(t, sig_out, "-b");
p(1).LineWidth = 2;
title('Определение наличия сигнала');
xlabel('t, сек');
ylabel('Наличие сигнала (0 - нет, 1 - есть)');
axis([0 tend 0 1]);

% Определение задержки определителя
k = find(sig_out, 1); % Находим первый ненулевой элемент
delay = t(k);         % Задержка определителя 4-го порядка
disp(['Задержка определителя: ', num2str(delay), ' сек']);
