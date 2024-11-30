% исходные данные
Fn = 360;    % Несущая частота
Fm = 10;     % частота модулируемого сигнала
Phin = 0;    % фаза несущей
Phim = 0;    % фаза модулируемого сигнала
m = 1;       % коэффициент модуляции 0 < m <= 1
I = 256;     % число уровней квантования
Kdiscr = 64;  % отношение частоты дискретизации к частоте несущей
tend = 0.05;  % время окончания модуляции сигнала
Fd = Fn*Kdiscr;   % частота дискретизации
Td = 1/Fd;   % период дискретизации
 
% Генерация сигнала
t = 0:Td:tend;
N = length(t);    % общее количество отсчетов
 
% Сгенерируем прямоугольный сигнал с разной скважностью
mod_sig_2 = generate_signal(Fn, 2, t, N);
mod_sig_4 = generate_signal(Fn, 4, t, N);
mod_sig_6 = generate_signal(Fn, 6, t, N);
 
% Построим графики
plot(t, mod_sig_2, t, mod_sig_4, t, mod_sig_6);
title('Модулированный сигнал');
xlabel('t, сек');
ylabel('Амплитуда');
legend('Скважность Q=2', 'Скважность Q=4', 'Скважность Q=6');
axis([-0.0005 0.01+0.0005 -1.2 1.2]);
 
% Дискретизация сигнала
y_2 = floor((mod_sig_2+m)*I/(2*m));
y_4 = floor((mod_sig_4+m)*I/(2*m));
y_6 = floor((mod_sig_6+m)*I/(2*m));
figure;
plot(t, y_2, t, y_4, t, y_6, '-');
title('Дискретизированный сигнал');
xlabel('t, сек');
ylabel('Уровень сигнала');
legend('Скважность Q=2', 'Скважность Q=4', 'Скважность Q=6');
axis([-0.0005 0.01+0.0005 -5 I+5]);
 
% первый перенос частоты - получим исходный сигнал, убрав несущую
detection_2 = first_shift_frequency(y_2, Fn, Fd, Fm, Kdiscr, N);
detection_4 = first_shift_frequency(y_4, Fn, Fd, Fm, Kdiscr, N);
detection_6 = first_shift_frequency(y_6, Fn, Fd, Fm, Kdiscr, N);
figure;
 
plot(t, detection_2, t, detection_4, t, detection_6);
title('Результат первого переноса частоты');
xlabel('t, сек');
ylabel('U, B');
legend('Скважность Q=2', 'Скважность Q=4', 'Скважность Q=6');
 
% второй перенос частоты
detection2 = second_shift_frequency(detection_2, Fd, Fm, Kdiscr, N);
detection4 = second_shift_frequency(detection_4, Fd, Fm, Kdiscr, N);
detection6 = second_shift_frequency(detection_6, Fd, Fm, Kdiscr, N);
 
cnt = length(t) / 4;
 
% Вычисление пороговое значения для разных скважностей
mid2 = sum(detection2(cnt:2*cnt));
mid4 = sum(detection4(cnt:2*cnt));
mid6 = sum(detection6(cnt:2*cnt));
 
mid2 = mid2 / cnt * 0.708;
mid4 = mid4 / cnt * 0.708;
mid6 = mid6 / cnt * 0.708;
 
figure;
plot(t, detection2, t, detection4, t, detection6);
yline(mid2,'-b','LineWidth',2);
yline(mid4,'-r','LineWidth',2);
yline(mid6,'-y','LineWidth',2);
title('Результат второго переноса частоты');
xlabel('t, сек');
ylabel('U, B');
legend('Скважность Q=2', 'Скважность Q=4', 'Скважность Q=6', 'Пороговое значение для скважности Q=2', 'Пороговое значение для скважности Q=4', 'Пороговое значение для скважности Q=6', 'Location', 'southeast');
 
% определение наличия сигнала 
sig_out2 = signal_check(mid2, detection2, N);
sig_out4 = signal_check(mid4, detection4, N);
sig_out6 = signal_check(mid6, detection6, N);
 
figure;
p = plot(t, sig_out2, "-b", t, sig_out4, "-r", t, sig_out6,  "-g");
p(1).LineWidth = 2;p(2).LineWidth = 2;p(3).LineWidth = 2;
title('Определение наличия сигнала');
xlabel('t, сек');
ylabel('Наличие сигнала (0 - нет, 1 - есть)');
legend('Скважность Q=2', 'Скважность Q=4', 'Скважность Q=6', 'Location', 'southeast');
axis([0.01 0.02 -0.1 1.1]);
 
% Определение задержки определителя
delay2 = find_delay(sig_out2, t);
delay4 = find_delay(sig_out4, t);
delay6 = find_delay(sig_out6, t);
 
% значения для ещё одной несущей частоты
Fn2 = Fn + 3 * Fm;
Phin2 = 0;
% Сложение сигналов
mod_sig2 = mod_sig_2 + 0.5*sin(2*pi*Fn2.*t+Phin2);
figure;
plot(t, mod_sig2);
title('Зашумленный сигнал');
xlabel('t, сек');
ylabel('U, B');
axis([0 tend -2 2]);
 
% анализ спектра сигнала - БПФ преобразование на fftL отсчетов
fftL = 2^nextpow2(N);
Y = abs(fft(mod_sig_2,fftL));       % Амплитуды преобразования Фурье сигнала
Y = 2*Y./N;                         % Нормировка спектра по амплитуде
Y(1) = Y(1)/2;                      % Нормировка постоянной составляющей в спектре
Y2 = abs(fft(mod_sig2,fftL));       % Амплитуды преобразования Фурье сигнала
Y2 = 2*Y2./N;                       % Нормировка спектра по амплитуде
Y2(1) = Y2(1)/2;                    % Нормировка постоянной составляющей в спектре
F=0:Fd/fftL:Fd/2-1/fftL;            % Массив частот вычисляемого спектра Фурье
 
 
% отобразим график только одной (левой) стороны спектра
figure;
plot(F, Y(1:length(F)));
title('Cпектр исходного сигнала')
xlabel('f, Гц')
ylabel('|Y(f)|')
axis([0 8000 0 1.3]);
figure;
plot(F, Y2(1:length(F)));
title('Cпектр зашумленного сигнала')
xlabel('f, Гц')
ylabel('|Yз(f)|')
axis([0 8000 0 1.3]);
 
% функция генерации прямоугольного сигнала
function [sig] = generate_signal(Fn, Q, t, N)
    sig = [];
    T = 1/Fn;
    Ti = T/Q;
    for i = 1:N
        phase = mod(t(i),T);
        if (phase < Ti)
            sig(end+1) = 1;
        else
            sig(end+1) = -1; 
        end
    end
end
 
% Функция для первого переноса частоты
function [detection] = first_shift_frequency(y, Fn, Fd, Fm, Kdiscr, N)
    % домножение на Sin и Cos
    for i = 1:N
        sin_out(i) = y(i)*sin((i-1) * 2*pi*Fn/Fd);
        cos_out(i) = y(i)*cos((i-1) * 2*pi*Fn/Fd);
    end
    % пропустим результаты домножения на Sin и Cos через БФ 2 порядка
    [b, a] = butter(6, Fm*Kdiscr/(2*Fd));
    sin_out_butt = filter(b, a, sin_out);
    cos_out_butt = filter(b, a, cos_out);
    % операция детектирования
    detection = sqrt(sin_out_butt.^2 + cos_out_butt.^2); 
end
 
% Функция для второго переноса частоты
function [detection] = second_shift_frequency(y, Fd, Fm, Kdiscr, N)
    % домножение на Sin и Cos
    for i = 1:N
        sin_out2(i) = y(i)*sin((i-1) * 2*pi*Fm/Fd);
        cos_out2(i) = y(i)*cos((i-1) * 2*pi*Fm/Fd);
    end 
 
    [b2, a2] = butter(6, Fm*Kdiscr/(2*8*Fd));
    sin_out2_butt = filter(b2, a2, sin_out2);
    cos_out2_butt = filter(b2, a2, cos_out2); 
    detection = sqrt(sin_out2_butt.^2 + cos_out2_butt.^2);
end
 
% Функция определения наличия сигнала
function[signal] = signal_check(mid, detection, N)
    lbound = mid; % нижний порог
    signal(1:N) = 1;
    for k = 1:N
      if detection(k) <= lbound
        signal(k) = 0;
      else
        break;
      end
    end
end
 
% Функция для определения задержки определителя
function[delay] = find_delay(signal, t)
    k = 1;
    while(signal(k) == 0)
        k = k + 1;
    end
    delay = t(k); % задержка определителя, выраженная в секундах
end
