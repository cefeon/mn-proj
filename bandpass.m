function bandpass(freq)
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e8;
    h = 1e-7;
    t = [ 0 : h : 0.1e-3]; 

    if freq>0
        e = @(t) sin(2*pi*t*freq);
    end
    
    if freq=='nosin'
        e = @(t) 1;
    end

    if freq=='cycle'    
       e = @(t) rectpulse(t,0.05e-3);
    end

    dy = @(t,y) ...
        [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
           1/C2 * ( (e(t) - y(1) - y(2))/R2 - y(2)/RL ) ];
    
    %u = euler(t, h, dy); %metoda Eulera
    u = beuler(t, h, dy); %ulepszona metoda Eulera
    plot(t, u(1,:)); %rysuj wykres u1
    %plot(t, u(2,:)); %rysuj wykres u2
    xlim([0 1e-4])
    grid on
end

%metoda Eulera
function y = euler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
       y(:, i+1) = y(:, i) + h * f(t(i), y(:, i));
    end
end

%ulepszona metoda Eulera
function y = beuler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
        prediction = y(:, i) + h/2 * f(t(i), y(:, i));
        y(:, i+1) = y(:, i) + h * f(t(i) +h/2, prediction);
    end
end

%funkcja okresowa kwad
function y = rectpulse(x,T)
    modulo = mod(x,T);
    if modulo<(T/2)
        y = 1;
    else
        y = 0;
    end
end