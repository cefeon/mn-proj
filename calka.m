function calka(freq)
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e8;
    h = 1e-8; %tutaj zmienia się dt
    t = [ 0 : h : 1]; 
    
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

    u = euler(t, h, dy);
    
    %obliczanie wartości dP ze wzoru
    dP = zeros(1,length(t));
    for i=1 : length(t)
        dP(i) = (e(t(i))-u(1,i))^2/R1 + (e(t(i))-u(1,i)-u(2,i))^2/R2;
    end
    parabole = int_simps (t,h,dP);
    prostokaty = int_rect (t,h,dP);
    
    fprintf('Metoda prostokatów: %e \n\n',prostokaty);
    fprintf('Metoda Simpsona: %e \n\n',parabole);
end

%złożona metoda parabol (Simpsona)
function calka = int_simps (t,h,df)
    simpson = zeros(1,(length(t)+1)/2);
    for i = 1 : 2 : length(t)-2
        simpson((i + 1) / 2) = h/3*(df(i)+4*df(i+1)+df(i+2));
    end
    calka = sum(simpson);
end

%złożona metoda prostokątów lewych
function calka = int_rect (t,h,df)
    dfdx = zeros(1,length(t)-1);
    for i =  1 : length(t)-1
        dfdx(i) = df(i) * h;
    end
    calka = sum(dfdx);
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

%funkcja okresowa kwadratowokształtna
function y = rectpulse(x,T)
    modulo = mod(x,T);
    if modulo<(T/2)
        y = 1;
    else
        y = 0;
    end
end