function calka()
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e8;
    h = 1e-7;
    t = [ 0 : h : 1e-4]; 
    freq = 50;
    
    e = @(t) sin(2*pi*t*freq);
    dy = @(t,y) ...
        [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
           1/C2 * ( (e(t) - y(1) - y(2))/R2 - y(2)/RL ) ];
    u = beuler(t, h, dy);
    
    %obliczanie wartości P ze wzoru
    for i=1 : length(t)
        dP(i) = (e(t(i))-u(1,i))^2/R1 + (e(t(i))-u(1,i)-u(2,i))^2/R2;
    end
    
    plot(t,dP);
    
    %złożona metoda prostokątów lewych
    prostokat = dP(1:end-1) * h;
    prostokaty = sum(prostokat);
    
    %złożona metoda parabol (Simpsona)
    for i = 1 : 2 : length(t)-2
        simpson((i + 1) / 2) = h/3*(dP(i)+4*dP(i+1)+dP(i+2));
    end
    parabole = sum(simpson);
    
    fprintf('Metoda prostokatów: %e \n\n',prostokaty);
    fprintf('Metoda Simpsona: %e \n\n',parabole);
end



function y = euler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
       y(:, i+1) = y(:, i) + h * f(t(i), y(:, i));
    end
end

function y = beuler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
        prediction = y(:, i) + h/2 * f(t(i), y(:, i));
        y(:, i+1) = y(:, i) + h * f(t(i) +h/2, prediction);
    end
end