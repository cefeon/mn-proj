function calka()
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e5;
    h = 1e-7;
    freq = 100000;
    t = [ 0 : h : 1 ]; 
    y = [0 0]';
    
    e = @(t) sin(2*pi*t*freq);
    dy = @(t,y) ...
        [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
           1/C2 * ( (e(t) - y(1) - y(2))/R2 + y(2)/RL ) ];
    u = euler(t, h, dy);

    u1 = u(1,:);
    u2 = u(2,:);
    
    for i=1 : length(t)
        dP(i) = (e(t(i))-u1(i))^2/R1 + (e(t(i))-u1(i)-u2(i))^2/R1;
    end
    
    
    %metoda prostokatów
    prostokat = dP(1:end-1) * h;
    prostokaty = sum(prostokat)
    
    %metoda parabol (Simpsona)
    for i = 1 : 2 : length(t)-2
        simpson((i + 1) / 2) = h/3*(dP(i)+4*dP(i+1)+dP(i+2));
    end
    parabole = sum(simpson)
end



function y = euler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
       y(:, i+1) = y(:, i) + h * f(t(i), y(:, i));
    end
end