function calka()
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e5;
    h = 1e-7;
    freq = 100000;
    t = [ 0 : h : 0.1e-5 ]; 
    y = [0 0]';
    e = @(t) sin(2*pi*t*freq);
    
    dy = @(t,y) ...
        [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
           1/C2 * ( (e(t) - y(1) - y(2))/R2 + y(2)/RL ) ];
    results = euler(t, h, dy);
    
    u1 = results(1,:);
    u2 = results(2,:);
    for i=1 : length(t)
        dP(i) = (e(t(i))-u1(i))^2/R1 + (e(t(i))-u1(i)-u2(i))^2/R1;
    end
    prostokat=dP*h;
    prostokaty = sum(prostokat)
    for a=1 : length(t)-2
        for i=1 : 3
            simpson(i)=h/3*(dP(a)+4*dP(a+1)+dP(a+2));
        end
    end
    parabole = sum(simpson)
end


function y = euler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
        y(:, i+1) = y(:, i) + h * f(t(i), y(:, i));
    end
end