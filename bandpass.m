function bandpass(freq)
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    RL = 1e5;
    h = 1e-7;
    t = [ 0 : h : 1]; 
    
    if (~isempty(freq)) 
        e = @(t) sin(2*pi*t*freq);
    else
        e = @(t) 1;
    end
    
    dy = @(t,y) ...
        [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
           1/C2 * ( (e(t) - y(1) - y(2))/R2 + y(2)/RL ) ];
    
    eulerV = euler(t, h, dy);
    plot(t, eulerV(2,:));
    xlim([0 1e-4])
    grid on
end

function y = euler(t,h,f)
    y = [0 0]';
    for i = 1 : length(t)-1
        y(:, i+1) = y(:, i) + h * f(t(i), y(:, i));
    end
end