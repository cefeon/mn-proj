	
function plot_g(R1, R2, C1, C2, fmin, fmax)
    F = @(x) 20 * log10 (abs((x / (C2 * R2)) / ...
         (x ^ 2 + x * (1 / (C1 * R1) + 1 / (C1 * R2) ...
         + 1 / (C2 * R2)) + 1 / (C1 * C2 * R1 * R2))));
    x = logspace (log10(fmin), log10(fmax));
    z = 2i * pi * x;

    y = zeros(1, length(x));
    for i = 1 : length(x)
        y(i) = F(z(i));
    end
    fc1=13811.801114;
    fc2=305659.838834;
    thr = max(y) - 3;
    
    semilogx(x, y);
    ylabel('Amplituda[dB]');
    xlabel('f[Hz]');
    hold on
    yline(thr,'-.','-3db');
    xline(fc1,'--','fc1');
    xline(fc2,'--','fc2');
    semilogx(fc1,thr,'or');
    semilogx(fc2,thr,'or');
    grid on
end