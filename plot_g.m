function plot_g(R1,R2,C1,C2,fmin,fmax)
    %NOT WORKING AT ALL, TODO%
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    fmin = 1;
    fmax = 10000;
    F = @(x) 20*log10(abs((x/(C2*R2)) / (x^2+x*(1/(C1*R1)+1/(C1*R2)+1/(C2*R2)+1/(C1*C2*R1*R2)))));
    x = [fmin : 0.1 : fmax];
    for i = 1 : length(x)
        x(i) = j*2*pi*x(i);
    end
    for i = 1 : length(x)
        y(i) = F(x(i));
    end
    semilogx(y);
    grid on
end