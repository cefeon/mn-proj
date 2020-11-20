function plot_g(R1,R2,C1,C2,fmin,fmax)
    F = @(x) 20*log10 (abs((x/(C2*R2)) / (x^2+x*(1/(C1*R1)+1/(C1*R2)+1/(C2*R2))+1/(C1*C2*R1*R2))));
    x = logspace (log10(fmin),log10(fmax));
    z = 2i*pi*x;
    
    y = zeros(1,length(x));
    for i = 1 : length(x)
        y(i) = F(z(i));
    end
    
    thr = zeros(1,length(x));
    for i = 1 : length(x)
        thr(i) = max(y)-3;
    end
    
    semilogx(x,y,x,thr);
    grid on
end