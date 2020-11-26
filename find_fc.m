function fc = find_fc(R1, R2, C1, C2)
    F = @(x) (0.7071 * 0.2727)...
       - abs(((2i * pi * x) * (1/(C2 * R2))) / (((2i * pi * x) ^ 2)...
       + (2i * pi * x) * (1 / (C1 * R1) + 1 / (C1 * R2) + 1 / (C2 * R2))...
       + 1 / (C1 * C2 * R1 * R2)));
    precision = 1e-7;
    a_dolna = 0;
    b_dolna = 35000;
    a_gorna = 300000;
    b_gorna = 500000;
    fprintf("Metoda bisekcji | dolna ");
    freq_bisf(F, a_dolna, b_dolna, precision);
    fprintf("Metoda bisekcji | górna ");
    freq_bisf(F, a_gorna, b_gorna, precision);
    fprintf("Metoda siecznych | dolna ");
    freq_secf(F, a_dolna, b_dolna, precision);
    fprintf("Metoda siecznych | górna ");
    freq_secf(F, b_gorna, b_gorna, precision);
    fprintf("Metoda quasi-Newtona | dolna ");
    freq_qn(F, a_dolna, b_dolna, 1e-7, precision);
    fprintf("Metoda quasi-Newtona | górna ");
    freq_qn(F, a_gorna, b_gorna, 1e-7, precision);

end

%metoda bisekcji
function x = freq_bisf(F, x0, x1, precision)
    iterator = 0;
    x = (x0 + x1) / 2;
    while abs(F(x)-F(x1)) > precision
        if (F(x) * F(x0) < 0)
            x1 = x;
        else
            x0 = x;
        end
        x = (x0 + x1) / 2;
        iterator = iterator + 1;
    end
    fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

%metoda siecznych
function x = freq_secf(F, x0, x1, precision)
    iterator = 0;
    x = x1 - F(x1) * (x1 - x0) / (F(x1) - F(x0));
    while abs(F(x)) > precision
        x = x1 - F(x1) * (x1 - x0) / (F(x1) - F(x0));
        x0 = x1;
        x1 = x;
        iterator = iterator + 1;
    end
    fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

%metoda quasi-Newtona
function x = freq_qn(F, x0,x1, delta, precision)
    dif = differ(F,x0,delta);
    if(F(x0)*differ(F, dif,delta)>0)
        x = x0;
    else
        x = x1;
    end
    iterator = 0;
    while abs(F(x)) > precision
        x = x - F(x) / differ(F, x, delta);
        iterator = iterator + 1;
    end
    fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

%obliczanie pochodnej 
function dif = differ (F, x, delta)
    dif = (F(x + delta) - F(x))/delta;
end