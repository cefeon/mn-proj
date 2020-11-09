function fc = find_fc(R1,R2,C1,C2)
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    F = @(x) (0.7071*0.2727)... 
        - abs(((2*pi*x*1i)/((C2)*R2)) / (((2*pi*x*1i)^2)...
        + (2*pi*x*1i)*(1/((C1)*R1)+1/((C1)*R2)+1/((C2)*R2))...
        + 1/((C1)*(C2)*R1*R2)));
    fprintf("Bisekcja | dolna ");
    freq_bisf(F,2,25000); 
    fprintf("Bisekcja | górna ");
    freq_bisf(F,300000,500000); 
    fprintf("Sieczne | dolna ");
    freq_secf(F,2,25000); 
    fprintf("Sieczne | górna ");
    freq_secf(F,300000,500000); 
    fprintf("quasi-Newtona | dolna ");
    freq_qn(F,10000,0.5e-11); 
    fprintf("quasi-Newtona | górna");
    freq_qn(F,300000,1e-9); 
end
function x = freq_bisf(F,x0,x1)
    iterator = 0;  
	x = (x0+x1)/2;
	precision = 1e-6;
	while abs(F(x))>precision
        if (F(x)*F(x0) < 0)
            x1=x;
        else
            x0=x;
        end
        x=(x0+x1)/2;
        iterator = iterator+1;
    end
	fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

function x = freq_secf(F,x0,x1)
    iterator = 0;
    x = x1-F(x1)*(x1-x0)/(F(x1)-F(x0));
    precision = 1e-6;
    while abs(F(x))>precision
        x = x1-F(x1)*(x1-x0)/(F(x1)-F(x0));
        x0 = x1;
        x1 = x;
        iterator = iterator+1;
    end
	fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

function x = freq_qn(F,x,delta)
    precision = 1e-6;
    iterator = 0;
    while abs(F(x))>precision
        x = x-F(x)*delta/(F(x+delta)-F(x));
        iterator = iterator+1;
    end
	fprintf("granica to %f. %d iteracji użyto \n",x,iterator);
end

