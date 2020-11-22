X = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
Y = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];

x = 1;
xx=[-1:0.01:1];
for i=1:length(xx)
    y(i) = splineit(xx(i),X,Y);
end
plot(X,Y,'o',xx,y)

%znajduje do którego sektora należy x
function sect = sector(x,X)
    N = length(X);
    R = N;
    L = 1;
    while 1
        if(R-L) <= 1
            sect = L;
            break;
        end

        sect = fix((L+R)/2);

        if x < X(sect)
            R = sect;
        else
            L = sect;
        end
    end
end
%oblicza macierz drugich pochodnych
function diff2 = difmatrix(X,Y,h)
    N = length(X);
    M = diag(ones(1, N - 3), 1) + diag(4 * ones(1, N - 2))...
      + diag(ones(1, N - 3), - 1);

    for i = 1:N - 2
        z(i) = (6 / h ^ 2) * (Y(i + 2) - 2 * Y(i + 1) + Y(i));
    end

    diff2 = [0; M \ z'; 0];
end
%oblicza wartość spline'a w punkcie x
function y = splineit(X,Y,x)
    p = sector(x,X);
    h = X(p)-X(p+1);
    diff2 = difmatrix(X,Y,h);
    y = ((x - X(p+1))^3/h - (x - X(p+1))*h)*diff2(p)/6 ...
        - ((x - X(p))^3/h - (x - X(p))*h)*diff2(p+1)/6 ...
        + Y(p)*(x - X(p+1))/h - Y(p+1)*(x - X(p))/h;
end

        

