function nonlinear()
    R1 = 30;
    R2 = 10;
    C1 = 0.1e-6;
    C2 = 0.2e-6;
    h = 1e-7;
    t = [ 0 : h : 0.01e-3]; 
    freq = 65000;
    
    e = @(t) sin(2*pi*t*freq);

    un = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
    in = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];

    c_aprox3 = aprox(un, in, 3);
    c_aprox5 = aprox(un, in, 5); 
    vander = vandermonde(un);
    diff2 = difmatrix(un, in, 0.25);
    i_inter = @(u) ctv(interpolate(in, vander), u);
    i_aprox3 = @(u) ctv(c_aprox3, u);
    i_aprox5 = @(u) ctv(c_aprox5, u);
    i_spline = @(u) splineit(un, in, diff2, u);

    dy = @(t,y,In) ...
            [  1/C1 * ( (e(t) - y(1) - y(2))/R2 + (e(t) - y(1))/R1 )
               1/C2 * ( (e(t) - y(1) - y(2))/R2 - In )];

    u = euler(t, h, dy, i_spline);

    p_rn = zeros(1,length(t));

    for v=1 : length(t)
        p_rn(v) = u(2,v)*i_spline(u(2,v));
    end
    
    p = int_simps(t,h,p_rn);
    fprintf('Moc wydzielana na rezystorze to %e', p);
    plot(t, u(2,:));
    xlim([0 1e-5])
end

%interpolacja
function A = interpolate(Y, M)
    A = M \ Y';
end

%aproksymacja stopnia n
function A = aprox(X, Y, n)
    N = length(X);
    M = zeros(N,n-1);
    for w = 1:N
        for k = 1:n+1
            M(w, k) = X(w) ^ (k - 1);
        end
    end
    MTM = M'*M;
    MTY = M'*Y';
    A = MTM \ MTY;
end

%macierz Vandermonde'a
function M = vandermonde(X)
    N = length(X);
    for w = 1:N
        for k = 1:N
            A(w, k) = X(w) ^ (k - 1);
        end
    end
    M = A;
end

%na podstawie macierzy współczynników i x wylicza wartość wielomianu W(x)
function y = ctv(A, x)
    s = 0;
    for i = 1:length(A)
        s = s + A(i) * x ^ (i - 1);
    end
    y = s;
end

%znajduje do którego sektora należy x do spline'a
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

%oblicza macierz drugich pochodnych do spline'a
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
function y = splineit(X,Y,diff2,x)
    p = sector(x,X);
    h = X(p)-X(p+1);
    y = ((x - X(p+1))^3/h - (x - X(p+1))*h)*diff2(p)/6 ...
        - ((x - X(p))^3/h - (x - X(p))*h)*diff2(p+1)/6 ...
        + Y(p)*(x - X(p+1))/h - Y(p+1)*(x - X(p))/h;
end

%metoda Eulera
function y = euler(t,h,f,In)
    y = [0 0]';
    for i = 1 : length(t)-1
       y(:, i+1) = y(:, i) + h * f(t(i), y(:, i), In(y(2,i)));
    end
end

%złożona metoda parabol (Simpsona)
function calka = int_simps (t,h,df)
    simpson = zeros(1,(length(t)+1)/2);
    for i = 1 : 2 : length(t)-2
        simpson((i + 1) / 2) = h/3*(df(i)+4*df(i+1)+df(i+2));
    end
    calka = sum(simpson);
end
