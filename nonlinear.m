un = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
in = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];

i = @(u) interpolate(u, un, in);
r = @(u) i(u) / u;

ud = [-1:0.01:1];
id = [];
ad3 = [];
ad5 = [];
spl = [];
for i = 1:length(ud)
    id(i) = ctv(interpolate(un, in), ud(i));
    ad3(i) = ctv(aprox(un, in, 3), ud(i));
    ad5(i) = ctv(aprox(un, in, 5), ud(i));
    hold on
end
    spl = spliner(un,in);
    plot(spl)
plot(un, in, 'o', ud, id, ud, ad3, ud, ad5)
grid on

%interpolacja
function A = interpolate(X, Y)
    M = vandermonde(X);
    A = M \ Y';
end

%aproksymacja stopnia n
function A = aprox(X, Y, n)
    M = vandermonde(X);
    M = M(:, 1:n);
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

%splajn 3 stopnia
function W = spliner(X, Y)
    N = length(X);
    h = (X(N) - X(1)) / (N - 1);
    M = diag(ones(1, N - 3), 1) + diag(4 * ones(1, N - 2))...
      + diag(ones(1, N - 3), - 1);
    
    for i = 1:N - 2
        z(i) = (6 / h ^ 2) * (Y(i + 2) - 2 * Y(i + 1) + Y(i));
    end
    w = M \ z';

    sigma = [0; w; 0];

    for i = 1 : N - 1
        a(i) = (sigma(i + 1) - sigma(i)) / (6 * h);
        b(i) = sigma(i) / 2;
        c(i) = (Y(i + 1) - Y(i)) / h ...
             - h / 6 * (2 * sigma(i) + sigma(i + 1));
        d(i) = Y(i);
    end

    r = 21;
    hh = h / r; %size of subinterval

    x = X(1):hh:X(N);

    W = [];
    for i = 1 : N - 1
        for j = r * (i - 1) + 1:r * i
            W(j) = a(i) * (x(j) - X(i)) ^ 3 ...
                 + b(i) * (x(j) - X(i)) ^ 2 ...
                 + c(i) * (x(j) - X(i)) + d(i);
        end
    end
    W (r * (N - 1) + 1) = Y(N);
end