un = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
in = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];
i = @(u) interpolate(u,un,in);
r = @(u) i(u)/u;

%interpolacja wielomianowa m. Vandermonde'a
function i = interpolate(u,un,in)
    n = length(un);
    A = zeros(n);

    for w=1:n
        for k=1:n
            A(w,k) = un(w)^(k-1);
        end
    end
    a = A\in';

    s = 0;
    for i=1:n
        s = s+a(i)*u^(i-1);
    end
    i = s;
end