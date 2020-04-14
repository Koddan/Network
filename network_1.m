clear all;

g = [1 0 1 1];
l = 4;
r = length(g) - 1;
n = l + r;
messages = createSet(l);
codes = zeros(2^l, n);
for i = 1:2^l
    codes(i,:) = encode(messages(i,:), r, g);
end
d = 0;
for i = 1:n
    d = d + codes(2,i);
end
for i = 3:2^l
    tmp = 0;
    for j = 1:n
        tmp = tmp + codes(i,j);
    end
    if (tmp < d)
        d = tmp;
    end
end
disp(d);

p = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
est = upperEstimate(p, d, n);
rightEst = rightEstimate(p, d, n, codes);

figure(1);
hold on;
grid on;
plot(p, est);
plot(p, rightEst);
xlabel('p');
ylabel('P_e');
legend('P_e^+','P_e');
hold off;

% c)

L = [2 3 4 5 6 7 8 9 10 11 12];
[estByK1, D1] = upperEstimateByL(L, 0.1, g);
[estByK2, D2] = upperEstimateByL(L, 0.2, g);
[estByK3, D3] = upperEstimateByL(L, 0.3, g);

[rightestByK1, D4] = rightEstimateByL(L, 0.1, g);
[rightestByK2, D5] = rightEstimateByL(L, 0.2, g);
[rightestByK3, D6] = rightEstimateByL(L, 0.3, g);

figure(2);
hold on;
grid on;
plot(L, estByK1, 'b');
plot(L, estByK2, 'b--');
plot(L, estByK3, 'b-.');
plot(L, rightestByK1, 'r');
plot(L, rightestByK2, 'r--');
plot(L, rightestByK3, 'r-.');
xlabel('l');
ylabel('P_e');
legend('P_e, p = 0.1', 'P_e, p = 0.2', 'P_e, p = 0.3', 'P_e^+, p = 0.1', 'P_e^+, p = 0.1', 'P_e^+, p = 0.1');
hold off;

figure(3);
hold on;
grid on;
plot(L, D1, 'b');
xlabel('l');
ylabel('d');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% функции

function polynom = product(polynom1, polynom2)
    matrix = zeros(length(polynom1), ((length(polynom1) - 1) + (length(polynom2) - 1)) + 1);
    for i = 1:length(polynom1)
        for j = 1:length(polynom2)
            matrix(i,(((i - 1) + (j - 1)) + 1)) = polynom1(i)*polynom2(j);
        end
    end
    for i = 1:((length(polynom1) - 1) + (length(polynom2) - 1)) + 1
        for j = 2:length(polynom1)
            matrix(1,i) = matrix(1,i) + matrix(j,i);
        end
    end
    polynom = matrix(1,:);
end

function polynom = division(polynom1, polynom2)
    for i = 1:length(polynom1)
        if ((length(polynom1) - (i - 1)) - length(polynom2)) >= 0
            degree = ((length(polynom1) - (i - 1)) - length(polynom2)) + 1;
            tmpPolynom = zeros(degree, 1);
            tmpPolynom(1) = polynom1(i);
            substraction = zeros(length(polynom1), 1);
            tmpSubstraction = product(polynom2, tmpPolynom);
            counter = 1;
            for k = (length(polynom1) - length(tmpSubstraction) + 1):length(polynom1)
                substraction(k) = tmpSubstraction(counter);
                counter = counter + 1;
            end
            
            for q = 1:length(substraction)
                substraction(q) = rem(substraction(q), 2);
                substraction(q) = substraction(q)^2;
            end
            for j = 1:length(polynom1)
                polynom1(j) = polynom1(j) - substraction(j);
            end
            for q = 1:length(polynom1)
                polynom1(q) = rem(polynom1(q), 2);
                polynom1(q) = polynom1(q)^2;
            end
        end
    end
    polynom = polynom1;
    
end

function set = createSet(length)
    set = zeros(2^length, length);
    for i = 1:2^length
        number = i - 1;
        for j = 1:length
            set(i, length - j + 1) = rem(number, 2);
            number = fix(number/2);
        end
    end
end

function a = encode(message, r, g)
    x = zeros(r + 1, 1);
    x(1) = 1;
    tmp = product(message, x);
    c = division(tmp, g);
    a = plus(tmp, c);
end

function est = upperEstimate(p, d, n)
    est = zeros(length(p), 1);
    for q = 1:length(p)
        sum = 0;
        for i = 0:d - 1
            C = factorial(n)/(factorial(i)*factorial(n - i));
            sum = sum + C*p(q)^i*(1 - p(q))^(n - i);
        end
        est(q) = 1 - sum;
    end
end

function [est, D] = upperEstimateByL(L, p, g)
    est = zeros(length(L), 1);
    counter = 1;
    D = zeros(length(L), 1);
    for h = 1:length(L)
        k = L(h);
        r = length(g) - 1;
        n = k + r;
        messages = createSet(k);
        codes = zeros(2^k, n);
        for i = 1:2^k
            codes(i,:) = encode(messages(i,:), r, g);
        end
        d = 0;
        for i = 1:n
            d = d + codes(2,i);
        end
        for i = 3:2^k
            tmp = 0;
            for j = 1:n
                tmp = tmp + codes(i,j);
            end
            if (tmp < d)
                d = tmp;
            end
        end
        D(h) = d;
        est(counter) = upperEstimate(p, d, n);
        counter = counter + 1;
    end
end

function est = rightEstimate(p, d, n, codes)
    est = zeros(length(p), 1);
    for q = 1:length(p)
        sum = 0;
        for i = d:n
            A = countByWeight(codes, n, i);
            sum = sum + A*p(q)^i*(1 - p(q))^(n - i);
        end
        est(q) = sum;
    end
end

function [est, D] = rightEstimateByL(L, p, g)
    est = zeros(length(L), 1);
    counter = 1;
    D = zeros(length(L), 1);
    for h = 1:length(L)
        k = L(h);
        r = length(g) - 1;
        n = k + r;
        messages = createSet(k);
        codes = zeros(2^k, n);
        for i = 1:2^k
            codes(i,:) = encode(messages(i,:), r, g);
        end
        d = 0;
        for i = 1:n
            d = d + codes(2,i);
        end
        for i = 3:2^k
            tmp = 0;
            for j = 1:n
                tmp = tmp + codes(i,j);
            end
            if (tmp < d)
                d = tmp;
            end
        end
        D(h) = d;
        est(counter) = rightEstimate(p, d, n, codes);
        counter = counter + 1;
    end
end

function number = countByWeight(codes, n, weight)
    number = 0;
    for i = 1:size(codes)
        tmp = 0;
        for j = 1:n
            tmp = tmp + codes(i,j);
        end
        if (tmp == weight)
            number = number + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%