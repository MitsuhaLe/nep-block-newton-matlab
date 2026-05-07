time = 0;
for iter = 1:10
    n=10000;
    A = rand(n);
    tic
    lu(A);
    time = time + toc;
end
avgTime = time / 10;