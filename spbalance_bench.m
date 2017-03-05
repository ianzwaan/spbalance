function spbalance_bench
    % SPBALANCE_BENCH  benchmark spbalance and mexspbalance.
    %
    %   Copyright Ian Zwaan <i.n.zwaan@tue.nl> 2013 - 2016.  Distributed
    %   under the Boost Software License, Version 1.0; see
    %   mexspbalance.c or copy at http://www.boost.org/LICENSE_1_0.txt.

    N = 10000;

    fprintf('Benchmark with %d iterations.\n', N);

    load mat/tols4000;

    tic;
    for i = 1:N
        [T, B] = spbalance(A);
    end
    t = toc;

    fprintf('\nspbalance(tols4000):\n');
    fprintf('    elapsed = %g\n', t);
    fprintf('    average = %g\n', t / N);
    fprintf('    bal/sec = %g\n', N / t);

    tic;
    for i = 1:N
        [T, B] = mexspbalance(A);
    end
    t = toc;

    fprintf('\nmexspbalance(tols4000):\n');
    fprintf('    elapsed = %g\n', t);
    fprintf('    average = %g\n', t / N);
    fprintf('    bal/sec = %g\n', N / t);
end
