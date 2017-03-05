function spbalance_demo
    % SPBALANCE_DEMO  demonstrate how to use spbalance.
    %
    %   Copyright Ian Zwaan <i.n.zwaan@tue.nl> 2013 - 2016.  Distributed
    %   under the Boost Software License, Version 1.0; see
    %   mexspbalance.c or copy at http://www.boost.org/LICENSE_1_0.txt.

    load mat/tols4000;

    problem_size = size(A, 1)

    tic
    [T, B] = spbalance(A);
    elapsed_time = toc

    subplot(131); spy(A); title('A');
    subplot(132); spy(T); title('T');
    subplot(133); spy(B); title('B');

    similarity_error = norm(B - T\A*T, 'fro')

    S = full(T(find(T)));
    scale_min = min(S)
    scale_max = max(S)
end
