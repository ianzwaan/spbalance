function spbalance_check
    % SPBALANCE_CHECK  check spbalance.
    %
    %   Copyright Ian Zwaan <i.n.zwaan@tue.nl> 2013 - 2016.  Distributed
    %   under the Boost Software License, Version 1.0; see
    %   mexspbalance.c or copy at http://www.boost.org/LICENSE_1_0.txt.

    rng(5489, 'twister');

    fibs = [0, 1 1 2 3 5 8 13 21 34 55 89 144 233 377 610 987];
    powers2 = [1 2 4 8 16 32 64 128 256 512 1024];
    powers3 = [1 3 9 27 81 243 729];
    powers5 = [1 5 25 125 625];
    powers7 = [1 7 49 343];
    powers10 = [1 10 100 1000];
    primes = [                                              ...
           2    3    5    7   11   13   17   19   23   29   ...
          31   37   41   43   47   53   59   61   67   71   ...
          73   79   83   89   97  101  103  107  109  113   ...
         127  131  137  139  149  151  157  163  167  173   ...
         179  181  191  193  197  199  211  223  227  229   ...
         233  239  241  251  257  263  269  271  277  281   ...
         283  293  307  311  313  317  331  337  347  349   ...
         353  359  367  373  379  383  389  397  401  409   ...
         419  421  431  433  439  443  449  457  461  463   ...
         467  479  487  491  499  503  509  521  523  541   ...
         547  557  563  569  571  577  587  593  599  601   ...
         607  613  617  619  631  641  643  647  653  659   ...
         661  673  677  683  691  701  709  719  727  733   ...
         739  743  751  757  761  769  773  787  797  809   ...
         811  821  823  827  829  839  853  857  859  863   ...
         877  881  883  887  907  911  919  929  937  941   ...
         947  953  967  971  977  983  991  997 1009 1013   ...
        1019 1021                                           ...
    ];

    sizes = [fibs powers2 powers3 powers5 powers7 powers10 primes];
    sizes = unique(sort(sizes));

    test_total = 0;
    test_passed = 0;
    test_failed = 0;

    fprintf('==============================================\n');
    fprintf('Test random permutations of diagonal matrices.\n');
    fprintf('==============================================\n');
    for s = 1:numel(sizes)
        N = sizes(s);
        Tperm = randperm(N);
        Aorig = speye(N);
        Aperm = Aorig(Tperm,Tperm);

        [T, B] = spbalance(Aperm);

        test_total = test_total + 1;
        if norm(Aorig - B, 'fro') <= eps
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', N);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', N);
        end
    end

    fprintf('\n');
    fprintf('======================================\n');
    fprintf('Test permuted block diagonal matrices.\n');
    fprintf('======================================\n');
    for s = 2:numel(sizes)
        N = max(sizes);
        M = sizes(s);
        if N / 2 < M
            break;
        end

        B = cell(floor(N/M), 1);
        for i = 1:M:N
            B{i} = sparse(unbalanced(rand(M)));
        end

        Aorig = blkdiag(B{:});
        Tperm = randperm(size(Aorig, 1));
        Aperm = Aorig(Tperm,Tperm);

        [T, B] = balance(full(Aorig));
        [TT, BB] = spbalance(Aperm);

        % NOTE: ordering of blocks can be different.
        simOk = norm(BB - TT\Aperm*TT, 'fro') <= eps;
        permOk = norm(spones(Aorig) - spones(BB), 'fro') <= eps;

        test_total = test_total + 1;
        if simOk && permOk
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', M);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', M);
        end
    end

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Test permuted upper triangular matrices.\n');
    fprintf('========================================\n');
    for s = 1:numel(sizes)
        N = sizes(s);
        Tperm = randperm(N);
        Aorig = unbalanced(sparse(triu(rand(N))));
        Aperm = Aorig(Tperm,Tperm);

        [T, B] = balance(full(Aorig));
        [TT, BB] = spbalance(Aperm);

        test_total = test_total + 1;
        if check(Aorig, T, B, Aperm, TT, BB)
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', N);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', N);
        end
    end

    fprintf('\n');
    fprintf('=================================================\n');
    fprintf('Test permuted strictly upper triangular matrices.\n');
    fprintf('=================================================\n');
    for s = 1:numel(sizes)
        N = sizes(s);
        Tperm = randperm(N);
        Aorig = unbalanced(sparse(triu(ones(N), 1)));
        Aperm = Aorig(Tperm,Tperm);

        [T, B] = balance(full(Aorig));
        [TT, BB] = spbalance(Aperm);

        test_total = test_total + 1;
        if check(Aorig, T, B, Aperm, TT, BB)
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', N);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', N);
        end
    end

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Test permuted lower triangular matrices.\n');
    fprintf('========================================\n');
    for s = 1:numel(sizes)
        N = sizes(s);
        Tperm = randperm(N);
        Aorig = unbalanced(sparse(tril(ones(N))));
        Aperm = Aorig(Tperm,Tperm);
        Aorig = fliplr(flipud(Aorig));

        [T, B] = balance(full(Aorig));
        [TT, BB] = spbalance(Aperm);

        test_total = test_total + 1;
        if check(Aorig, T, B, Aperm, TT, BB)
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', N);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', N);
        end
    end

    fprintf('\n');
    fprintf('=================================================\n');
    fprintf('Test permuted strictly lower triangular matrices.\n');
    fprintf('=================================================\n');
    for s = 1:numel(sizes)
        N = sizes(s);
        Tperm = randperm(N);
        Aorig = unbalanced(sparse(tril(ones(N), -1)));
        Aperm = Aorig(Tperm,Tperm);
        Aorig = fliplr(flipud(Aorig));

        [T, B] = balance(full(Aorig));
        [TT, BB] = spbalance(Aperm);

        test_total = test_total + 1;
        if check(Aorig, T, B, Aperm, TT, BB)
            test_passed = test_passed + 1;
            fprintf('Size %d => passed.\n', N);
        else
            test_failed = test_failed + 1;
            fprintf('Size %d => FAILED!\n', N);
        end
    end

    fprintf('\n');
    fprintf('Summary: total = %d, passed = %d, failed = %d\n',...
        test_total, test_passed, test_failed);
end

function U = unbalanced(A)
    U = A;
    N = size(U, 2);
    for j = 1:N
        I = find(U(:,j));
        if ~isempty(I)
            l = min(U(I,j));
            h = max(U(I,j));
            g = 1000*rand*1000*rand;
            if l == h
                U(I,j) = g;
            else
                U(I,j) = exp((U(I,j)-l)/(h-l)*log(g));
            end
        end
    end
end

function ok = check(Aorig, T, B, Aperm, TT, BB)
    balOk = norm(B - BB, 'fro') <= eps;
    simOk = norm(BB - TT\Aperm*TT, 'fro') <= eps;
    permOk = norm(spones(Aorig) - spones(BB), 'fro') <= eps;
    scalOk = norm(sort(T(find(T))) - sort(TT(find(TT)))) <= eps;

    ok = balOk && simOk && permOk && scalOk;
end
