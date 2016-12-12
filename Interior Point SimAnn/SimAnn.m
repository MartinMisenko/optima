function [x, itr_x, itr_y] = SimAnn(fun, x0, kmax, T)
    % SimAnn searchs a given function using Simulated Annealing method to
    % help to find its global optimum. Performance of the searching depends on number
    % of iterations, temperature and starting point. 
    %
    % Parameters kmax - maximum number of iteration is optional as well as
    % parameter T - temperature. Default number of iterations is 2000 and
    % default value of temperature is 150.
    %
    %  x = SimAnn(fun, x0, kmax, T)
    %
    % [x, itr_x, itr_y] = SimAnn(fun, x0, kmax, T)
    %       The function also returns values of all itrations


    if ~exist('T','var')
          T0 = 150;
    end
    
    if ~exist('kmax','var')
          kmax = 2000;
    end

    %% Iterations initializing
    itr_y = 0;
    i = 1; % Number of acceptances

    T0 = 150;
    x = x0;
    e = fun(x);

    for k=1:kmax
        %% Appropriate decrease of the temperature
        T = T0*0.95^k;
        if T < 1e-4
            break
        end

        %% Choosing a new neighbour x
        ex = rand/log(k);   % Gain (Magic)
        if rand > 0.5       % Random decision about sign of the gain
            ex = ex*(-1);
        end
        x_n = ex + x;

        e_n = fun(x_n);  % Compute a new function value

        %% Assessment of the new function value 
        if e >= e_n
            acp = 1;
        else
            switch e_n
                case Inf
                    acp = 0;
                case -Inf
                    disp('Function has no minimum');
                    break
                otherwise
                    acp = exp(-(e_n - e)/T);
            end
        end

        if rand<acp && min(abs(x) ~= Inf)
            x = x_n;
            for j=1:length(x)
                itr_x{i,j} = x(j);   % Log x iteration
            end
            e = e_n;
            itr_y(i) = e;   % Log y iteration

            i = i + 1;      % Increase number of acceptances 
        end
    end
    itr_x = cell2mat(itr_x);
end