function Y = FUNC_LegendreCoeff2Data(X,LC)
    N = length(LC)-1; % degree of polynomial
    m = length(X); % length of data
    if N > 1
        coeff = zeros(N+1);
        coeff([1 N+3]) = 1;
        for ii = 3:N+1
            coeff(ii,:) = (2-1/(ii-1))*coeff(ii-1,[end 1:end-1]) - (1-1/(ii-1))*coeff(ii-2,:);
        end
    else
        coeff = eye(N+1);
    end
    M = -1:2/(m-1):1;
    M = M(:);
    D = cumprod([ones(m,1) M(:,ones(1,N))], 2) * coeff.';
    Y = D*LC;     
end