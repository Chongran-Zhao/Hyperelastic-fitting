function out = Hill_Q_proj(lambda, strainValues, strainDerivatives, V)
out = zeros(3, 3, 3, 3);

d = strainDerivatives ./ lambda;
theta = zeros(3, 3);
for ii = 1:3
    for jj = 1:3
        samePrincipalValue = abs(lambda(ii) - lambda(jj)) <= ...
            1.0e-12 .* max([1.0, abs(lambda(ii)), abs(lambda(jj))]);
        if ii == jj || samePrincipalValue
            theta(ii, jj) = d(ii);
        else
            theta(ii, jj) = 2.0 .* (strainValues(ii) - strainValues(jj)) ./ ...
                (lambda(ii) .* lambda(ii) - lambda(jj) .* lambda(jj));
        end
    end
end

for ii = 1:3
    Ni = V(:, ii);
    out = out + d(ii) .* cross_otimes_1d_to_4d(Ni, Ni, Ni, Ni);
end

for ii = 1:3
    for jj = (ii + 1):3
        Ni = V(:, ii);
        Nj = V(:, jj);
        out = out + theta(ii, jj) .* dot_otimes(Ni, Nj) + ...
            theta(ii, jj) .* dot_otimes(Nj, Ni);
    end
end
end

function out = dot_otimes(Na, Nb)
out = 0.5 .* cross_otimes_1d_to_4d(Na, Nb, Na, Nb) + ...
    0.5 .* cross_otimes_1d_to_4d(Na, Nb, Nb, Na);
end

function out = cross_otimes_1d_to_4d(Na, Nb, Nc, Nd)
out = zeros(3, 3, 3, 3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii, jj, kk, ll) = Na(ii) .* Nb(jj) .* Nc(kk) .* Nd(ll);
            end
        end
    end
end
end
