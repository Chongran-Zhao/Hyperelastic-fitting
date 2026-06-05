function strain = generalized_strain(family, lambda, parameters, V)
% Return principal values, derivatives, and inverse for a strain family.
%
% The scale function E(lambda) is evaluated componentwise on the principal
% stretch vector lambda. The derivative is dE/dlambda.
% The inverse handle strain.inverse(xi) returns lambda = E^{-1}(xi).
% If the principal direction matrix V is provided, the output also contains
% Hill's fourth-order projection tensor Q = 2*dE/dC.
%
% Supported families:
%   SH, Seth-Hill: E = (lambda^m - 1)/m, with Hencky limit at m = 0
%   Hencky:        E = log(lambda)
%   Biot:          E = lambda - 1
%   BI:            E = (lambda^m - lambda^(-m))/(2*m)
%   CR:            E = (lambda^m - lambda^(-n))/(m + n)
%   CZ:            E = (2+m)/8*lambda^2 - (2-m)/8*lambda^(-2) - m/4
%   DN:            E = (exp(m*(lambda-1)) - exp(n*(1/lambda-1)))/(m+n)

if nargin < 3
    parameters = [];
end
if nargin < 4
    V = [];
end

lambda = lambda(:);
parameters = parameters(:).';
familyId = normalize_family(family);

switch familyId
    case 'sh'
        require_parameter_count(parameters, 1, 'Seth-Hill');
        m = parameters(1);
        if abs(m) < 1.0e-12
            strain.values = log(lambda);
            strain.derivatives = 1.0 ./ lambda;
        else
            strain.values = (lambda .^ m - 1.0) ./ m;
            strain.derivatives = lambda .^ (m - 1.0);
        end
        strain.family = 'SH';
        strain.family_label = 'Seth-Hill';
        strain.parameter_names = {'m'};
        strain.inverse = @(xi) inverse_seth_hill(xi, m);

    case 'hencky'
        require_parameter_count(parameters, 0, 'Hencky');
        strain.values = log(lambda);
        strain.derivatives = 1.0 ./ lambda;
        strain.family = 'Hencky';
        strain.family_label = 'Hencky';
        strain.parameter_names = {};
        strain.inverse = @(xi) exp(xi);

    case 'biot'
        require_parameter_count(parameters, 0, 'Biot');
        strain.values = lambda - 1.0;
        strain.derivatives = ones(size(lambda));
        strain.family = 'Biot';
        strain.family_label = 'Biot';
        strain.parameter_names = {};
        strain.inverse = @(xi) 1.0 + xi;

    case 'bi'
        require_parameter_count(parameters, 1, 'Bazant-Itskov');
        m = parameters(1);
        if abs(m) < 1.0e-12
            strain.values = log(lambda);
            strain.derivatives = 1.0 ./ lambda;
        else
            strain.values = (lambda .^ m - lambda .^ (-m)) ./ (2.0 .* m);
            strain.derivatives = 0.5 .* ...
                (lambda .^ (m - 1.0) + lambda .^ (-m - 1.0));
        end
        strain.family = 'BI';
        strain.family_label = 'Bazant-Itskov';
        strain.parameter_names = {'m'};
        strain.inverse = @(xi) inverse_bazant_itskov(xi, m);

    case 'cr'
        require_parameter_count(parameters, 2, 'Curnier-Rakotomanana');
        m = parameters(1);
        n = parameters(2);
        if m .* n <= 0.0
            error('generalized_strain:InvalidParameter', ...
                'Curnier-Rakotomanana strain requires m*n > 0.');
        end
        strain.values = (lambda .^ m - lambda .^ (-n)) ./ (m + n);
        strain.derivatives = m ./ (m + n) .* lambda .^ (m - 1.0) + ...
            n ./ (m + n) .* lambda .^ (-n - 1.0);
        strain.family = 'CR';
        strain.family_label = 'Curnier-Rakotomanana';
        strain.parameter_names = {'m', 'n'};
        strain.inverse = @(xi) inverse_curnier_rakotomanana(xi, m, n);

    case 'cz'
        require_parameter_count(parameters, 1, 'Curnier-Zysset');
        m = parameters(1);
        if m < -2.0 || m > 2.0
            error('generalized_strain:InvalidParameter', ...
                'Curnier-Zysset strain requires -2 <= m <= 2.');
        end
        strain.values = (2.0 + m) ./ 8.0 .* lambda .^ 2.0 - ...
            (2.0 - m) ./ 8.0 .* lambda .^ (-2.0) - m ./ 4.0;
        strain.derivatives = (2.0 + m) ./ 4.0 .* lambda + ...
            (2.0 - m) ./ 4.0 .* lambda .^ (-3.0);
        strain.family = 'CZ';
        strain.family_label = 'Curnier-Zysset';
        strain.parameter_names = {'m'};
        strain.inverse = @(xi) inverse_curnier_zysset(xi, m);

    case 'dn'
        require_parameter_count(parameters, 2, 'Darijani-Naghdabadi');
        m = parameters(1);
        n = parameters(2);
        if m <= 0.0 || n <= 0.0
            error('generalized_strain:InvalidParameter', ...
                'Darijani-Naghdabadi strain requires m > 0 and n > 0.');
        end
        strain.values = (exp(m .* (lambda - 1.0)) - ...
            exp(n .* (lambda .^ (-1.0) - 1.0))) ./ (m + n);
        strain.derivatives = (m .* exp(m .* (lambda - 1.0)) + ...
            n .* lambda .^ (-2.0) .* ...
            exp(n .* (lambda .^ (-1.0) - 1.0))) ./ (m + n);
        strain.family = 'DN';
        strain.family_label = 'Darijani-Naghdabadi';
        strain.parameter_names = {'m', 'n'};
        strain.inverse = @(xi) inverse_darijani_naghdabadi(xi, m, n);

    otherwise
        error('generalized_strain:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end

strain.parameters = parameters;
if ~isempty(V)
    strain.Q = hill_Q_proj(lambda, strain.values, strain.derivatives, V);
end
end

function familyId = normalize_family(family)
familyId = lower(char(family));
familyId = strrep(familyId, '-', '_');
familyId = strrep(familyId, ' ', '_');

switch familyId
    case {'seth_hill', 'sethhill'}
        familyId = 'sh';
    case {'he', 'log', 'logarithmic'}
        familyId = 'hencky';
    case {'bazant_itskov', 'bazantitskov'}
        familyId = 'bi';
    case {'curnier_rakotomanana', 'curnierrakotomanana'}
        familyId = 'cr';
    case {'curnier_zysset', 'curnierzysset'}
        familyId = 'cz';
    case {'darijani_naghdabadi', 'darijaninaghdabadi'}
        familyId = 'dn';
end
end

function lambda = inverse_seth_hill(xi, m)
if abs(m) < 1.0e-12
    lambda = exp(xi);
else
    lambda = (1.0 + m .* xi) .^ (1.0 ./ m);
end
end

function lambda = inverse_bazant_itskov(xi, m)
if abs(m) < 1.0e-12
    lambda = exp(xi);
else
    lambda = exp(asinh(m .* xi) ./ m);
end
end

function lambda = inverse_curnier_rakotomanana(xi, m, n)
strainFunction = @(lambda) (lambda .^ m - lambda .^ (-n)) ./ (m + n);
lambda = inverse_monotone_strain(xi, strainFunction);
end

function lambda = inverse_curnier_zysset(xi, m)
if abs(m + 2.0) < 1.0e-12
    lambda = (1.0 - 2.0 .* xi) .^ (-0.5);
elseif abs(m - 2.0) < 1.0e-12
    lambda = sqrt(1.0 + 2.0 .* xi);
else
    numerator = 2.0 .* m + 8.0 .* xi + ...
        sqrt((2.0 .* m + 8.0 .* xi) .^ 2.0 + ...
        4.0 .* (2.0 + m) .* (2.0 - m));
    denominator = 2.0 .* (2.0 + m);
    lambda = sqrt(numerator ./ denominator);
end
end

function lambda = inverse_darijani_naghdabadi(xi, m, n)
strainFunction = @(lambda) (exp(m .* (lambda - 1.0)) - ...
    exp(n .* (lambda .^ (-1.0) - 1.0))) ./ (m + n);
lambda = inverse_monotone_strain(xi, strainFunction);
end

function lambda = inverse_monotone_strain(xi, strainFunction)
lambda = arrayfun(@(target) inverse_monotone_scalar(target, strainFunction), xi);
end

function lambda = inverse_monotone_scalar(target, strainFunction)
if ~isfinite(target)
    lambda = NaN;
    return;
end

if target == 0.0
    lambda = 1.0;
    return;
end

if target > 0.0
    lower = 1.0;
    upper = 2.0;
    while strainFunction(upper) < target
        upper = 2.0 .* upper;
    end
else
    lower = 0.5;
    upper = 1.0;
    while strainFunction(lower) > target
        lower = 0.5 .* lower;
    end
end

for iteration = 1:100
    midpoint = exp(0.5 .* (log(lower) + log(upper)));
    value = strainFunction(midpoint);

    if value < target
        lower = midpoint;
    else
        upper = midpoint;
    end
end

lambda = exp(0.5 .* (log(lower) + log(upper)));
end

function out = hill_Q_proj(lambda, strainValues, strainDerivatives, V)
% Hill fourth-order projection tensor Q = 2*dE/dC.
%
% E = sum_a E_a * M_a, C = sum_a lambda_a^2 * M_a,
% M_a = N_a * N_a', and V = [N_1, N_2, N_3].
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

function require_parameter_count(parameters, expectedCount, familyName)
if length(parameters) ~= expectedCount
    error('generalized_strain:InvalidParameterCount', ...
        '%s strain expects %d parameter(s), got %d.', ...
        familyName, expectedCount, length(parameters));
end
end
