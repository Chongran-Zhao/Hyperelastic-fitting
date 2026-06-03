function strain = generalized_strain(family, lambda, parameters)
% Return principal values and derivatives for a generalized strain family.
%
% The scale function E(lambda) is evaluated componentwise on the principal
% stretch vector lambda. The derivative is dE/dlambda.
%
% Supported families:
%   SH, Seth-Hill: E = (lambda^m - 1)/m, with Hencky limit at m = 0
%   Hencky:        E = log(lambda)
%   CR:            E = (lambda^m - lambda^(-n))/(m + n)
%   CZ:            E = (2+m)/8*lambda^2 - (2-m)/8*lambda^(-2) - m/4
%   DN:            E = (exp(m*(lambda-1)) - exp(n*(1/lambda-1)))/(m+n)

if nargin < 3
    parameters = [];
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

    case 'hencky'
        require_parameter_count(parameters, 0, 'Hencky');
        strain.values = log(lambda);
        strain.derivatives = 1.0 ./ lambda;
        strain.family = 'Hencky';
        strain.family_label = 'Hencky';
        strain.parameter_names = {};

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

    otherwise
        error('generalized_strain:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end

strain.parameters = parameters;
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
    case {'curnier_rakotomanana', 'curnierrakotomanana'}
        familyId = 'cr';
    case {'curnier_zysset', 'curnierzysset'}
        familyId = 'cz';
    case {'darijani_naghdabadi', 'darijaninaghdabadi'}
        familyId = 'dn';
end
end

function require_parameter_count(parameters, expectedCount, familyName)
if length(parameters) ~= expectedCount
    error('generalized_strain:InvalidParameterCount', ...
        '%s strain expects %d parameter(s), got %d.', ...
        familyName, expectedCount, length(parameters));
end
end
