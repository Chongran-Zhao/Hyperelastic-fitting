function residual = objective(parameters, models, cases)
% Residual vector for simultaneous fitting of one combined material model.
%
% This is intended for lsqnonlin:
%   objectiveFunction = @(parameters) objective(parameters, models, cases);
%   parameters = lsqnonlin(objectiveFunction, models.parameters, ...
%       models.lower_bounds, models.upper_bounds);

models = models.set_parameters(parameters);
caseList = normalize_cases(cases);

residual = [];
for caseIndex = 1:length(caseList)
    caseData = caseList{caseIndex};
    residual = [residual; case_residual(models, caseData)];
end
end

function caseList = normalize_cases(cases)
if iscell(cases)
    caseList = cases;
elseif isstruct(cases)
    caseList = num2cell(cases);
else
    error('objective:InvalidCases', ...
        'cases must be a cell array or struct array.');
end
end

function residual = case_residual(models, caseData)
residual = [];
for pointIndex = 1:caseData.num_points
    F = caseData.F_list(:, :, pointIndex);
    P_fit = models.P(F);

    if strcmpi(caseData.stress_type, 'scalar')
        [row, col] = component_indices(caseData.stress_component);
        residual(end + 1, 1) = P_fit(row, col) - caseData.P_exp(pointIndex);
    elseif strcmpi(caseData.stress_type, 'relative_cauchy')
        sigmaValue = relative_cauchy_stress(P_fit, F, ...
            caseData.relative_cauchy_component);
        residual(end + 1, 1) = sigmaValue - caseData.sigma_exp(pointIndex);
    else
        P_exp = caseData.P_list(:, :, pointIndex);
        mask = caseData.P_mask(:, :, pointIndex);
        difference = P_fit - P_exp;
        residual = [residual; difference(mask)];
    end
end
end

function value = relative_cauchy_stress(P, F, name)
J = det(F);
sigma = P * F' ./ J;

switch upper(char(name))
    case 'SIGMA11_MINUS_SIGMA22'
        value = sigma(1, 1) - sigma(2, 2);
    otherwise
        error('objective:UnsupportedRelativeCauchyComponent', ...
            'Unsupported relative Cauchy component: %s.', name);
end
end

function [row, col] = component_indices(name)
switch upper(char(name))
    case 'P11'
        row = 1; col = 1;
    case 'P22'
        row = 2; col = 2;
    case 'P33'
        row = 3; col = 3;
    case 'P12'
        row = 1; col = 2;
    case 'P21'
        row = 2; col = 1;
    otherwise
        error('objective:UnsupportedComponent', ...
            'Unsupported PK1 component: %s.', name);
end
end
