% This function is used to calculate the 1st principal value
% of 1st PK stress
% lambda is a scalar double value here
% tag in {UT, ET, PS}
CR_PK_handle = @(paras, lambda, tag) ( ...
    result = 0;
    switch tag
        case 'UT'
            for ii = 1:size(paras, 1)
                result = result + (dCR_dLambda(paras(ii,:), lambda) - lambda.^(-1.5) .* dCR_dLambda(paras(ii,:), lambda.^(-0.5)));
            end
        case 'ET'
            for ii = 1:size(paras, 1)
                result = result + (dCR_dLambda(paras(ii,:), lambda) - lambda.^(-3.0) .* dCR_dLambda(paras(ii,:), lambda.^(-2.0)));
            end
        case 'PS'
            for ii = 1:size(paras, 1)
                result = result + (dCR_dLambda(paras(ii,:), lambda) - lambda.^(-2.0) .* dCR_dLambda(paras(ii,:), lambda.^(-1.0)));
            end
        otherwise
            error('Invalid tag value. Expected "UT", "ET", or "PS".');
    end
    result
);