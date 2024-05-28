% paras = [mu1, alpha1, mu2, alpha2, mu3, alpha3]

%    loading_type    [lambda1,     lambda2,         lambda3    ]
%  Uniaxial tensile  [lambda1, lambda1.^(-0.5), lambda1.^(-0.5)]
% Equbiaxial tensile [lambda1,     lambda1,     lambda1.^(-2.0)]
%      Pure shear    [lambda1,       1.0,       lambda1.^(-1.0)]
classdef Ogden

    properties
        num_paras;
        lambda1_UT;
        lambda2_UT;
        lambda3_UT;
        lambda1_ET;
        lambda2_ET;
        lambda3_ET;
        lambda1_PS;
        lambda2_PS;
        lambda3_PS;
        P1_exp_UT;
        P1_exp_ET;
        P1_exp_PS;
        % paras;
    end

    methods
        function model = Ogden(input_num_paras, input_lambda1_UT, input_lambda1_ET, input_lambda1_PS,...
                input_P1_exp_UT, input_P1_exp_ET, input_P1_exp_PS)
            model.num_paras = input_num_paras;

            model.lambda1_UT = input_lambda1_UT;
            model.lambda2_UT = input_lambda1_UT.^(-0.5);
            model.lambda3_UT = input_lambda1_UT.^(-0.5);

            model.lambda1_ET = input_lambda1_ET;
            model.lambda2_ET = input_lambda1_ET;
            model.lambda3_ET = input_lambda1_ET.^(-2.0);

            model.lambda1_PS = input_lambda1_PS;
            model.lambda2_PS = ones(length(input_lambda1_PS), 1);
            model.lambda3_PS = input_lambda1_PS.^(-1.0);

            model.P1_exp_UT = input_P1_exp_UT;
            model.P1_exp_ET = input_P1_exp_ET;
            model.P1_exp_PS = input_P1_exp_PS;
            % model.paras = ones(num_paras);
        end

        function out = Psi(paras, lambda1, lambda2, lambda3)
            out = (paras(1)./paras(2)) .* (lambda1.^paras(2) + lambda2.^paras(2) + lambda3.^paras(2) - 3)...
                + (paras(3)./paras(4)) .* (lambda1.^paras(4) + lambda2.^paras(4) + lambda3.^paras(4) - 3)...
                + (paras(5)./paras(6)) .* (lambda1.^paras(6) + lambda2.^paras(6) + lambda3.^paras(6) - 3);
        end

        function out = dPsi_dlambda(paras, lambda)
            out = paras(1) .* lambda.^( paras(2)-1.0 )...
                + paras(3) .* lambda.^( paras(4)-1.0 )...
                + paras(5) .* lambda.^( paras(6)-1.0 );
        end

        function out = P1(paras, loading_type)
            switch loading_type
                case 'UT'
                    out = dPsi_dlambda(paras, model.lambda1_UT) - model.lambda1_UT.^(-1.5) .* dPsi_dlambda(paras, model.lambda2_UT);
                case 'ET'
                    out = dPsi_dlambda(paras, model.lambda1_ET) - model.lambda1_ET.^(-3.0) .* dPsi_dlambda(paras, model.lambda3_ET);
                case 'PS'
                    out = dPsi_dlambda(paras, model.lambda1_PS) - model.lambda1_PS.^(-2.0) .* dPsi_dlambda(paras, model.lambda3_PS);
            end
        end

        function out = objective(paras)
            out = sum( ( P1(paras, 'UT') - model.P1_exp_UT).^2 )...
                + sum( ( P1(paras, 'ET') - model.P1_exp_ET).^2 )...
                + sum( ( P1(paras, 'PS') - model.P1_exp_PS).^2 );
        end

    end

end