clc;clear;
addpath('src/');
% loading data
Treloar_UT_strain = importdata("./Treloar-UT/strain.txt");
Treloar_UT_stress = importdata("./Treloar-UT/stress.txt");

Treloar_ET_strain = importdata("./Treloar-ET/strain.txt");
Treloar_ET_stress = importdata("./Treloar-ET/stress.txt");

Treloar_PS_strain = importdata("./Treloar-PS/strain.txt");
Treloar_PS_stress = importdata("./Treloar-PS/stress.txt");

paras0 = [1.0, 1.0, 2.0; 1.0, 1.0, 2.0];
objectiveFunction = @(paras) objective(paras, @CR_PK, Treloar_UT_strain, Treloar_UT_stress, Treloar_ET_strain, Treloar_ET_stress, Treloar_PS_strain, Treloar_PS_stress);
options = optimoptions('lsqnonlin', 'Algorithm', 'interior-point', 'MaxIterations', 5000);
[paras, res] = lsqnonlin( objectiveFunction, paras0, [], [], options);

function out = objective(paras, CR_PK, lambda_UT, PK_UT, lambda_ET, PK_ET, lambda_PS, PK_PS)
    out = sum( (CR_PK(paras, @dCR_dLambda, lambda_UT, 'UT') - PK_UT).^2 ) / length(lambda_UT) ...
        + sum( (CR_PK(paras, @dCR_dLambda, lambda_ET, 'ET') - PK_ET).^2 ) / length(lambda_ET) ...
        + sum( (CR_PK(paras, @dCR_dLambda, lambda_PS, 'PS') - PK_PS).^2 ) / length(lambda_PS);
end