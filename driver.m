clc; clear; close all;
addpath('src/')
addpath('material_models/')
addpath('tools/')

% Fitting data cases.
fitting_cases = {};

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Treloar', 'UT');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Treloar', 'ET');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Treloar', 'PS');

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawabata', 'UT');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawabata', 'ET');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawabata', 'PS');

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Meunier', 'UT');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Meunier', 'ET');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Meunier', 'PS');

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawamura', 'UT');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawamura', 'UE');
% fitting_cases = Add_exp_data_sets(fitting_cases, 'Kawamura', 'ET');

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Jones', 'UT');

fitting_cases = Add_exp_data_sets(fitting_cases, 'James', 'UT');

% fitting_cases = Add_exp_data_sets(fitting_cases, 'Katashima', 'ET');

models = [];

% models = Add_material_model(models, Neo_Hookean([0.1]), 0.0, Inf);
% models = Add_material_model(models, Mooney_Rivlin([0.1, 0.01]), [0.0, -Inf], [Inf, Inf]);
% models = Add_material_model(models, Yeoh([0.2, 0.01, 0.001]), [0.0, -Inf, -Inf], [Inf, Inf, Inf]);
% models = Add_material_model(models, Ogden([0.1, 1.3]), [0.0, -Inf], [Inf, Inf]);
% models = Add_material_model(models, Arruda_Boyce([0.1, 20.0]), [0.0, 0.0], [Inf, Inf]);

% models = Add_material_model(models, Zhan_Gaussian([1.0]), 0.0, Inf);
models = Add_material_model(models, Zhan_NonGaussian([1.0, 100.0]), [0.0, 0.0], [Inf, Inf]);

% models = Add_material_model(models, Hill_SH([0.01, 1.0]), [0.0, -Inf], [Inf, Inf]);
% models = Add_material_model(models, Hill_Hencky([0.01]), 0.0, Inf);
% models = Add_material_model(models, Hill_CR([0.01, 1.0, 1.0]), [0.0, 0.0, 0.0], [Inf, Inf, Inf]);
% models = Add_material_model(models, Hill_CZ([0.01, 0.0]), [0.0, -2.0], [Inf, 2.0]);
% models = Add_material_model(models, Hill_DN([0.01, 1.0, 1.0]), [0.0, 0.0, 0.0], [Inf, Inf, Inf]);

[models, fit] = start_fit(models, fitting_cases);

plot_simultaneous_fit(models, fitting_cases);

% Prediction data cases.
prediction_cases = {};

% BT data are used for prediction/validation, not simultaneous fitting.
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'UE');
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'ET');
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.7);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.5);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.3);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.1);

% prediction_cases = Add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.0);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.502);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.984);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Jones', 'BT', 2.295);
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Jones', 'BT', 2.623);

prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 1.3);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 1.5);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 1.7);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 2.0);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 2.5);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 3.0);
prediction_cases = Add_exp_data_sets(prediction_cases, 'James', 'BT', 3.5);

% prediction_cases = Add_exp_data_sets(prediction_cases, 'Katashima', 'BT');
% prediction_cases = Add_exp_data_sets(prediction_cases, 'Katashima', 'PS');

if ~isempty(prediction_cases)
    plot_prediction(models, prediction_cases);
end
