clc; clear; close all;
addpath('src/')
addpath('material_models/')
addpath('tools/')

% Fitting data cases.
fitting_cases = {};
% Prediction data cases.
prediction_cases = {};

% fitting_cases = add_exp_data_sets(fitting_cases, 'Treloar', 'UT');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Treloar', 'ET');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Treloar', 'PS');

% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawabata', 'UT');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawabata', 'ET');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawabata', 'PS');

% fitting_cases = add_exp_data_sets(fitting_cases, 'Meunier', 'UT');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Meunier', 'ET');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Meunier', 'PS');

% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawamura', 'UT');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawamura', 'UE');
% fitting_cases = add_exp_data_sets(fitting_cases, 'Kawamura', 'ET');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'UE');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'ET');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.7);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.5);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.3);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Kawamura', 'BT', 1.1);

% fitting_cases = add_exp_data_sets(fitting_cases, 'Jones', 'UT');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.0);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.502);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Jones', 'BT', 1.984);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Jones', 'BT', 2.295);
% prediction_cases = add_exp_data_sets(prediction_cases, 'Jones', 'BT', 2.623);

fitting_cases = add_exp_data_sets(fitting_cases, 'James', 'UT');
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 1.3);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 1.5);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 1.7);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 2.0);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 2.5);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 3.0);
prediction_cases = add_exp_data_sets(prediction_cases, 'James', 'BT', 3.5);

% fitting_cases = add_exp_data_sets(fitting_cases, 'Katashima', 'ET');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Katashima', 'BT');
% prediction_cases = add_exp_data_sets(prediction_cases, 'Katashima', 'PS');

models = [];

% models = add_material_model(models, Neo_Hookean([0.1]), 0.0, Inf);
% models = add_material_model(models, Mooney_Rivlin([0.1, 0.01]), [0.0, -Inf], [Inf, Inf]);
% models = add_material_model(models, Yeoh([0.2, 0.01, 0.001]), [0.0, -Inf, -Inf], [Inf, Inf, Inf]);
% models = add_material_model(models, Ogden([0.1, 1.3]), [0.0, -Inf], [Inf, Inf]);
% models = add_material_model(models, Arruda_Boyce([0.1, 20.0]), [0.0, 0.0], [Inf, Inf]);
% models = add_material_model(models, Zhan_Gaussian([1.0]), 0.0, Inf);
% models = add_material_model(models, Zhan_NonGaussian([1.0, 100.0]), [0.0, 0.0], [Inf, Inf]);
models = add_material_model(models, Micro_GenStrain('Gaussian', 'Biot', 'Biot'));
% models = add_material_model(models, Hill_GenStrain('CR'));

[models, fit] = start_fit(models, fitting_cases);

plot_simultaneous_fit(models, fitting_cases);

if ~isempty(prediction_cases)
    plot_prediction(models, prediction_cases);
end
