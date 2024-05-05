% This SH-energy is based on the quadratic form of SH strain
% paras = {mu, mm}
% lambda = {lambda_1, lambda_2, lambda_3}
function out = SH(paras, lambda)

mu = paras(1);
mm = paras(2);
out = mu ./ (mm .* mm) .* (...
        (lambda(1).^mm - 1.0) .* (lambda(1).^mm - 1.0)...
    +   (lambda(2).^mm - 1.0) .* (lambda(2).^mm - 1.0)...
    +   (lambda(3).^mm - 1.0) .* (lambda(3).^mm - 1.0));
end



% UT_x = linspace(1.0, max(Treloar_UT_strain), 25);
% ET_x = linspace(1.0, max(Treloar_ET_strain), 25);
% PS_x = linspace(1.0, max(Treloar_PS_strain), 25);
% 
% figure;
% hold on;
% 
% plot(Treloar_UT_strain, Treloar_UT_stress, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
% plot(UT_x, CR_PK(paras, UT_x, 'UT'), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2, 'LineStyle', '-');
% 
% plot(Treloar_ET_strain, Treloar_ET_stress, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); 
% plot(ET_x, CR_PK(paras, ET_x, 'ET'), 'Color', [0.85, 0.33, 0], 'LineWidth', 2); % 橙色
% 
% plot(Treloar_PS_strain, Treloar_PS_stress, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8); 
% plot(PS_x, CR_PK(paras, PS_x, 'PS'), 'Color', [0, 0.5, 0.5], 'LineWidth', 2); % 青色
% 
% hold off;
% grid off;