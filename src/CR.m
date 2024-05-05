% This CR-energy is based on the quadratic form of CR strain
% paras = {mu, mm, nn}
% lambda = {lambda_1, lambda_2, lambda_3}
% function out = CR(paras, lambda)
% 
% mu = paras(1);
% mm = paras(2);
% nn = paras(3);
% out = mu ./ ( (mm+nn) .* (mm+nn) ) .* (...
%         (lambda(1).^mm - lambda(1).^(-nn)) .* (lambda(1).^mm - lambda(1).^(-nn))...
%     +   (lambda(2).^mm - lambda(2).^(-nn)) .* (lambda(2).^mm - lambda(2).^(-nn))...
%     +   (lambda(3).^mm - lambda(3).^(-nn)) .* (lambda(3).^mm - lambda(3).^(-nn)) );
% end