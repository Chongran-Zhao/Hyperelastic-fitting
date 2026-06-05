function metrics = evaluation_metrics(seriesList)
%   Averaged NMAD and R^2 over plotted data series.
%
%   metrics = EVALUATION_METRICS(seriesList) computes two fitting-quality
%   indicators for each plotted experimental/predicted data series and then
%   reports their arithmetic mean over all valid series.
%
%   The normalized mean absolute difference is computed as
%
%       NMAD = 100 * mean(abs(p - e)) / max(mean(abs(e)), mean(abs(p))),
%
%   where e denotes the experimental values and p denotes the predicted
%   values. This implementation follows the max-normalized NMAD definition
%   used in MCalibration-style polymer/rubber material calibration studies,
%   including:
%
%       - Gonzalez-Vega et al. (2022), Experimental viscoelastic properties
%         evaluation, under impact loads and large strain conditions, of
%         coated & uncoated rubber from end-of-life tires, Polymer Testing
%         107, 107468.
%
%       - Gonzalez-Vega et al. (2024), Impact properties of an end of life
%         tires' rubber. Numerical validation considering large strain and
%         strain rate conditions.
%
%   The general use of NMAD as a material-parameter extraction objective
%   function is also discussed in:
%
%       - Bergstrom (2015), Mechanics of Solid Polymers: Theory and
%         Computational Modeling, Chapter 9, Eq. (9.5).
%
%   Note that Bergstrom's Eq. (9.5) normalizes by mean(abs(e)) only, whereas
%   the present implementation uses the symmetric max-normalized form.
%
%   The coefficient of determination is computed as
%
%       R^2 = 1 - SS_res / SS_tot,
%       SS_res = sum((e - p).^2),
%       SS_tot = sum((e - mean(e)).^2).
%
%   This follows the standard residual-sum-of-squares definition; see also
%   Bergstrom (2015), Chapter 9, Eq. (9.6).

metrics = struct();
metrics.NMAD = NaN;
metrics.R_square = NaN;
metrics.num_series = 0;

if isempty(seriesList)
    return;
end

nmadValues = [];
rSquareValues = [];

for seriesIndex = 1:length(seriesList)
    expValues = seriesList(seriesIndex).exp(:);
    fitValues = seriesList(seriesIndex).fit(:);

    validMask = isfinite(expValues) & isfinite(fitValues);
    expValues = expValues(validMask);
    fitValues = fitValues(validMask);

    if isempty(expValues)
        continue;
    end

    nmadValues(end + 1) = evaluation_nmad(fitValues, expValues);
    rSquareValues(end + 1) = evaluation_r_square(expValues, fitValues);
end

if isempty(nmadValues)
    return;
end

metrics.NMAD = mean(nmadValues, 'omitnan');

finiteRSquare = rSquareValues(isfinite(rSquareValues));
if ~isempty(finiteRSquare)
    metrics.R_square = mean(finiteRSquare);
end

metrics.num_series = sum(isfinite(nmadValues));
end

function out = evaluation_nmad(preValues, expValues)
meanAbsError = mean(abs(preValues - expValues));
meanAbsExp = mean(abs(expValues));
meanAbsPre = mean(abs(preValues));

denominator = max(meanAbsExp, meanAbsPre);

if denominator <= eps
    if meanAbsError <= eps
        out = 0.0;
    else
        out = NaN;
    end
else
    out = 100.0 * meanAbsError / denominator;
end
end

function out = evaluation_r_square(expValues, preValues)
residualSum = sum((expValues - preValues).^2);
totalSum = sum((expValues - mean(expValues)).^2);

if totalSum <= eps
    if residualSum <= eps
        out = 1.0;
    else
        out = NaN;
    end
else
    out = 1.0 - residualSum / totalSum;
end
end