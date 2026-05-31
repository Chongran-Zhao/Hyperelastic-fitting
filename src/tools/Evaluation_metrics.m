function metrics = Evaluation_metrics(seriesList)
%EVALUATION_METRICS Average NMAD and R-square over plotted data series.

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

metrics.NMAD = mean(nmadValues);
finiteRSquare = rSquareValues(isfinite(rSquareValues));
if ~isempty(finiteRSquare)
    metrics.R_square = mean(finiteRSquare);
end
metrics.num_series = length(nmadValues);
end

function out = evaluation_nmad(preValues, expValues)
meanAbsError = mean(abs(preValues - expValues));
meanAbsExp = mean(abs(expValues));
meanAbsPre = mean(abs(preValues));
out = 100.0 * meanAbsError / max(meanAbsExp, meanAbsPre);
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
