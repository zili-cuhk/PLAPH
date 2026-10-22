function result = run_PLAPH_real_data(analysis_label, varargin)
%RUN_PLAPH_REAL_DATA Fit the proposed copula-PLAPH-MIC method to TCGA-LUAD.
%
% This wrapper uses the same MPL/MIC/block-CD implementation as the
% simulation code. It only supplies real-data matrices, tuning options, and
% reproducible output files.

if nargin < 1 || strlength(string(analysis_label)) == 0
    analysis_label = "primary_p500_q200";
end

opts = parse_real_data_options(varargin{:});
analysis_label = string(analysis_label);
root = string(opts.WorkingSetDir);
[source_p, source_q] = source_dimensions_for_label(analysis_label, opts);

scenario = struct();
scenario.real_data = true;
scenario.scenario_id = "TCGA_LUAD_" + analysis_label;
scenario.fit_copula = string(opts.FitCopula);
scenario.tau = opts.Tau;
scenario.profile_tau = opts.EstimateTau;
scenario.p = opts.P;
scenario.q = opts.Q;
scenario.clinical_file = fullfile(root, analysis_label + "_clinical.csv");
scenario.z_file = fullfile(root, analysis_label + "_Z_RNA_p" + string(source_p) + ".csv");
scenario.x_file = fullfile(root, analysis_label + "_X_CNA_q" + string(source_q) + ".csv");
scenario.z_feature_file = fullfile(root, analysis_label + "_Z_RNA_p" + string(source_p) + "_features.csv");
scenario.x_feature_file = fullfile(root, analysis_label + "_X_CNA_q" + string(source_q) + "_features.csv");

if ~exist(opts.OutputDir, "dir")
    mkdir(opts.OutputDir);
end

result = plaph_run_one_scenario_impl(scenario, opts);
result.fit_copula = string(opts.FitCopula);
result.estimate_tau = opts.EstimateTau;
result.mean_loglik = result.diagnostics.mpl_final_mean_loglik;
result.total_loglik = result.diagnostics.mpl_final_total_loglik;
result.spline_group_df = max(result.num_basis-1,0);
% v2: rank-based degrees of freedom (selected joint design rank in each
% margin, both baseline step functions, and tau when profiled).
result.effective_df = result.df_rank + double(result.estimate_tau);
result.effective_df_nominal = approximate_real_data_df(result);
result.aic = -2 * result.total_loglik + 2 * result.effective_df;
result.bic = -2 * result.total_loglik + log(result.n) * result.effective_df;
result.options = opts;
result.estimator_definition = struct( ...
    'linear_penalty',string(opts.LinearMICForm), ...
    'linear_kink_coefficient',opts.LinearKinkCoefficient, ...
    'mic_penalty_scale',opts.MICPenaltyScale, ...
    'mic_penalty_normalization',string(opts.MICPenaltyNormalization), ...
    'additive_penalty',string(opts.SplineGroupPenaltyType), ...
    'group_iht_threshold_failure',opts.SplineGroupIHTThresholdFailure, ...
    'group_iht_threshold_censoring',opts.SplineGroupIHTThresholdCensoring, ...
    'selection_support',"frozen_before_postselection_refit", ...
    'post_selection_refit',string(opts.PostSelectionRefitType), ...
    'fmincon_polish',opts.MPLOracleRefitFminconPolish, ...
    'baseline_target_bins',opts.BaselineBinCount);
save(fullfile(opts.OutputDir, analysis_label + "_PLAPH_result.mat"), "result", "-v7.3");

writetable(result.linear_failure, fullfile(opts.OutputDir, analysis_label + "_linear_failure.csv"));
writetable(result.linear_censoring, fullfile(opts.OutputDir, analysis_label + "_linear_censoring.csv"));
writetable(result.additive_failure, fullfile(opts.OutputDir, analysis_label + "_additive_failure.csv"));
writetable(result.additive_censoring, fullfile(opts.OutputDir, analysis_label + "_additive_censoring.csv"));
writetable(result.function_curves_failure, fullfile(opts.OutputDir, analysis_label + "_function_curves_failure.csv"));
writetable(result.function_curves_censoring, fullfile(opts.OutputDir, analysis_label + "_function_curves_censoring.csv"));

summary = table( ...
    string(result.scenario.scenario_id), string(result.fit_copula), result.estimate_tau, ...
    result.n, result.p, result.q, result.num_basis, result.spline_group_df, ...
    result.n_failure_events, result.n_censoring_events, result.censor_rate, ...
    result.tau_hat, result.alpha_hat, result.mean_loglik, result.total_loglik, ...
    result.effective_df, result.aic, result.bic, result.elapsed_seconds, ...
    sum(result.linear_failure.selected), sum(result.linear_censoring.selected), ...
    sum(result.additive_failure.selected), sum(result.additive_censoring.selected), ...
    'VariableNames', {'scenario_id', 'fit_copula', 'estimate_tau', ...
    'n', 'p', 'q', 'num_basis', 'spline_group_df', ...
    'n_failure_events', 'n_censoring_events', 'censor_rate', ...
    'tau_hat', 'alpha_hat', 'mean_loglik', 'total_loglik', ...
    'effective_df', 'aic', 'bic', 'elapsed_seconds', ...
    'selected_linear_failure', 'selected_linear_censoring', ...
    'selected_additive_failure', 'selected_additive_censoring'});
writetable(summary, fullfile(opts.OutputDir, analysis_label + "_summary.csv"));

disp(summary);
end


function df = approximate_real_data_df(result)
% Approximate effective degrees of freedom for copula-family comparison.
% The count includes selected linear coefficients, selected spline
% coefficients by group, both baseline step functions, and the copula
% association parameter when it is estimated.
selected_linear = sum(result.linear_failure.selected) + sum(result.linear_censoring.selected);
selected_additive = sum(result.additive_failure.selected) + sum(result.additive_censoring.selected);
baseline_df = numel(result.seta) + numel(result.gama);
copula_df = double(result.estimate_tau);
df = selected_linear + result.spline_group_df * selected_additive + ...
    baseline_df + copula_df;
end


function opts = parse_real_data_options(varargin)
parser = inputParser;
default_working_dir = "D:\科研项目\课题9\real_data\plaph_luad_noimage_pilot_screening\working_sets";
default_output_dir = "D:\科研项目\课题9\real_data\plaph_luad_noimage_plaph_results";
parser.addParameter("WorkingSetDir", default_working_dir, @(x) ischar(x) || isstring(x));
parser.addParameter("OutputDir", default_output_dir, @(x) ischar(x) || isstring(x));
parser.addParameter("P", 500, @(x) isscalar(x) && x > 0);
parser.addParameter("Q", 200, @(x) isscalar(x) && x > 0);
parser.addParameter("SourceP", [], @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("SourceQ", [], @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("FitCopula", "Frank", @(x) any(strcmpi(string(x), ["Frank" "Gaussian" "Clayton" "Gumbel"])));
parser.addParameter("Tau", 0.2, ...
    @(x) isscalar(x) && isfinite(x) && abs(x) < 1);
parser.addParameter("EstimateTau", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("BaseSeed", 20260707, @(x) isscalar(x) && x > 0);
parser.addParameter("MaxOuterIter", 50, @(x) isscalar(x) && x > 0);
parser.addParameter("NumInternalKnots", 8, @(x) isscalar(x) && x > 0);
parser.addParameter("ThresholdLinear", 5e-4, @(x) isscalar(x) && x >= 0);
parser.addParameter("ThresholdLinearFailure", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("ThresholdLinearCensoring", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("ThresholdSpline", 1e-3, @(x) isscalar(x) && x >= 0);
parser.addParameter("MICPenaltyScale", 10, @(x) isscalar(x) && x > 0);
parser.addParameter("MICPenaltyScaleFailure", [], ...
    @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("MICPenaltyScaleCensoring", [], ...
    @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("MICPenaltyNormalization", "sum", @(x) any(strcmpi(string(x), ["mean" "sum"])));
parser.addParameter("MICSharpness", 90, @(x) isscalar(x) && x > 0);
parser.addParameter("LinearMICForm", "smooth", @(x) any(strcmpi(string(x), ...
    ["smooth" "kinked" "scad" "mcp"])));
parser.addParameter("LinearKinkCoefficient", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("LinearL1Type", "none", @(x) any(strcmpi(string(x), ["none" "lasso" "adaptive"])));
parser.addParameter("LinearL1Scale", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("LinearRidgeScale", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("AdaptiveWeightGamma", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("AdaptiveWeightEpsilon", 0.05, @(x) isscalar(x) && x > 0);
parser.addParameter("AdaptiveWeightCap", 100, @(x) isscalar(x) && x >= 1);
parser.addParameter("InitialLinearStrategy", "legacy", @(x) any(strcmpi(string(x), ["legacy" "zero" "continuation"])));
parser.addParameter("InitialLinearPilotScale", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("InitialLinearCap", Inf, @(x) isscalar(x) && x > 0);
parser.addParameter("InitialSplineWeighting", "legacy_tail", ...
    @(x) any(strcmpi(string(x), ["none" "legacy_tail"])));
parser.addParameter("InitialSplineStrategy", "legacy", ...
    @(x) any(strcmpi(string(x), ["legacy" "group_lasso"])));
parser.addParameter("InitialSplinePilotScale", 1, ...
    @(x) isscalar(x) && x > 0);
parser.addParameter("MPLStepScale", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("MPLDecomposablePilotIterations", 0, ...
    @(x) isscalar(x) && x >= 0 && fix(x) == x);
parser.addParameter("MPLDecomposablePilotLinearScale", 1, ...
    @(x) isscalar(x) && x > 0);
parser.addParameter("MPLDecomposablePilotGroupScale", 1, ...
    @(x) isscalar(x) && x > 0);
parser.addParameter("MPLDecomposablePilotStart", "zero", ...
    @(x) any(strcmpi(string(x), ["zero" "initial"])));
parser.addParameter("MPLDecomposablePilotWeighting", "none", ...
    @(x) any(strcmpi(string(x), ["none" "adaptive"])));
parser.addParameter("MPLDecomposablePilotGroupWeightGamma", 1, ...
    @(x) isscalar(x) && x > 0);
parser.addParameter("MPLStepStrategy", "fixed", @(x) any(strcmpi(string(x), ...
    ["fixed" "backtracking" "auto"])));
parser.addParameter("MPLBacktrackingFactor", 2, @(x) isscalar(x) && x > 1);
parser.addParameter("MPLMaxBacktracks", 12, @(x) isscalar(x) && x >= 1 && fix(x) == x);
parser.addParameter("MPLAutoExpansion", 1.25, @(x) isscalar(x) && x > 1);
parser.addParameter("MPLStepMinScale", 0.5, @(x) isscalar(x) && x > 0);
parser.addParameter("MPLStepMaxScale", 16, @(x) isscalar(x) && x >= 1);
parser.addParameter("MPLOracleRefitFminconPolish", false, ...
    @(x) islogical(x) && isscalar(x));
parser.addParameter("CollectDiagnostics", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("CollectFunctionFits", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("FunctionFitGridSize", 101, @(x) isscalar(x) && x >= 20);
parser.addParameter("BaselineBinCount", [], @(x) isempty(x) || ...
    (isscalar(x) && x >= 2 && fix(x) == x));
parser.addParameter("BaselineRoughnessScale", 0.1, ...
    @(x) isscalar(x) && x >= 0);
parser.addParameter("SplineBasisScaling", "center", @(x) any(strcmpi(string(x), ["zscore" "center" "whiten" "none"])));
parser.addParameter("SplineThresholdStage", "iterative", @(x) any(strcmpi(string(x), ["iterative" "final"])));
parser.addParameter("SplineSelectionMetric", "effect", @(x) any(strcmpi(string(x), ["coef" "effect"])));
parser.addParameter("SplineGroupPenaltyScale", 15, @(x) isscalar(x) && x >= 0);
parser.addParameter("SplineGroupPenaltyScaleFailure", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("SplineGroupPenaltyScaleCensoring", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("SplineGroupPenaltyType", "group_mic", @(x) any(strcmpi(string(x), ["group_lasso" "group_mic" "group_mcp" "group_iht"])));
parser.addParameter("SplineGroupPenaltyNormalization", "sum", @(x) any(strcmpi(string(x), ["mean" "sum"])));
parser.addParameter("SplineGroupMICSharpness", 90, @(x) isscalar(x) && x > 0);
parser.addParameter("SplineGroupMICContinuation", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("SplineGroupMICSharpnessPath", [], @(x) isempty(x) || (isnumeric(x) && all(x > 0)));
parser.addParameter("SplineGroupMICForm", "smooth", @(x) any(strcmpi(string(x), ["smooth" "kinked"])));
parser.addParameter("SplineGroupKinkCoefficient", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("SplineGroupMCPConcavity", 3, @(x) isscalar(x) && x > 1);
parser.addParameter("SplineGroupIHTThreshold", 0.1, @(x) isscalar(x) && x >= 0);
parser.addParameter("SplineGroupIHTThresholdFailure", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("SplineGroupIHTThresholdCensoring", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("SplineGroupIHTWarmupIterations", 10, @(x) isscalar(x) && x >= 0 && fix(x) == x);
parser.addParameter("SplineGroupIHTContinuationIterations", 0, ...
    @(x) isscalar(x) && x >= 0 && fix(x) == x);
parser.addParameter("SplineGroupIHTStartFraction", 1, ...
    @(x) isscalar(x) && x > 0 && x <= 1);
parser.addParameter("PostSelectionRefitType", "spline", ...
    @(x) any(strcmpi(string(x), ...
    ["spline" "joint_coefficients" "joint_oracle" "joint_oracle_alpha"])));
parser.addParameter("PostSelectionRefitMaxIter", 30, ...
    @(x) isscalar(x) && x >= 1 && fix(x) == x);
% v2 options
parser.addParameter("SplineOrder", 3, @(x) isscalar(x) && any(x == [3 4]));
parser.addParameter("UnpenalizedClinical", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("MergeDuplicateAdditive", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("SubsampleSeed", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("SubsampleFraction", 0.8, @(x) isscalar(x) && x > 0 && x <= 1);
parser.parse(varargin{:});
opts = parser.Results;
if any(strcmpi(opts.LinearMICForm,["scad" "mcp"])) && ...
        opts.LinearKinkCoefficient <= 0
    error('PLAPH:LinearPenaltySlope', ...
        ['LinearKinkCoefficient must be positive for SCAD/MCP because it ' ...
        'sets their local slope relative to the MIC complexity weight.']);
end
if opts.Tau < 0 && any(strcmpi(string(opts.FitCopula), ...
        ["Clayton" "Gumbel"]))
    error('PLAPH:NegativeTauUnsupportedFamily', ...
        ['Negative fixed association requires a two-sided copula family. ' ...
        'Use Frank or Gaussian, or implement and validate a rotated family.']);
end
margin_specific_mic = ~isempty(opts.MICPenaltyScaleFailure) || ...
    ~isempty(opts.MICPenaltyScaleCensoring);
if isempty(opts.MICPenaltyScaleFailure)
    opts.MICPenaltyScaleFailure = opts.MICPenaltyScale;
end
if isempty(opts.MICPenaltyScaleCensoring)
    opts.MICPenaltyScaleCensoring = opts.MICPenaltyScale;
end
opts.MICPenaltyMarginSpecific = margin_specific_mic;
opts.WorkingSetDir = char(opts.WorkingSetDir);
opts.OutputDir = char(opts.OutputDir);
opts.MICPenaltyNormalization = lower(string(opts.MICPenaltyNormalization));
opts.LinearMICForm = lower(string(opts.LinearMICForm));
opts.LinearL1Type = lower(string(opts.LinearL1Type));
opts.InitialLinearStrategy = lower(string(opts.InitialLinearStrategy));
opts.InitialSplineWeighting = lower(string(opts.InitialSplineWeighting));
opts.InitialSplineStrategy = lower(string(opts.InitialSplineStrategy));
opts.MPLDecomposablePilotStart = lower(string(opts.MPLDecomposablePilotStart));
opts.MPLDecomposablePilotWeighting = lower(string(opts.MPLDecomposablePilotWeighting));
opts.MPLStepStrategy = lower(string(opts.MPLStepStrategy));
opts.SplineBasisScaling = lower(string(opts.SplineBasisScaling));
opts.SplineThresholdStage = lower(string(opts.SplineThresholdStage));
opts.SplineSelectionMetric = lower(string(opts.SplineSelectionMetric));
opts.SplineGroupPenaltyType = lower(string(opts.SplineGroupPenaltyType));
opts.SplineGroupPenaltyNormalization = lower(string(opts.SplineGroupPenaltyNormalization));
opts.SplineGroupMICForm = lower(string(opts.SplineGroupMICForm));
opts.PostSelectionRefitType = lower(string(opts.PostSelectionRefitType));
if isempty(opts.ThresholdLinearFailure)
    opts.ThresholdLinearFailure = opts.ThresholdLinear;
end
if isempty(opts.ThresholdLinearCensoring)
    opts.ThresholdLinearCensoring = opts.ThresholdLinear;
end
if isempty(opts.SplineGroupPenaltyScaleFailure)
    opts.SplineGroupPenaltyScaleFailure = opts.SplineGroupPenaltyScale;
end
if isempty(opts.SplineGroupPenaltyScaleCensoring)
    opts.SplineGroupPenaltyScaleCensoring = opts.SplineGroupPenaltyScale;
end
if isempty(opts.SplineGroupIHTThresholdFailure)
    opts.SplineGroupIHTThresholdFailure = opts.SplineGroupIHTThreshold;
end
if isempty(opts.SplineGroupIHTThresholdCensoring)
    opts.SplineGroupIHTThresholdCensoring = opts.SplineGroupIHTThreshold;
end
if opts.SplineGroupPenaltyType == "group_iht"
    if opts.SplineBasisScaling ~= "center"
        error('PLAPH:GroupIHTRequiresCenteredBasis', ...
            'SplineGroupPenaltyType="group_iht" requires SplineBasisScaling="center".');
    end
    if opts.ThresholdSpline ~= 0
        error('PLAPH:GroupIHTDoubleThreshold', ...
            'Set ThresholdSpline=0 when using group-IHT; kappa is the sole spline-selection threshold.');
    end
end
if opts.SplineGroupMICContinuation && isempty(opts.SplineGroupMICSharpnessPath)
    opts.SplineGroupMICSharpnessPath = unique([2 4 8 16 32 opts.SplineGroupMICSharpness], "stable");
end
end


function [source_p, source_q] = source_dimensions_for_label(label, opts)
if ~isempty(opts.SourceP) && ~isempty(opts.SourceQ)
    source_p = opts.SourceP;
    source_q = opts.SourceQ;
    return;
end
label = string(label);
if label == "primary_p500_q200"
    source_p = 500;
    source_q = 200;
elseif label == "fallback_p500_q100"
    source_p = 500;
    source_q = 100;
elseif label == "sensitivity_p1000_q100"
    source_p = 1000;
    source_q = 100;
else
    source_p = opts.P;
    source_q = opts.Q;
end
if ~isempty(opts.SourceP)
    source_p = opts.SourceP;
end
if ~isempty(opts.SourceQ)
    source_q = opts.SourceQ;
end
end
