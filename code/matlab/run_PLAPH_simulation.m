function sim = run_PLAPH_simulation(experiment, varargin)
%RUN_PLAPH_SIMULATION Main entry point for the revised PLAPH simulations.
%
% This function is the single formal entry point for all simulation modules
% in the paper. It should only orchestrate experiments: define scenarios,
% set seeds, call the single-scenario runner, and save outputs. The lower
% level data-generation, MPL, PL, oracle, and summary functions should stay
% outside this file.
%
% Examples
%   sim = run_PLAPH_simulation("calibration", "DryRun", true);
%   sim = run_PLAPH_simulation("main", "R", 500, "DryRun", true);
%   sim = run_PLAPH_simulation("copula_correct", "DryRun", true);
%   sim = run_PLAPH_simulation("all", "DryRun", true);
%
% Inputs
%   experiment:
%       "calibration"    Experiment 0, calibration and QC.
%       "main"           Experiment 1, Frank correct-copula main study.
%       "copula_correct" Experiment 2, correct family across copulas.
%       "tau_profile"    Experiment 3, fixed tau versus profile tau.
%       "copula_mis"     Experiment 4, copula misspecification sensitivity.
%       "highdim"        Experiment 5, high-dimensional stress test.
%       "ultrahighdim_main_smoke"
%                        Moderate p,q>n candidate, excluded from "all".
%       "ultrahighdim_balanced_smoke"
%                        Moderate p,q>n candidate with balanced f1 signal.
%       "ultrahighdim_strong_smoke"
%                        Strong-signal partner for the weak E8 setting.
%       "sample_size_feasibility"
%                        Reduced n=400/600, p=q=800 validation block.
%       "informative_discrimination"
%                        Distinct event/censoring supports for Gate 3.
%       "informative_discrimination_calibration"
%                        Low-dimensional data-only censoring calibration.
%       "all"            Build all formal scenario blocks.
%
% Name-value options
%   Experiment-specific defaults
%                  "main" uses MICPenaltyScale=10, InitialLinearCap=Inf,
%                  and SplineGroupPenaltyScale=15. "calibration" uses
%                  MICPenaltyScale=50, InitialLinearCap=0.75, and
%                  SplineGroupPenaltyScale=50 unless overridden by caller.
%   R              Override the default number of Monte Carlo replications.
%   GenerateCopulaOverride
%                  Replace the scenario generating family after construction.
%   FitCopulaOverride
%                  Replace the scenario working family after construction.
%   ProfileTauOverride
%                  Replace the scenario fixed/profile association indicator.
%   BaseSeed       First random seed used by scenario/replication runners.
%   OutputDir      Directory for output tables and MAT files.
%   DryRun         If true, only build the scenario table.
%   UseParallel    If true, the future scenario runner may use parfor.
%   SaveCheckpoints
%                  If true, save scenario-level MAT/CSV checkpoints after
%                  each completed scenario in non-dry runs.
%   MaxScenarios   Optional limit for smoke tests.
%   ScenarioIndices
%                  Optional row indices of the scenario table to run.
%   MaxOuterIter   Optional MPL outer-iteration limit for smoke tests.
%   ThresholdLinear
%                  Small numerical cleanup threshold for linear selection.
%   ThresholdLinearFailure
%                  Optional hard threshold for failure-margin linear selection.
%   ThresholdLinearCensoring
%                  Optional hard threshold for censoring-margin linear selection.
%   MICPenaltyScale
%                  Multiplier for the linear MIC/SELO penalty. The current
%                  main-study default is 10; Experiment 0 may use a stronger
%                  value as a small-sample sensitivity/calibration setting.
%   MICPenaltyNormalization
%                  "sum" uses the legacy original-code scale log(n0);
%                  "mean" uses the theory-aligned average-likelihood scale
%                  log(n0)/n.
%   MICSharpness   SELO/MIC sharpness parameter a in tanh(a*x^2).
%   LinearMICForm  "smooth" uses tanh(a*x^2); "kinked" uses
%                  tanh(a*x^2 + b*abs(x)).
%   LinearKinkCoefficient
%                  Internal absolute-value coefficient b for kinked MIC.
%   LinearL1Type   Prototype hybrid linear penalty: "none", "lasso", or
%                  "adaptive". Defaults to "none" and therefore leaves all
%                  formal archived results unchanged.
%   LinearL1Scale  L1 penalty level on the mean-score scale. Zero disables it.
%   LinearRidgeScale
%                  Optional ridge level in the prototype hybrid penalty.
%   AdaptiveWeightGamma, AdaptiveWeightEpsilon, AdaptiveWeightCap
%                  Adaptive-L1 weight controls based on the warm-start
%                  coefficients.
%   InitialLinearStrategy
%                  "legacy" reproduces the original unpenalized linear
%                  warmup; "zero" starts MIC/SELO from zero; "lasso" uses
%                  margin-specific Cox-L1 pilot estimates; "continuation"
%                  uses a staged MIC/SELO path for diagnostics.
%   InitialLinearPilotScale
%                  Score-rate multiplier for the diagnostic Cox-L1 pilot.
%   InitialLinearCap
%                  Optional cap on warmup linear-coefficient magnitudes before
%                  MIC/SELO. Inf reproduces the original warmup.
%   InitialSplineWeighting
%                  "none" uses ordinary margin-specific partial-likelihood
%                  spline warm starts. "legacy_tail" reproduces the archived
%                  DGP-specific tail weights on X1 and X2 and should be used
%                  only to reproduce historical simulation output.
%   InitialSplineStrategy
%                  "legacy" preserves the unpenalized spline warm start.
%                  "group_lasso" uses a data-driven margin-specific Cox
%                  fitted-effect group-lasso pilot on the screened set.
%   InitialSplinePilotScale
%                  Score-rate multiplier for the optional spline group-lasso
%                  pilot. It has no effect under the default legacy strategy.
%   DiscriminatingFailureScale, DiscriminatingCensoringScale
%                  Baseline time scales used only by the distinct-support
%                  informative-discrimination DGP. The censoring scale can be
%                  scalar or one predeclared value per generated scenario.
%   DiscriminatingCensoringFunctionDesign
%                  "current" retains the distinct-support E10 censoring
%                  functions. "shared_opposite" is a development diagnostic
%                  with the four failure functions entering censoring with
%                  opposite signs. "shared_opposite_bounded" additionally
%                  replaces the unbounded Gaussian quadratic component by a
%                  bounded, normalized quadratic component.
%                  "shared_balanced_bounded" uses the same bounded convex
%                  quadratic in both margins while retaining opposite signs
%                  for the other shared active functions. None of these
%                  diagnostic designs is part of the frozen R=200 archive.
%   MPLStepScale   Multiplier for MPL coefficient step denominators. Values
%                  above 1 damp the MPL outer update for diagnostics.
%   MPLDecomposablePilotIterations
%                  Number of initial joint observed-likelihood iterations
%                  using scalar L1 and fitted-effect group-lasso penalties.
%                  Zero, the default, preserves the existing MPL path.
%   MPLDecomposablePilotLinearScale, MPLDecomposablePilotGroupScale
%                  Score-rate multipliers for the optional joint pilot.
%   MPLDecomposablePilotStart
%                  "zero" starts the optional joint pilot at zero regression
%                  effects; "initial" uses the existing margin warm starts.
%   MPLDecomposablePilotWeighting
%                  "none" uses uniform pilot penalties; "adaptive" weights
%                  scalar and function groups from the existing warm starts.
%   CollectDiagnostics
%                  If true, collect one-rep low-level diagnostics for data,
%                  initialization, and the first MPL baseline update.
%   CollectFunctionFits
%                  If true, store fitted additive-function curves on fixed
%                  grids for manuscript recovery plots.
%   FunctionFitGridSize
%                  Number of grid points used for stored function curves.
%   PreScreeningType
%                  "none" preserves the archived path. The recommended
%                  p,q>n diagnostic is "dual_margin_rank_fused_information",
%                  which rank-fuses margin-specific Cox information with
%                  model-free observed-outcome screening.
%   ScreenLinearSize, ScreenAdditiveSize
%                  Predeclared common working-set sizes. Both are required
%                  when screening is enabled and must not be tuned by truth.
%   ScreenAdditiveRank
%                  Rank of the screening-only spline projection. The final
%                  PLAPH fit still uses the full estimation basis.
%   PostSelectionRefitType
%                  "spline" preserves archived behavior. "joint_coefficients"
%                  refits selected linear and spline coefficients without
%                  penalties while holding the final baseline and copula fixed.
%                  "joint_oracle" additionally refits both nonnegative
%                  baseline hazards while holding the working copula fixed.
%                  "joint_oracle_alpha" also updates the copula association
%                  inside the accepted oracle-refit iterations.
%   SplineGroupIHTContinuationIterations, SplineGroupIHTStartFraction
%                  Optional finite threshold-continuation initialization.
%                  Defaults 0 and 1 preserve the archived immediate switch.
%   SplineBasisScaling
%                  "center" only centers B-spline columns and is the formal
%                  default; "zscore" reproduces the original unstable scaling;
%                  "whiten" standardizes each centered spline block by its
%                  empirical Gram matrix; "none" leaves raw bases.
%   SplineThresholdStage
%                  "iterative" thresholds spline groups in each CD iteration;
%                  "final" thresholds only before post-selection refit.
%   SplineSelectionMetric
%                  "coef" selects by coefficient group norm; "effect" selects
%                  by fitted-effect RMS norm sqrt(mean((B_j gamma_j)^2)).
%   ThresholdSpline
%                  Small numerical cleanup threshold for spline groups after
%                  group-MIC shrinkage. It is not used as the main selection
%                  device in formal simulations.
%   SplineGroupPenaltyScale
%                  Scale for the selected group penalty on spline blocks.
%                  Zero disables the group penalty for diagnostics.
%   SplineGroupPenaltyType
%                  "group_mic" uses an adaptive MIC/SELO-style group LQA
%                  shrinkage; "group_mcp" uses the exact block MCP proximal
%                  map; "group_iht" applies exact block hard thresholding on
%                  whitened spline blocks; "group_lasso" is retained for
%                  diagnostics.
%   SplineGroupPenaltyNormalization
%                  "sum" uses the legacy summed-likelihood scale;
%                  "mean" divides by n for the average-likelihood objective.
%   SplineGroupMICSharpness
%                  Sharpness parameter for the group-MIC spline penalty.
%   SplineGroupMICContinuation
%                  If true, apply group-MIC through a staged sharpness path.
%   SplineGroupMICSharpnessPath
%                  Numeric path for continuation; empty uses an automatic path.
%   SplineGroupMCPConcavity
%                  Concavity gamma for group MCP; the default is 3.
%   SplineGroupIHTThreshold
%                  Common group-IHT threshold kappa. Failure- and censoring-
%                  margin values can be supplied separately. The group-IHT
%                  branch uses the centered-basis empirical Gram metric and
%                  requires SplineBasisScaling="center" and
%                  ThresholdSpline=0.
%   SplineGroupIHTWarmupIterations
%                  Number of group-MIC continuation iterations used only to
%                  initialize the nonconvex group-IHT search; zero uses a direct
%                  zero-start IHT path.
%   InitialLinearWarmup
%                  If true, reproduce the original unpenalized linear warmup
%                  before MIC/SELO. Retained for backward compatibility.
%
% Output
%   sim.scenarios  Scenario table.
%   sim.options    Parsed options.
%   sim.results    Scenario-level results when DryRun is false.

if nargin < 1 || strlength(string(experiment)) == 0
    experiment = "main";
end

raw_options = varargin;
opts = parse_options(raw_options{:});
experiment = lower(string(experiment));
opts = apply_experiment_defaults(opts, experiment, raw_options);

scenarios = make_scenarios(experiment, opts);
if ~isempty(opts.ScenarioIndices)
    scenarios = scenarios(opts.ScenarioIndices, :);
end

if ~isempty(opts.MaxScenarios)
    scenarios = scenarios(1:min(opts.MaxScenarios, height(scenarios)), :);
end
if strlength(opts.GenerateCopulaOverride) > 0
    scenarios.generate_copula(:) = opts.GenerateCopulaOverride;
end
if strlength(opts.FitCopulaOverride) > 0
    scenarios.fit_copula(:) = opts.FitCopulaOverride;
end
if ~isempty(opts.ProfileTauOverride)
    scenarios.profile_tau(:) = opts.ProfileTauOverride;
end
if ~isempty(opts.FitTauOverride)
    % v2 ablation: fit with a fixed working association (e.g. 0 =
    % independence) while data are generated at scenarios.tau.
    scenarios.generate_tau = scenarios.tau;
    scenarios.fit_tau = repmat(opts.FitTauOverride,height(scenarios),1);
end
if strlength(opts.GenerateCopulaOverride) > 0 || ...
        strlength(opts.FitCopulaOverride) > 0
    for scenario_index = 1:height(scenarios)
        family_tag = "gen" + scenarios.generate_copula(scenario_index) + ...
            "_fit" + scenarios.fit_copula(scenario_index);
        scenarios.scenario_id(scenario_index) = regexprep( ...
            scenarios.scenario_id(scenario_index), ...
            "gen[^_]+_fit[^_]+",family_tag);
    end
end
if strlength(opts.MethodsOverride) > 0
    scenarios.methods(:) = {char(opts.MethodsOverride)};
end
sim = struct();
sim.experiment = experiment;
sim.options = opts;
sim.scenarios = scenarios;
sim.results = [];

if opts.DryRun
    fprintf("Dry run: built %d scenario(s) for experiment '%s'.\n", height(scenarios), experiment);
    disp(scenarios);
    return;
end

if ~exist(opts.OutputDir, "dir")
    mkdir(opts.OutputDir);
end

results = cell(height(scenarios), 1);
for s = 1:height(scenarios)
    scenario = table2struct(scenarios(s, :));
    fprintf("Running scenario %d/%d: %s\n", s, height(scenarios), scenario.scenario_id);
    results{s} = plaph_run_one_scenario(scenario, opts);
    % v2: plain-text per-replication records for aggregation without HDF5.
    writetable(results{s}.replicates, fullfile(opts.OutputDir, ...
        sprintf("scenario_%03d_replicates.csv", s)));
    if isfield(results{s},"linear_penalty_tuning")
        tuning_rows = results{s}.linear_penalty_tuning;
        filled = find(~cellfun(@isempty,tuning_rows));
        tuning_rows = tuning_rows(filled);
        if ~isempty(tuning_rows)
            for t = 1:numel(tuning_rows)
                tuning_rows{t}.replication = repmat(filled(t),height(tuning_rows{t}),1);
                tuning_rows{t}.seed = repmat(results{s}.replicates.seed(filled(t)), ...
                    height(tuning_rows{t}),1);
            end
            writetable(vertcat(tuning_rows{:}), fullfile(opts.OutputDir, ...
                sprintf("scenario_%03d_tuning_path.csv", s)));
        end
    end
    if opts.SaveCheckpoints
        save_scenario_checkpoint(opts.OutputDir, s, scenarios(s, :), results{s});
        sim.results = results;
        save(fullfile(opts.OutputDir, "PLAPH_simulation_checkpoint.mat"), "sim", "-v7.3");
    end
end

sim.results = results;
save(fullfile(opts.OutputDir, "PLAPH_simulation_results.mat"), "sim", "-v7.3");
end


function save_scenario_checkpoint(output_dir, scenario_index, scenario_row, result)
checkpoint_name = sprintf("scenario_%03d_checkpoint.mat", scenario_index);
scenario_checkpoint = scenario_row;
save(fullfile(output_dir, checkpoint_name), "scenario_checkpoint", "result", "-v7.3");

summary_path = fullfile(output_dir, sprintf("scenario_%03d_summary.csv", scenario_index));
if istable(result.summary)
    summary_table = result.summary;
else
    summary_table = cell2table(struct2cell(result.summary)', ...
        "VariableNames", string(fieldnames(result.summary))');
end
writetable(summary_table, summary_path);
end


function opts = parse_options(varargin)
parser = inputParser;
parser.addParameter("R", [], @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("BaseSeed", 20260705, @(x) isscalar(x) && x > 0);
parser.addParameter("OutputDir", fullfile(pwd, "simulation_outputs"), @(x) ischar(x) || isstring(x));
parser.addParameter("DryRun", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("UseParallel", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("SaveCheckpoints", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("MaxScenarios", [], @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("ScenarioIndices", [], @(x) isempty(x) || isnumeric(x));
parser.addParameter("MethodsOverride", "", @(x) ischar(x) || isstring(x));
parser.addParameter("GenerateCopulaOverride", "", ...
    @(x) ischar(x) || isstring(x));
parser.addParameter("FitCopulaOverride", "", ...
    @(x) ischar(x) || isstring(x));
parser.addParameter("ProfileTauOverride", [], ...
    @(x) isempty(x) || (islogical(x) && isscalar(x)));
parser.addParameter("MaxOuterIter", 300, @(x) isscalar(x) && x > 0);
parser.addParameter("NumInternalKnots", 8, @(x) isscalar(x) && x > 0);
parser.addParameter("ThresholdLinear", 5e-4, @(x) isscalar(x) && x >= 0);
parser.addParameter("ThresholdLinearFailure", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("ThresholdLinearCensoring", [], @(x) isempty(x) || (isscalar(x) && x >= 0));
parser.addParameter("ThresholdSpline", 1e-3, @(x) isscalar(x) && x >= 0);
parser.addParameter("LinearPenaltyScale", 10, @(x) isscalar(x) && x > 0);
parser.addParameter("MICPenaltyScale", [], @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("MICPenaltyScaleFailure", [], ...
    @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("MICPenaltyScaleCensoring", [], ...
    @(x) isempty(x) || (isscalar(x) && x > 0));
parser.addParameter("MICPenaltyNormalization", "sum", @(x) any(strcmpi(string(x), ["mean" "sum"])));
parser.addParameter("MICSharpness", 90, @(x) isscalar(x) && x > 0);
parser.addParameter("LinearMICForm", "smooth", @(x) any(strcmpi(string(x), ["smooth" "kinked"])));
parser.addParameter("LinearKinkCoefficient", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("LinearL1Type", "none", @(x) any(strcmpi(string(x), ["none" "lasso" "adaptive"])));
parser.addParameter("LinearL1Scale", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("LinearRidgeScale", 0, @(x) isscalar(x) && x >= 0);
parser.addParameter("AdaptiveWeightGamma", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("AdaptiveWeightEpsilon", 0.05, @(x) isscalar(x) && x > 0);
parser.addParameter("AdaptiveWeightCap", 100, @(x) isscalar(x) && x >= 1);
parser.addParameter("InitialLinearStrategy", "legacy", @(x) any(strcmpi(string(x), ...
    ["legacy" "zero" "lasso" "continuation"])));
parser.addParameter("InitialLinearPilotScale", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("InitialLinearCap", Inf, @(x) isscalar(x) && x > 0);
parser.addParameter("InitialLinearWarmup", true, @(x) islogical(x) && isscalar(x));
parser.addParameter("InitialSplineWeighting", "none", @(x) any(strcmpi(string(x), ...
    ["none" "legacy_tail"])));
parser.addParameter("InitialSplineStrategy", "legacy", @(x) any(strcmpi(string(x), ...
    ["legacy" "group_lasso"])));
parser.addParameter("InitialSplinePilotScale", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("DiscriminatingFailureScale", 1, @(x) isscalar(x) && x > 0);
parser.addParameter("DiscriminatingCensoringScale", [], ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && all(x > 0)));
parser.addParameter("DiscriminatingTauValues", [0 0.2 0.5 0.8], ...
    @(x) isnumeric(x) && isvector(x) && all(abs(x) < 1));
parser.addParameter("DiscriminatingSampleSizes", [500 700], ...
    @(x) isnumeric(x) && isvector(x) && ~isempty(x) && ...
    all(x >= 1) && all(fix(x) == x));
parser.addParameter("DiscriminatingDimensions", 800, ...
    @(x) isnumeric(x) && isvector(x) && ~isempty(x) && ...
    all(x >= 1) && all(fix(x) == x));
parser.addParameter("DiscriminatingCensorRates", [0.30 0.50], ...
    @(x) isnumeric(x) && isvector(x) && ~isempty(x) && ...
    all(x > 0 & x < 1));
parser.addParameter("DiscriminatingAdditiveSignalScale", 1.5, ...
    @(x) isscalar(x) && x > 0);
parser.addParameter("DiscriminatingCensoringFunctionDesign", "current", ...
    @(x) any(strcmpi(string(x), ...
    ["current" "shared_opposite" "shared_opposite_bounded" ...
    "shared_balanced_bounded"])));
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
parser.addParameter("MPLOracleRefitFminconPolish", false, ...
    @(x) islogical(x) && isscalar(x));
parser.addParameter("MPLOracleRefitFminconAlgorithm", "sqp", ...
    @(x) any(strcmpi(string(x), ["sqp" "interior-point"])));
parser.addParameter("MPLOracleRefitFminconScaleProblem", false, ...
    @(x) islogical(x) && isscalar(x));
parser.addParameter("MPLOracleRefitGradientAudit", false, ...
    @(x) islogical(x) && isscalar(x));
parser.addParameter("MPLCollectModelBasedInference", false, ...
    @(x) islogical(x) && isscalar(x));
parser.addParameter("MPLAutoExpansion", 1.25, @(x) isscalar(x) && x > 1);
parser.addParameter("MPLStepMinScale", 0.5, @(x) isscalar(x) && x > 0);
parser.addParameter("MPLStepMaxScale", 16, @(x) isscalar(x) && x >= 1);
parser.addParameter("CollectDiagnostics", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("CollectFunctionFits", false, @(x) islogical(x) && isscalar(x));
parser.addParameter("FunctionFitGridSize", 101, @(x) isscalar(x) && x >= 20);
parser.addParameter("BaselineBinCount", [], @(x) isempty(x) || ...
    (isscalar(x) && x >= 2 && fix(x) == x));
parser.addParameter("BaselineRoughnessScale", 0.1, ...
    @(x) isscalar(x) && x >= 0);
parser.addParameter("SplineBasisScaling", "center", @(x) any(strcmpi(string(x), ["zscore" "center" "whiten" "none"])));
parser.addParameter("SplineThresholdStage", "iterative", @(x) any(strcmpi(string(x), ["iterative" "final"])));
parser.addParameter("SplineSelectionMetric", "coef", @(x) any(strcmpi(string(x), ["coef" "effect"])));
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
parser.addParameter("SplineGroupIHTContinuationIterations", 0, @(x) isscalar(x) && x >= 0 && fix(x) == x);
parser.addParameter("SplineGroupIHTStartFraction", 1, @(x) isscalar(x) && x > 0 && x <= 1);
parser.addParameter("PostSelectionRefitType", "spline", @(x) any(strcmpi(string(x), ...
    ["spline" "joint_coefficients" "joint_oracle" "joint_oracle_alpha"])));
parser.addParameter("PostSelectionRefitMaxIter", 30, ...
    @(x) isscalar(x) && x >= 1 && fix(x) == x);
parser.addParameter("PreScreeningType", "none", @(x) any(strcmpi(string(x), ...
    ["none" "dual_margin_cox_score" "dual_margin_cox_information" ...
    "dual_margin_rank_fused_information" ...
    "joint_distance_correlation" "joint_fused_observed" ...
    "oracle_support_diagnostic"])));
parser.addParameter("ScreenLinearSize", [], @(x) isempty(x) || (isscalar(x) && x >= 1 && fix(x) == x));
parser.addParameter("ScreenAdditiveSize", [], @(x) isempty(x) || (isscalar(x) && x >= 1 && fix(x) == x));
parser.addParameter("ScreenAdditiveRank", 4, @(x) isscalar(x) && x >= 1 && fix(x) == x);
% v2 options
parser.addParameter("DiscriminatingRho", 0.25, @(x) isscalar(x) && x >= 0 && x < 1);
parser.addParameter("SplineOrder", 3, @(x) isscalar(x) && any(x == [3 4]));
parser.addParameter("TuningMode", "fixed", @(x) any(strcmpi(string(x), ["fixed" "grid_ebic" "grid_bic"])));
parser.addParameter("TuningGridCL", [0.5 1.0 1.5 2.0 2.5], @(x) isnumeric(x) && all(x > 0));
parser.addParameter("TuningGridKappa", [0.03 0.05 0.075 0.10 0.15], @(x) isnumeric(x) && all(x > 0));
parser.addParameter("PLVersion", "v2", @(x) any(strcmpi(string(x), ["v2" "legacy"])));
parser.addParameter("ExportDataDir", "", @(x) ischar(x) || isstring(x));
parser.addParameter("FitTauOverride", [], @(x) isempty(x) || (isscalar(x) && abs(x) < 1));
parser.addParameter("LinearPenaltyTuningCriterion", "none", @(x) any(strcmpi(string(x), ["none" "bic" "ebic"])));
parser.addParameter("LinearPenaltyMultiplierGrid", 25.6685, @(x) isnumeric(x));
parser.addParameter("LinearPenaltyEBICGamma", 1, @(x) isscalar(x) && x >= 0);
parser.parse(varargin{:});
opts = parser.Results;
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
opts.OutputDir = char(opts.OutputDir);
opts.MethodsOverride = string(opts.MethodsOverride);
opts.GenerateCopulaOverride = canonical_copula_override( ...
    opts.GenerateCopulaOverride);
opts.FitCopulaOverride = canonical_copula_override(opts.FitCopulaOverride);
margin_specific_mic = ~isempty(opts.MICPenaltyScaleFailure) || ...
    ~isempty(opts.MICPenaltyScaleCensoring);
if isempty(opts.MICPenaltyScale)
    opts.MICPenaltyScale = opts.LinearPenaltyScale;
end
if isempty(opts.MICPenaltyScaleFailure)
    opts.MICPenaltyScaleFailure = opts.MICPenaltyScale;
end
if isempty(opts.MICPenaltyScaleCensoring)
    opts.MICPenaltyScaleCensoring = opts.MICPenaltyScale;
end
opts.MICPenaltyMarginSpecific = margin_specific_mic;
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
opts.PostSelectionRefitType = lower(string(opts.PostSelectionRefitType));
opts.PreScreeningType = lower(string(opts.PreScreeningType));
if ~any(opts.PreScreeningType == ["none" "oracle_support_diagnostic"]) && ...
        (isempty(opts.ScreenLinearSize) || isempty(opts.ScreenAdditiveSize))
    error('PLAPH:MissingScreeningSize', ...
        'ScreenLinearSize and ScreenAdditiveSize must be prespecified when screening is enabled.');
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
if ~opts.InitialLinearWarmup && opts.InitialLinearStrategy == "legacy"
    opts.InitialLinearStrategy = "zero";
end
if isempty(opts.ThresholdLinearFailure)
    opts.ThresholdLinearFailure = opts.ThresholdLinear;
end
if isempty(opts.ThresholdLinearCensoring)
    opts.ThresholdLinearCensoring = opts.ThresholdLinear;
end
end


function opts = apply_experiment_defaults(opts, experiment, raw_options)
%APPLY_EXPERIMENT_DEFAULTS Override defaults only when not set by the caller.
if experiment ~= "calibration"
    return;
end

has_mic = has_name_value(raw_options, "MICPenaltyScale");
has_linear = has_name_value(raw_options, "LinearPenaltyScale");
if ~has_mic && ~has_linear
    opts.MICPenaltyScale = 50;
end
if ~has_name_value(raw_options, "InitialLinearCap")
    opts.InitialLinearCap = 0.75;
end
if ~has_name_value(raw_options, "SplineGroupPenaltyScale")
    opts.SplineGroupPenaltyScale = 50;
end
end


function tf = has_name_value(raw_options, name)
if isempty(raw_options)
    tf = false;
    return;
end
option_names = string(raw_options(1:2:end));
tf = any(strcmpi(option_names, string(name)));
end


function scenarios = make_scenarios(experiment, opts)
switch experiment
    case "calibration"
        scenarios = scenarios_calibration(opts);
    case "main"
        scenarios = scenarios_main(opts);
    case "copula_correct"
        scenarios = scenarios_copula_correct(opts);
    case "tau_profile"
        scenarios = scenarios_tau_profile(opts);
    case "copula_mis"
        scenarios = scenarios_copula_mis(opts);
    case "highdim"
        scenarios = scenarios_highdim(opts);
    case "ultrahighdim_smoke"
        scenarios = scenarios_ultrahighdim_smoke(opts);
    case "ultrahighdim_main_smoke"
        scenarios = scenarios_ultrahighdim_main_smoke(opts);
    case "ultrahighdim_balanced_smoke"
        scenarios = scenarios_ultrahighdim_balanced_smoke(opts);
    case "ultrahighdim_strong_smoke"
        scenarios = scenarios_ultrahighdim_strong_smoke(opts);
    case "sample_size_feasibility"
        scenarios = scenarios_sample_size_feasibility(opts);
    case "informative_discrimination"
        scenarios = scenarios_informative_discrimination(opts);
    case "informative_discrimination_calibration"
        scenarios = scenarios_informative_discrimination_calibration(opts);
    case "all"
        scenarios = [
            scenarios_calibration(opts);
            scenarios_main(opts);
            scenarios_copula_correct(opts);
            scenarios_tau_profile(opts);
            scenarios_copula_mis(opts);
            scenarios_highdim(opts)
        ];
    otherwise
        error("Unknown experiment: %s", experiment);
end
end


function T = scenarios_calibration(opts)
R = default_R(opts, 200);
T = scenario_grid("E0_calibration", R, 300, 20, 20, 0.25, [0.20 0.40], ...
    "Frank", "Frank", [0 0.2 0.5], false, ["MPL_fixed_tau"]);
end


function T = scenarios_main(opts)
R = default_R(opts, 200);
T1 = scenario_grid("E1_main_n500", R, 500, 50, 50, [0.25 0.50], [0.20 0.40], ...
    "Frank", "Frank", [0 0.2 0.5], false, ["MPL_fixed_tau" "PL" "Oracle"]);
T2 = scenario_grid("E1_main_n800", R, 800, 100, 100, [0.25 0.50], [0.20 0.40], ...
    "Frank", "Frank", [0 0.2 0.5], false, ["MPL_fixed_tau" "PL" "Oracle"]);
T = [T1; T2];
end


function T = scenarios_copula_correct(opts)
R = default_R(opts, 200);
families = ["Frank" "Gaussian" "Clayton" "Gumbel"];
T = table();
for f = families
    tau_values = [0.2 0.5];
    T = [T; scenario_grid("E2_correct_" + f, R, 500, 50, 50, 0.25, 0.40, ...
        f, f, tau_values, false, ["MPL_fixed_tau" "PL" "Oracle"])];
end
end


function T = scenarios_tau_profile(opts)
R = default_R(opts, 200);
T_fixed = scenario_grid("E3_tau_fixed", R, 500, 50, 50, 0.25, 0.40, ...
    "Frank", "Frank", [0.2 0.5], false, ["MPL_fixed_tau" "PL" "Oracle"]);
T_profile = scenario_grid("E3_tau_profile", R, 500, 50, 50, 0.25, 0.40, ...
    "Frank", "Frank", [0.2 0.5], true, ["MPL_profile_tau" "PL" "Oracle"]);
T = [T_fixed; T_profile];
end


function T = scenarios_copula_mis(opts)
R = default_R(opts, 200);
generate_families = ["Frank" "Gaussian" "Clayton" "Gumbel"];
fit_families = ["Frank" "Gaussian" "Clayton" "Gumbel"];
T = table();
for g = generate_families
    for f = fit_families
        T = [T; scenario_grid("E4_mis_" + g + "_fit_" + f, R, 500, 50, 50, ...
            0.25, 0.40, g, f, 0.5, false, ["MPL_fixed_tau" "PL" "Oracle"])];
    end
end
end


function T = scenarios_highdim(opts)
R = default_R(opts, 200);
T = scenario_grid("E5_highdim", R, 500, 250, 250, 0.25, 0.40, ...
    "Frank", "Frank", 0.5, false, ["MPL_fixed_tau" "PL" "Oracle"]);
end


function T = scenarios_ultrahighdim_smoke(opts)
% Separate p,q>n diagnostic. It is intentionally excluded from "all" until
% its numerical and statistical behavior has been validated.
R = default_R(opts, 200);
T = scenario_grid("E6_ultrahighdim_smoke", R, 300, 500, 500, 0.25, 0.40, ...
    "Frank", "Frank", 0.5, false, "MPL_fixed_tau");
end


function T = scenarios_ultrahighdim_main_smoke(opts)
% Moderate p,q>n candidate with the original signal specification. It is
% excluded from "all" until screening and second-stage pilots are complete.
R = default_R(opts, 200);
T = scenario_grid("E7_ultrahighdim_main_smoke", R, 500, 600, 600, 0.25, 0.40, ...
    "Frank", "Frank", 0.5, false, "MPL_fixed_tau");
end


function T = scenarios_ultrahighdim_balanced_smoke(opts)
% Preserve the original DGP except for the weakest additive linear component.
R = default_R(opts, 200);
T = scenario_grid("E8_ultrahighdim_balanced_smoke", R, 500, 600, 600, ...
    0.25, 0.40, "Frank", "Frank", 0.5, false, "MPL_fixed_tau");
T.additive_linear_coefficient(:) = 0.8;
end


function T = scenarios_ultrahighdim_strong_smoke(opts)
% Strong beta-min partner fixed before generating its pilot data.
R = default_R(opts, 200);
T = scenario_grid("E9_ultrahighdim_strong_smoke", R, 500, 600, 600, ...
    0.25, 0.40, "Frank", "Frank", 0.5, false, "MPL_fixed_tau");
T.additive_linear_coefficient(:) = 0.8;
T.linear_signal_scale(:) = 1.5;
T.additive_signal_scale(:) = 1.5;
end


function T = scenarios_sample_size_feasibility(opts)
% Reduced validation of the proposed n and ambient-dimension settings.
% This block is excluded from "all" and is not a formal simulation result.
R = default_R(opts, 5);
T = scenario_grid("E10_sample_size_feasibility", R, [400 600], 800, 800, ...
    0.25, 0.50, "Frank", "Frank", 0.5, false, "MPL_fixed_tau");
T.additive_linear_coefficient(:) = 0.8;
T.linear_signal_scale(:) = 1.5;
T.additive_signal_scale(:) = 1.5;
end


function T = scenarios_informative_discrimination(opts)
% Gate 3: the event and censoring margins have different sparse structures.
R = default_R(opts,20);
T = table();
for dimension = opts.DiscriminatingDimensions(:)'
    dimension_grid = scenario_grid("E10_informative_discrimination",R, ...
        opts.DiscriminatingSampleSizes,dimension,dimension, ...
        opts.DiscriminatingRho,opts.DiscriminatingCensorRates,"Frank","Frank", ...
        opts.DiscriminatingTauValues,false, ...
        ["MPL_fixed_tau" "PL" "Oracle"]);
    T = [T; dimension_grid]; %#ok<AGROW>
end
T.additive_linear_coefficient(:) = 0.8;
T.linear_signal_scale(:) = 1.5;
T.additive_signal_scale(:) = opts.DiscriminatingAdditiveSignalScale;
T.dgp_type = repmat("discriminating",height(T),1);
T.censoring_function_design = repmat( ...
    lower(string(opts.DiscriminatingCensoringFunctionDesign)),height(T),1);
T.failure_scale = repmat(opts.DiscriminatingFailureScale,height(T),1);
scales = opts.DiscriminatingCensoringScale(:);
if isempty(scales)
    calibrated = [9.6 9.4 8.9 8.4 1.3 1.3 1.3 1.3]';
    scales = repmat(calibrated,2,1);
elseif isscalar(scales)
    scales = repmat(scales,height(T),1);
elseif numel(scales) ~= height(T)
    error('PLAPH:DiscriminatingCensorScaleSize', ...
        'DiscriminatingCensoringScale must be scalar or have %d values.',height(T));
end
T.censoring_scale = scales;
end


function T = scenarios_informative_discrimination_calibration(opts)
% Large-n, low-dimensional data-only version of the exact Gate 3 DGP.
R = default_R(opts,20);
T = scenario_grid("E10_informative_calibration",R,5000,7,6,0.25, ...
    opts.DiscriminatingCensorRates,"Frank","Frank", ...
    opts.DiscriminatingTauValues,false,"DataOnly");
T.additive_linear_coefficient(:) = 0.8;
T.linear_signal_scale(:) = 1.5;
T.additive_signal_scale(:) = opts.DiscriminatingAdditiveSignalScale;
T.dgp_type = repmat("discriminating",height(T),1);
T.censoring_function_design = repmat( ...
    lower(string(opts.DiscriminatingCensoringFunctionDesign)),height(T),1);
T.failure_scale = repmat(opts.DiscriminatingFailureScale,height(T),1);
scales = opts.DiscriminatingCensoringScale(:);
if isempty(scales)
    error('PLAPH:MissingCalibrationCensorScale', ...
        'The data-only calibration experiment requires an explicit censoring scale.');
elseif isscalar(scales)
    scales = repmat(scales,height(T),1);
elseif numel(scales) ~= height(T)
    error('PLAPH:DiscriminatingCensorScaleSize', ...
        'DiscriminatingCensoringScale must be scalar or have %d values.',height(T));
end
T.censoring_scale = scales;
end


function family = canonical_copula_override(value)
family = string(value);
if strlength(family) == 0
    return;
end
allowed = ["Frank","Gaussian","Clayton","Gumbel"];
match = find(strcmpi(family,allowed),1);
if isempty(match)
    error("PLAPH:CopulaOverride", ...
        "Copula override must be Frank, Gaussian, Clayton, or Gumbel.");
end
family = allowed(match);
end


function R = default_R(opts, fallback)
if isempty(opts.R)
    R = fallback;
else
    R = opts.R;
end
end


function T = scenario_grid(prefix, R, n_values, p_values, q_values, rho_values, cr_values, gen_families, fit_families, tau_values, profile_tau, methods)
rows = {};
for n = n_values
    for p = p_values
        for q = q_values
            for rho = rho_values
                for cr = cr_values
                    for gen_family = gen_families
                        for fit_family = fit_families
                            for tau = tau_values
                                scenario_id = sprintf("%s_n%d_p%d_q%d_rho%.2f_cr%.2f_gen%s_fit%s_tau%.2f", ...
                                    prefix, n, p, q, rho, cr, gen_family, fit_family, tau);
                                method_text = strjoin(cellstr(string(methods)), "|");
                                rows(end+1, :) = {scenario_id, prefix, R, n, p, q, rho, cr, ...
                                    gen_family, fit_family, tau, profile_tau, method_text}; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end
end
T = cell2table(rows, "VariableNames", ...
    ["scenario_id" "experiment" "R" "n" "p" "q" "rho" "censor_rate" ...
     "generate_copula" "fit_copula" "tau" "profile_tau" "methods"]);
end


function out = plaph_run_one_scenario(varargin)
out = plaph_run_one_scenario_impl(varargin{:});
end
