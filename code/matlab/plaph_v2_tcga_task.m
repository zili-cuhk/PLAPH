function out = plaph_v2_tcga_task(task, working_set_dir, out_root, varargin)
%PLAPH_V2_TCGA_TASK One task of the v2 TCGA-LUAD reanalysis.
%
% Tasks (see plaph_v2_tcga_tasks):
%   "grid_<Family>_<tau>"  fixed working association, EBIC over (c_L, kappa)
%   "profile_<Family>"     descriptive profile fit at the tau=0 EBIC tuning
%   "subsample_<k>"        k-th stratified 80% subsample at fixed Frank tau
%                          (tuning fixed at the full-data tau=0 EBIC choice)
%   "kpath_<kappa>"        additive sensitivity at tau=0: c_L at the tau=0
%                          EBIC choice, group threshold kappa in both margins
%   "wexp"                 additive sensitivity at tau=0: tau=0 EBIC tuning,
%                          adaptive group-weight exponent 0.25 instead of 1
%   "dim_<P>_<Q>"          screened-dimension sensitivity: EBIC grid at tau=0
%                          on the first P RNA and Q CNA screened candidates
% All fits: age, male, stage II/III/IV unpenalized in both margins;
% identical CNA profiles merged; same estimator settings as the simulations.

parser = inputParser;
parser.addParameter("SplineOrder", 3, @(x) any(x == [3 4]));
parser.addParameter("GridCL", [0.5 1.0 1.5 2.0 2.5 3.0 4.0 5.0 6.0]);  % extended upward: the
% simulation grid put the TCGA EBIC optimum at its upper edge c_L = 2.5
parser.addParameter("GridKappa", [0.03 0.05 0.075 0.10 0.15]);
parser.addParameter("SubsampleTau", 0, @(x) isscalar(x) && abs(x) < 1);
parser.addParameter("Label", "primary_p500_q200");
parser.addParameter("EBICGamma", 1);
parser.addParameter("GroupWeightGamma", 1);
parser.addParameter("P", 500);
parser.addParameter("Q", 200);
parser.parse(varargin{:});
o = parser.Results;
task = string(task);
out_root = string(out_root);
parts = split(task,"_");
kind = parts(1);

switch kind
    case "grid"
        family = parts(2);
        tau = str2double(parts(3));
        out = ebic_grid(family, tau, false, 0, working_set_dir, ...
            fullfile(out_root,"grid",family + "_tau" + sprintf("%+.2f",tau)), o);
    case "profile"
        family = parts(2);
        chosen = read_choice(fullfile(out_root,"grid","Frank_tau+0.00"));
        out = fit_once(family, 0.2, true, 0, chosen.c_L, chosen.kappa, ...
            working_set_dir, fullfile(out_root,"profile",family), o);
    case "subsample"
        k = str2double(parts(2));
        chosen = read_choice(fullfile(out_root,"grid", ...
            "Frank_tau" + sprintf("%+.2f",o.SubsampleTau)));
        out = fit_once("Frank", o.SubsampleTau, false, 20261020 + k, ...
            chosen.c_L, chosen.kappa, working_set_dir, ...
            fullfile(out_root,"subsample",sprintf("sub_%03d",k)), o);
    case "kpath"
        kappa = str2double(parts(2));
        chosen = read_choice(fullfile(out_root,"grid","Frank_tau+0.00"));
        out = fit_once("Frank", 0, false, 0, chosen.c_L, kappa, working_set_dir, ...
            fullfile(out_root,"additive","kpath_" + sprintf("%.4f",kappa)), o);
    case "wexp"
        chosen = read_choice(fullfile(out_root,"grid","Frank_tau+0.00"));
        o.GroupWeightGamma = 0.25;
        out = fit_once("Frank", 0, false, 0, chosen.c_L, chosen.kappa, working_set_dir, ...
            fullfile(out_root,"additive","weight_exponent_0.25"), o);
    case "dim"
        % First P RNA and Q CNA columns of the primary working set, which is
        % ordered by the preprocessing screen rank.
        o.P = str2double(parts(2)); o.Q = str2double(parts(3));
        out = ebic_grid("Frank", 0, false, 0, working_set_dir, ...
            fullfile(out_root,"dimension",sprintf("p%d_q%d",o.P,o.Q)), o);
    otherwise
        error("PLAPH:TCGATask","Unknown task %s.",task);
end
end


function best = ebic_grid(family, tau, profile, subsample_seed, wsd, out_dir, o)
[gc,gk] = ndgrid(o.GridCL(:), o.GridKappa(:));
rows = cell(numel(gc),1);
for g = 1:numel(gc)
    res = fit_once(family, tau, profile, subsample_seed, gc(g), gk(g), wsd, ...
        fullfile(out_dir, sprintf("cL%.2f_kappa%.2f",gc(g),gk(g))), o);
    rows{g} = res;
end
T = vertcat(rows{:});
valid = T.accepted == 1 & isfinite(T.ebic);
T.selected = false(height(T),1);
if any(valid)
    v = T.ebic; v(~valid) = Inf;
else
    v = T.ebic; % no accepted fit: flagged in the output
end
% Same tie rule as the simulation engine: within max(1e-4,1e-6|min|) of the
% minimum, take the largest c_L, then the largest kappa.
m = min(v);
tie = find(v <= m + max(1e-4, 1e-6*abs(m)));
[~,order] = sortrows([T.c_L(tie) T.kappa(tie)], [-1 -2]);
b = tie(order(1));
T.tie_count = repmat(numel(tie), height(T), 1);
T.at_boundary = repmat(~any(T.c_L(tie) > min(o.GridCL) & T.c_L(tie) < max(o.GridCL) & ...
    T.kappa(tie) > min(o.GridKappa) & T.kappa(tie) < max(o.GridKappa)), height(T), 1);
T.selected(b) = true;
writetable(T, fullfile(out_dir, "ebic_grid.csv"));
best = T(b,:);
writetable(best, fullfile(out_dir, "ebic_choice.csv"));
end


function row = fit_once(family, tau, profile, subsample_seed, cL, kappa, wsd, out_dir, o)
if ~isfolder(out_dir)
    mkdir(out_dir);
end
result = run_PLAPH_real_data(o.Label, ...
    'WorkingSetDir',wsd,'OutputDir',out_dir,'P',o.P,'Q',o.Q, ...
    'FitCopula',family,'Tau',tau,'EstimateTau',profile, ...
    'BaseSeed',20260717,'MaxOuterIter',6000,'NumInternalKnots',4, ...
    'ThresholdLinear',5e-4,'ThresholdSpline',0, ...
    'MICPenaltyScale',cL,'MICPenaltyScaleFailure',cL,'MICPenaltyScaleCensoring',cL, ...
    'MICPenaltyNormalization','mean','MICSharpness',90, ...
    'LinearMICForm','kinked','LinearKinkCoefficient',25.6685, ...
    'InitialLinearStrategy','legacy','InitialLinearCap',Inf, ...
    'InitialSplineStrategy','group_lasso','InitialSplinePilotScale',0.25, ...
    'InitialSplineWeighting','none', ...
    'MPLDecomposablePilotIterations',100,'MPLDecomposablePilotLinearScale',1, ...
    'MPLDecomposablePilotGroupScale',0.5,'MPLDecomposablePilotStart','initial', ...
    'MPLDecomposablePilotWeighting','adaptive','MPLDecomposablePilotGroupWeightGamma',o.GroupWeightGamma, ...
    'MPLStepScale',0.25,'MPLStepStrategy','fixed', ...
    'SplineBasisScaling','center','SplineSelectionMetric','effect', ...
    'SplineGroupPenaltyType','group_iht','SplineGroupPenaltyNormalization','mean', ...
    'SplineGroupIHTThresholdFailure',kappa,'SplineGroupIHTThresholdCensoring',kappa, ...
    'SplineGroupIHTWarmupIterations',0,'SplineGroupIHTContinuationIterations',0, ...
    'SplineGroupIHTStartFraction',1, ...
    'PostSelectionRefitType',ternary(profile,'joint_oracle_alpha','joint_oracle'), ...
    'PostSelectionRefitMaxIter',1000,'MPLOracleRefitFminconPolish',true, ...
    'BaselineBinCount',10,'BaselineRoughnessScale',0, ...
    'CollectDiagnostics',true,'CollectFunctionFits',true, ...
    'SplineOrder',o.SplineOrder,'UnpenalizedClinical',true, ...
    'MergeDuplicateAdditive',true,'SubsampleSeed',subsample_seed);
n_clin = numel(result.unpenalized_clinical);
p_rna = result.p - n_clin;
sT = sum(result.selection_support.linear_failure(n_clin+1:end));
sC = sum(result.selection_support.linear_censoring(n_clin+1:end));
aT = sum(result.selection_support.additive_failure);
aC = sum(result.selection_support.additive_censoring);
model_space = log_binom(p_rna,sT) + log_binom(p_rna,sC) + ...
    log_binom(result.q,aT) + log_binom(result.q,aC);
row = table(string(family), tau, profile, subsample_seed, cL, kappa, ...
    result.n, result.n_failure_events, result.tau_hat, result.total_loglik, ...
    result.effective_df, result.aic, result.bic, ...
    result.bic + 2*o.EBICGamma*model_space, sT, sC, aT, aC, ...
    double(result.selection_gate_pass), double(result.refit_stationarity_pass), ...
    double(result.accepted), result.selection_gate_maximum, ...
    result.refit_stationarity_inf, string(out_dir), ...
    'VariableNames',{'family','tau','profile','subsample_seed','c_L','kappa', ...
    'n','deaths','tau_hat','loglik','df','aic','bic','ebic', ...
    'rna_failure','rna_censoring','cna_failure','cna_censoring', ...
    'selection_gate','refit_pass','accepted','selection_residual', ...
    'refit_residual','output_dir'});
writetable(row, fullfile(out_dir,"fit_row.csv"));
end


function chosen = read_choice(grid_dir)
f = fullfile(grid_dir,"ebic_choice.csv");
if ~isfile(f)
    error("PLAPH:TCGAMissingChoice", ...
        "Run the corresponding grid task first: %s is missing.", f);
end
chosen = readtable(f);
end


function v = log_binom(n,k)
v = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
end


function out = ternary(c,a,b)
if c
    out = a;
else
    out = b;
end
end
