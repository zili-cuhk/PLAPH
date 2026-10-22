function result = plaph_run_one_scenario_impl(scenario, opts)
%PLAPH_RUN_ONE_SCENARIO_IMPL Run one PLAPH simulation scenario.
%
% PLAPH engine v2 (2026-10-20). Features: configurable spline order
% (opts.SplineOrder, default 3 = quadratic), selection-stage baseline
% preconditioner on the mean scale, a first-order selection gate for every
% copula family, step-doubling backtracking, an iteration budget that excludes
% the pilot, the aligned PL comparator, rank-based degrees of freedom, KKT
% residuals in estimable spline coordinates, trapezoidal IMSE, and optional
% EBIC grid tuning over (c_L, kappa).

tic;
if ~isfield(opts,"SplineOrder") || isempty(opts.SplineOrder)
    opts.SplineOrder = 3;
end
plaph_spline_order(opts.SplineOrder);
opts = plaph_v2_default_options(opts);
if isfield(scenario, "real_data") && logical(scenario.real_data)
    result = plaph_run_real_data_impl(scenario, opts);
    return;
end
scenario = normalize_scenario(scenario);
R = scenario.R;
n = scenario.n;
p = scenario.p;
q = scenario.q;
rho = scenario.rho;
generate_tau = scenario.tau;
fit_tau = scenario.tau;
if isfield(scenario,"generate_tau")
    generate_tau = scenario.generate_tau;
end
if isfield(scenario,"fit_tau")
    fit_tau = scenario.fit_tau;
end
cr_target = scenario.censor_rate;
data_copula_family = char(scenario.generate_copula);
fit_copula_family = char(scenario.fit_copula);
methods = string(scenario.methods);
data_only = any(methods == "DataOnly");
estimate_tau = isfield(scenario, "profile_tau") && logical(scenario.profile_tau);
additive_linear_coefficient = scenario.additive_linear_coefficient;
linear_signal_scale = scenario.linear_signal_scale;
additive_signal_scale = scenario.additive_signal_scale;
dgp_type = lower(string(scenario.dgp_type));
censoring_function_design = lower(string(scenario.censoring_function_design));

n0 = choose_bin_count(n,opts.BaselineBinCount);
block_size = 10;
q_size = 5;
a = opts.MICSharpness;
mu = zeros(p,1);
c = (1:p);
ama = bsxfun(@minus,c,c');
sigma = rho.^(abs(ama));
if dgp_type == "discriminating"
    if p < 7 || q < 6
        error('PLAPH:DiscriminatingDGPDimension', ...
            'The discriminating DGP requires p>=7 and q>=6.');
    end
    Beta = zeros(p,1);
    PHI = zeros(p,1);
    Beta([1 3 6]) = linear_signal_scale.*[1.0;-0.8;0.6];
    PHI([2 4 7]) = linear_signal_scale.*[-0.8;1.0;0.6];
    active_additive_T = [1 2 3 4];
    if any(censoring_function_design == ...
            ["shared_opposite" "shared_opposite_bounded" ...
            "shared_balanced_bounded"])
        active_additive_C = [1 2 3 4];
    else
        active_additive_C = [1 3 5 6];
    end
    lamt = scenario.failure_scale;
    censoring_scale = scenario.censoring_scale;
else
    Beta = linear_signal_scale.*[1.0;0;-0.8;0;0;0.6;zeros(p-6,1)];
    PHI = linear_signal_scale.*[0;1.0;0;0;0;-0.8;zeros(p-6,1)];
    active_additive_T = 1:min(4,q);
    active_additive_C = active_additive_T;
    lamt = lambda_guess(cr_target, generate_tau);
    censoring_scale = 5;
end
active_beta = find(Beta ~= 0);
active_phi = find(PHI ~= 0);
inactive_beta_reference = find(Beta == 0,1,"first");
active_additive_union = union(active_additive_T,active_additive_C,'stable');

num_internal_knots = opts.NumInternalKnots;
threshold_linear_beta = opts.ThresholdLinearFailure;
threshold_linear_phi = opts.ThresholdLinearCensoring;
threshold_spline = opts.ThresholdSpline;
max_outer_iter = opts.MaxOuterIter;
mic_penalty_scale = opts.MICPenaltyScale;
mic_penalty_scale_T = opts.MICPenaltyScaleFailure;
mic_penalty_scale_C = opts.MICPenaltyScaleCensoring;
mic_penalty_margin_specific = opts.MICPenaltyMarginSpecific;
mic_penalty_normalization = opts.MICPenaltyNormalization;
linear_mic_form = opts.LinearMICForm;
linear_kink_coefficient = opts.LinearKinkCoefficient;
linear_penalty_tuning_criterion = opts.LinearPenaltyTuningCriterion;
linear_penalty_multiplier_grid = opts.LinearPenaltyMultiplierGrid;
linear_penalty_ebic_gamma = opts.LinearPenaltyEBICGamma;
linear_l1_type = opts.LinearL1Type;
linear_l1_scale = opts.LinearL1Scale;
linear_ridge_scale = opts.LinearRidgeScale;
adaptive_weight_gamma = opts.AdaptiveWeightGamma;
adaptive_weight_epsilon = opts.AdaptiveWeightEpsilon;
adaptive_weight_cap = opts.AdaptiveWeightCap;
initial_linear_strategy = opts.InitialLinearStrategy;
initial_linear_pilot_scale = opts.InitialLinearPilotScale;
initial_linear_cap = opts.InitialLinearCap;
initial_spline_weighting = opts.InitialSplineWeighting;
initial_spline_strategy = opts.InitialSplineStrategy;
initial_spline_pilot_scale = opts.InitialSplinePilotScale;
mpl_step_scale = opts.MPLStepScale;
mpl_step_options = mpl_step_options_from_opts(opts);
decomposable_pilot_options = decomposable_pilot_options_from_opts(opts);
collect_diagnostics = opts.CollectDiagnostics;
spline_basis_scaling = opts.SplineBasisScaling;
spline_threshold_stage = opts.SplineThresholdStage;
spline_selection_metric = opts.SplineSelectionMetric;
spline_group_penalty_scale = opts.SplineGroupPenaltyScale;
spline_group_penalty_scale_T = opts.SplineGroupPenaltyScaleFailure;
spline_group_penalty_scale_C = opts.SplineGroupPenaltyScaleCensoring;
spline_group_penalty_type = opts.SplineGroupPenaltyType;
spline_group_penalty_normalization = opts.SplineGroupPenaltyNormalization;
spline_group_mic_sharpness = opts.SplineGroupMICSharpness;
spline_group_mic_continuation = opts.SplineGroupMICContinuation;
spline_group_mic_sharpness_path = opts.SplineGroupMICSharpnessPath;
spline_group_mic_form = opts.SplineGroupMICForm;
spline_group_kink_coefficient = opts.SplineGroupKinkCoefficient;
spline_group_mcp_concavity = opts.SplineGroupMCPConcavity;
spline_group_iht_threshold_T = opts.SplineGroupIHTThresholdFailure;
spline_group_iht_threshold_C = opts.SplineGroupIHTThresholdCensoring;
spline_group_iht_warmup_iterations = opts.SplineGroupIHTWarmupIterations;
spline_group_iht_continuation_iterations = opts.SplineGroupIHTContinuationIterations;
spline_group_iht_start_fraction = opts.SplineGroupIHTStartFraction;
post_selection_refit_type = opts.PostSelectionRefitType;
collect_function_fits = opts.CollectFunctionFits;
pre_screening_type = opts.PreScreeningType;

rep = repmat(empty_replication_result(), R, 1);
diagnostics = repmat(empty_diagnostic_result(), R, 1);
screening = cell(R, 1);
linear_penalty_tuning = cell(R,1);
selected_supports = cell(R,1);
if collect_function_fits
    function_fits = initialize_function_fit_storage( ...
        opts.FunctionFitGridSize, active_additive_union, R, ...
        additive_linear_coefficient, additive_signal_scale,dgp_type, ...
        censoring_function_design);
end

for iter = 1:R
    seed = opts.BaseSeed + iter;
    % Process workers default to Threefry, whereas the historical serial
    % campaign uses Twister. Pin the generator so a seed identifies the same
    % Monte Carlo sample independently of worker scheduling.
    rng(seed,"twister");
    fprintf('  Replication %d/%d, seed=%d\n', iter, R, seed);
    rep(iter).replication = iter;
    rep(iter).seed = seed;
    rep(iter).generate_tau = generate_tau;
    rep(iter).fit_tau = fit_tau;

    try
        [T,Z,X,B,Delta,T1,Z1,X1,B1,Delta1,num_basis] = survival_data( ...
            n,q,Beta,PHI,mu,sigma,data_copula_family,generate_tau,lamt,iter, ...
            num_internal_knots,spline_basis_scaling,additive_linear_coefficient, ...
            additive_signal_scale,dgp_type,rho,censoring_scale, ...
            censoring_function_design,~data_only);
        [rep(iter).data_time_signature,rep(iter).data_linear_signature, ...
            rep(iter).data_additive_signature,rep(iter).data_basis_signature, ...
            rep(iter).data_event_signature] = paired_data_signatures( ...
            T,Z,X,B,Delta);
        if strlength(string(opts.ExportDataDir)) > 0
            % v2: export the generated sample for external comparators
            % (grpreg::grpsurv). Columns: time, status, Z_1..Z_p, X_1..X_q.
            export_dir = char(opts.ExportDataDir);
            if ~isfolder(export_dir)
                mkdir(export_dir);
            end
            export_file = fullfile(export_dir,sprintf('data_seed%d.csv',seed));
            writematrix([T(:),double(Delta(:)),Z,X],export_file);
            gzip(export_file);
            delete(export_file);
        end
        rep(iter).censor_rate = 1 - mean(Delta);
        rep(iter).n_failure_events = sum(Delta);
        rep(iter).n_censoring_events = sum(1-Delta);
        x2_distance = abs(X(:,2)-median(X(:,2)));
        x2_center = x2_distance <= quantile(x2_distance,0.50);
        x2_tail = x2_distance >= quantile(x2_distance,0.80);
        rep(iter).censor_rate_x2_center = mean(1-Delta(x2_center));
        rep(iter).censor_rate_x2_tail = mean(1-Delta(x2_tail));
        rep(iter).x2_tail_to_center_censor_ratio = ...
            rep(iter).censor_rate_x2_tail ./ ...
            max(rep(iter).censor_rate_x2_center,eps);
        if data_only
            rep(iter).ok = true;
            rep(iter).message = "";
            continue;
        end
        if collect_diagnostics
            diagnostics(iter) = data_diagnostics(iter, seed, n, p, q, Z, X, B, Delta, T, num_basis);
        end

        if pre_screening_type == "oracle_support_diagnostic"
            if dgp_type ~= "discriminating"
                error('PLAPH:OracleScreenDiagnosticDGP', ...
                    'Oracle-support screening is restricted to the discriminating diagnostic DGP.');
            end
            screen = struct( ...
                'type', pre_screening_type, ...
                'linear_T', active_beta(:), 'linear_C', active_phi(:), ...
                'linear_union', union(active_beta(:),active_phi(:),'stable'), ...
                'additive_T', active_additive_T(:), ...
                'additive_C', active_additive_C(:), ...
                'additive_union', union(active_additive_T(:), ...
                active_additive_C(:),'stable'), ...
                'linear_score_T', NaN(p,1), 'linear_score_C', NaN(p,1), ...
                'additive_score_T', NaN(q,1), 'additive_score_C', NaN(q,1));
        else
            screen = build_pre_screening(pre_screening_type, T, Z, X, B, Delta, ...
                num_basis, opts.ScreenLinearSize, opts.ScreenAdditiveSize, ...
                opts.ScreenAdditiveRank);
        end
        screening{iter} = screen;
        rep(iter) = assign_screening_fields(rep(iter),screen,active_beta,active_phi, ...
            active_additive_T,active_additive_C);
        linear_idx = screen.linear_union;
        additive_idx = screen.additive_union;
        linear_mask_T = ismember(linear_idx,screen.linear_T);
        linear_mask_C = ismember(linear_idx,screen.linear_C);
        additive_mask_T = ismember(additive_idx,screen.additive_T);
        additive_mask_C = ismember(additive_idx,screen.additive_C);
        spline_cols = spline_group_columns(additive_idx, num_basis);
        q_fit = numel(additive_idx);
        Z_fit = Z(:, linear_idx);
        Z1_fit = Z1(:, linear_idx);
        X_fit = X(:, additive_idx);
        X1_fit = X1(:, additive_idx);
        B_fit = B(:, spline_cols);
        B1_fit = B1(:, spline_cols);

        [ini_beta,ini_phi,ini_spline_beta,ini_spline_phi] = initial_value( ...
            n,Z_fit,X_fit,B_fit,1e-5,Delta,a,block_size,q_size, ...
            mic_penalty_scale,mic_penalty_normalization,initial_linear_strategy, ...
            initial_linear_cap,initial_spline_weighting,initial_linear_pilot_scale, ...
            initial_spline_strategy,initial_spline_pilot_scale, ...
            mic_penalty_scale_T,mic_penalty_scale_C,mic_penalty_margin_specific);
        ini_beta(~linear_mask_T) = 0;
        ini_phi(~linear_mask_C) = 0;
        ini_spline_beta(~repelem(additive_mask_T(:),num_basis)) = 0;
        ini_spline_phi(~repelem(additive_mask_C(:),num_basis)) = 0;
        if spline_group_penalty_type == "group_iht" && ...
                spline_group_iht_warmup_iterations == 0 && ...
                decomposable_pilot_options.iterations == 0
            ini_spline_beta(:) = 0;
            ini_spline_phi(:) = 0;
        end
        ini_beta_full = expand_linear_coefficients(ini_beta, linear_idx, p);
        ini_phi_full = expand_linear_coefficients(ini_phi, linear_idx, p);
        ini_spline_beta_full = expand_spline_coefficients(ini_spline_beta, additive_idx, num_basis, q);
        ini_spline_phi_full = expand_spline_coefficients(ini_spline_phi, additive_idx, num_basis, q);
        rep(iter).ini_beta_selected = count_selected_linear(ini_beta_full, threshold_linear_beta);
        rep(iter).ini_phi_selected = count_selected_linear(ini_phi_full, threshold_linear_phi);
        rep(iter).ini_beta_min_signal = min(abs(ini_beta_full(active_beta)));
        rep(iter).ini_phi_min_signal = min(abs(ini_phi_full(active_phi)));
        rep(iter).ini_beta_max_noise = max_abs_noise(ini_beta_full, active_beta);
        rep(iter).ini_phi_max_noise = max_abs_noise(ini_phi_full, active_phi);
        if collect_diagnostics
            diagnostics(iter).ini_beta_selected = rep(iter).ini_beta_selected;
            diagnostics(iter).ini_phi_selected = rep(iter).ini_phi_selected;
            diagnostics(iter).ini_beta_min_signal = rep(iter).ini_beta_min_signal;
            diagnostics(iter).ini_phi_min_signal = rep(iter).ini_phi_min_signal;
            diagnostics(iter).ini_beta_max_noise = rep(iter).ini_beta_max_noise;
            diagnostics(iter).ini_phi_max_noise = rep(iter).ini_phi_max_noise;
            diagnostics(iter).ini_spline_beta_max_abs = max(abs(ini_spline_beta_full));
            diagnostics(iter).ini_spline_phi_max_abs = max(abs(ini_spline_phi_full));
            ini_effect_T = active_effect_norms(B, ini_spline_beta_full, ...
                num_basis, q, active_additive_T);
            ini_effect_C = active_effect_norms(B, ini_spline_phi_full, ...
                num_basis, q, active_additive_C);
            for jj = 1:min(4, numel(ini_effect_T))
                diagnostics(iter).(sprintf('ini_spline_beta_eff_g%d', jj)) = ini_effect_T(jj);
                diagnostics(iter).(sprintf('ini_spline_phi_eff_g%d', jj)) = ini_effect_C(jj);
            end
        end

        if any(methods == "Oracle")
            oracle_spline_cols = spline_group_columns(active_additive_T,num_basis);
            [oracle_beta, oracle_bspline_beta] = oracle_value( ...
                n,Z(:,active_beta),X(:,active_additive_T),B(:,oracle_spline_cols),1e-5,Delta);
            for jj = 1:min(3,numel(oracle_beta))
                rep(iter).(sprintf('oracle_beta_active%d',jj)) = oracle_beta(jj);
            end
            rep(iter).mse_oracle = (oracle_beta-Beta(active_beta))' * ...
                sigma(active_beta,active_beta) * (oracle_beta-Beta(active_beta));
            rep(iter).oracle_selected_additive = count_selected_groups_metric( ...
                B(:,oracle_spline_cols),oracle_bspline_beta,num_basis, ...
                numel(active_additive_T),threshold_spline,spline_selection_metric);
            rep(iter).imse_oracle_add = additive_imse(B, X, oracle_bspline_beta, ...
                num_basis,q,active_additive_T,additive_linear_coefficient, ...
                additive_signal_scale,dgp_type,"T",true, ...
                censoring_function_design, ...
                opts.FunctionFitGridSize);
            if collect_function_fits
                oracle_coef_full = zeros(q*num_basis, 1);
                for jj = 1:numel(active_additive_T)
                    full_idx = ((active_additive_T(jj)-1)*num_basis+1):(active_additive_T(jj)*num_basis);
                    oracle_idx = ((jj-1)*num_basis+1):(jj*num_basis);
                    oracle_coef_full(full_idx) = oracle_bspline_beta(oracle_idx);
                end
                function_fits = store_function_fit_curves( ...
                    function_fits,"oracle",iter,B,X,oracle_coef_full,num_basis, ...
                    active_additive_union);
            end
        end

        if any(methods == "MPL_fixed_tau") || any(methods == "MPL_profile_tau")
            tune_linear_penalty = any(linear_mic_form == ...
                ["kinked" "scad" "mcp"]) && ...
                linear_penalty_tuning_criterion ~= "none";
            tuning_mode = lower(string(opts.TuningMode));
            % v2 candidate grid: columns are [xi, c_L(T), c_L(C), kappa_T, kappa_C].
            if tuning_mode == "fixed"
                candidate_multipliers = linear_kink_coefficient;
                if tune_linear_penalty
                    candidate_multipliers = sort( ...
                        unique(linear_penalty_multiplier_grid(:)),'descend');
                end
                candidate_grid = [candidate_multipliers(:), ...
                    repmat([mic_penalty_scale_T mic_penalty_scale_C ...
                    spline_group_iht_threshold_T spline_group_iht_threshold_C], ...
                    numel(candidate_multipliers),1)];
            else
                tune_linear_penalty = true;
                linear_penalty_tuning_criterion = erase(tuning_mode,"grid_");
                [grid_cL,grid_kappa] = ndgrid(opts.TuningGridCL(:), ...
                    opts.TuningGridKappa(:));
                candidate_grid = [repmat(linear_kink_coefficient,numel(grid_cL),1), ...
                    grid_cL(:),grid_cL(:),grid_kappa(:),grid_kappa(:)];
            end
            candidate_multipliers = candidate_grid(:,1);
            candidate_rows = cell(size(candidate_grid,1),1);
            candidate_store = cell(size(candidate_grid,1),1);
            best_value = Inf;
            best_index = NaN;
            best_fit = struct();
            warm_beta = ini_beta;
            warm_phi = ini_phi;
            warm_spline_beta = ini_spline_beta;
            warm_spline_phi = ini_spline_phi;
            for candidate_index = 1:size(candidate_grid,1)
                candidate_multiplier = candidate_grid(candidate_index,1);
                candidate_cL_T = candidate_grid(candidate_index,2);
                candidate_cL_C = candidate_grid(candidate_index,3);
                candidate_kappa_T = candidate_grid(candidate_index,4);
                candidate_kappa_C = candidate_grid(candidate_index,5);
                if tuning_mode ~= "fixed"
                    % Cold start at every grid point: the estimator does not
                    % depend on the order in which the grid is visited.
                    warm_beta = ini_beta;
                    warm_phi = ini_phi;
                    warm_spline_beta = ini_spline_beta;
                    warm_spline_phi = ini_spline_phi;
                end
                [candidate_beta,candidate_phi,candidate_spline_beta, ...
                    candidate_spline_phi,~,~,~,~,candidate_tau, ...
                    candidate_alpha,candidate_diag,candidate_pilot, ...
                    candidate_inference] = MPL( ...
                    n,q_fit,n0,warm_beta,warm_phi,warm_spline_beta, ...
                    warm_spline_phi,Z1_fit,X1_fit,B1_fit,T1,1e-5,Delta1, ...
                    a,fit_copula_family,fit_tau,block_size,q_size, ...
                    threshold_linear_beta,threshold_linear_phi, ...
                    threshold_spline,max_outer_iter,mic_penalty_scale, ...
                    mic_penalty_normalization,mpl_step_scale, ...
                    collect_diagnostics,spline_threshold_stage, ...
                    spline_selection_metric,spline_group_penalty_scale, ...
                    spline_group_penalty_type, ...
                    spline_group_penalty_normalization, ...
                    spline_group_mic_sharpness, ...
                    spline_group_mic_continuation, ...
                    spline_group_mic_sharpness_path,estimate_tau, ...
                    linear_mic_form,candidate_multiplier, ...
                    linear_l1_type,linear_l1_scale,linear_ridge_scale, ...
                    adaptive_weight_gamma,adaptive_weight_epsilon, ...
                    adaptive_weight_cap,spline_group_mic_form, ...
                    spline_group_kink_coefficient, ...
                    spline_group_penalty_scale_T, ...
                    spline_group_penalty_scale_C, ...
                    spline_group_mcp_concavity, ...
                    candidate_kappa_T, ...
                    candidate_kappa_C, ...
                    spline_group_iht_warmup_iterations, ...
                    post_selection_refit_type, ...
                    spline_group_iht_continuation_iterations, ...
                    spline_group_iht_start_fraction,mpl_step_options, ...
                    decomposable_pilot_options,linear_mask_T,linear_mask_C, ...
                    additive_mask_T,additive_mask_C, ...
                    opts.PostSelectionRefitMaxIter, ...
                    opts.BaselineRoughnessScale,candidate_cL_T, ...
                    candidate_cL_C,mic_penalty_margin_specific);

                [candidate_bic,candidate_ebic,candidate_df, ...
                    selected_linear_T,selected_linear_C, ...
                    selected_additive_T,selected_additive_C] = ...
                    linear_penalty_information_criteria( ...
                    candidate_diag.mpl_final_total_loglik,n,n0,num_basis, ...
                    candidate_pilot,linear_mask_T,linear_mask_C, ...
                    additive_mask_T,additive_mask_C, ...
                    linear_penalty_ebic_gamma,Z_fit,B_fit);
                candidate_valid = isfinite( ...
                    candidate_diag.mpl_final_total_loglik) && ...
                    ~candidate_diag.mpl_armijo_rejected && ...
                    candidate_diag.mpl_selection_gate_pass;
                if any(post_selection_refit_type == ...
                        ["joint_oracle" "joint_oracle_alpha"])
                    candidate_valid = candidate_valid && ...
                        candidate_diag.mpl_oracle_refit_stationarity_pass;
                end
                if linear_penalty_tuning_criterion == "bic"
                    candidate_value = candidate_bic;
                elseif linear_penalty_tuning_criterion == "ebic"
                    candidate_value = candidate_ebic;
                else
                    candidate_value = 0;
                end
                candidate_rows{candidate_index} = table( ...
                    candidate_index,candidate_multiplier, ...
                    candidate_cL_T,candidate_kappa_T, ...
                    candidate_diag.mpl_linear_mic_lambda_T*candidate_multiplier, ...
                    candidate_diag.mpl_linear_mic_lambda_C*candidate_multiplier, ...
                    candidate_diag.mpl_final_total_loglik,candidate_df, ...
                    candidate_bic,candidate_ebic,candidate_value, ...
                    selected_linear_T,selected_linear_C,selected_additive_T, ...
                    selected_additive_C, ...
                    candidate_diag.mpl_armijo_rejected, ...
                    candidate_diag.mpl_oracle_refit_stationarity_pass, ...
                    candidate_diag.mpl_oracle_refit_final_stationarity_inf, ...
                    candidate_diag.mpl_inference_valid,candidate_valid,false, ...
                    'VariableNames',{'candidate_index','multiplier', ...
                    'c_L','kappa', ...
                    'lambda_T','lambda_C','total_loglik','degrees_of_freedom', ...
                    'bic','ebic','criterion_value','selected_linear_T', ...
                    'selected_linear_C','selected_additive_T', ...
                    'selected_additive_C','armijo_rejected', ...
                    'stationarity_pass','stationarity_inf','inference_valid', ...
                    'numerically_valid','selected'});
                candidate_store{candidate_index} = struct('beta',candidate_beta, ...
                    'phi',candidate_phi,'spline_beta',candidate_spline_beta, ...
                    'spline_phi',candidate_spline_phi,'tau',candidate_tau, ...
                    'alpha',candidate_alpha,'diag',candidate_diag, ...
                    'pilot',candidate_pilot,'inference',candidate_inference, ...
                    'bic',candidate_bic,'ebic',candidate_ebic, ...
                    'df',candidate_df,'multiplier',candidate_multiplier, ...
                    'cL',candidate_cL_T,'kappa',candidate_kappa_T, ...
                    'valid',candidate_valid); %#ok<AGROW>
                selection_eligible = candidate_valid || ~tune_linear_penalty;
                if selection_eligible && candidate_value < best_value
                    best_value = candidate_value;
                    best_index = candidate_index;
                    best_fit = struct('beta',candidate_beta, ...
                        'phi',candidate_phi,'spline_beta',candidate_spline_beta, ...
                        'spline_phi',candidate_spline_phi,'tau',candidate_tau, ...
                        'alpha',candidate_alpha,'diag',candidate_diag, ...
                        'pilot',candidate_pilot,'inference',candidate_inference, ...
                        'bic',candidate_bic,'ebic',candidate_ebic, ...
                        'df',candidate_df,'multiplier',candidate_multiplier, ...
                        'cL',candidate_cL_T,'kappa',candidate_kappa_T, ...
                        'valid',candidate_valid);
                end
                if candidate_valid && candidate_index < numel(candidate_multipliers)
                    warm_beta = candidate_beta;
                    warm_phi = candidate_phi;
                    warm_spline_beta = candidate_spline_beta;
                    warm_spline_phi = candidate_spline_phi;
                end
            end
            tuning_path = vertcat(candidate_rows{:});
            tuning_path.criterion = repmat( ...
                linear_penalty_tuning_criterion,height(tuning_path),1);
            linear_penalty_tuning{iter} = tuning_path;
            if isnan(best_index)
                % v2: no candidate passed both gates. Select by the criterion
                % among all finite candidates and record the fit as not accepted.
                criterion_values = tuning_path.criterion_value;
                criterion_values(~isfinite(criterion_values)) = Inf;
                [best_value,best_index] = min(criterion_values);
                if ~isfinite(best_value)
                    error('PLAPH:LinearPenaltyTuningNoValidCandidate', ...
                        'No finite %s candidate was available.', linear_mic_form);
                end
                warning('PLAPH:NoAcceptedCandidate', ...
                    'No tuning candidate passed both numerical gates; refitting the best finite one.');
                best_fit = candidate_store{best_index};
            end
            tuning_tie_count = NaN;
            tuning_on_boundary = false;
            if tuning_mode ~= "fixed" && any(tuning_path.numerically_valid)
                % v2: EBIC ties (in practice the same selected model at several
                % grid points) are broken deterministically in favour of the
                % most regularized point, so that the reported tuning does not
                % depend on refit rounding.
                [best_index,tuning_tie_count,tuning_on_boundary] = ...
                    plaph_tuning_tie_break(tuning_path.criterion_value, ...
                    tuning_path.numerically_valid,tuning_path.c_L, ...
                    tuning_path.kappa,opts.TuningGridCL,opts.TuningGridKappa);
                best_fit = candidate_store{best_index};
            end
            tuning_path.selected(best_index) = true;
            linear_penalty_tuning{iter} = tuning_path;
            mpl_beta = best_fit.beta;
            mpl_phi = best_fit.phi;
            bspline_beta = best_fit.spline_beta;
            bspline_phi = best_fit.spline_phi;
            tau_hat = best_fit.tau;
            alph_hat = best_fit.alpha;
            mpl_diag = best_fit.diag;
            pilot_state = best_fit.pilot;
            mpl_inference = best_fit.inference;
            % v2 acceptance and tuning bookkeeping.
            rep(iter).spline_order = plaph_spline_order();
            rep(iter).mpl_selection_gate_pass = mpl_diag.mpl_selection_gate_pass;
            rep(iter).mpl_selection_status_code = mpl_diag.mpl_selection_status_code;
            rep(iter).mpl_selection_gate_maximum = mpl_diag.mpl_selection_gate_maximum;
            rep(iter).mpl_refit_stationarity_pass = ...
                mpl_diag.mpl_oracle_refit_stationarity_pass;
            rep(iter).mpl_refit_stationarity_inf = ...
                mpl_diag.mpl_oracle_refit_final_stationarity_inf;
            rep(iter).mpl_accepted = rep(iter).mpl_selection_gate_pass && ...
                rep(iter).mpl_refit_stationarity_pass && ~mpl_diag.mpl_armijo_rejected;
            rep(iter).mpl_step_doublings = mpl_diag.mpl_step_doublings;
            rep(iter).mpl_target_iterations = mpl_diag.mpl_target_iterations;
            rep(iter).mpl_df_rank = best_fit.df;
            rep(iter).tuning_cL = best_fit.cL;
            rep(iter).tuning_kappa = best_fit.kappa;
            rep(iter).tuning_ebic = best_fit.ebic;
            rep(iter).tuning_valid_count = sum(tuning_path.numerically_valid);
            if tuning_mode ~= "fixed"
                % True when no EBIC-optimal grid point (ties included) is
                % interior to the grid.
                rep(iter).tuning_at_boundary = tuning_on_boundary;
                rep(iter).tuning_tie_count = tuning_tie_count;
            end
            final_mask_beta = expand_support_mask( ...
                pilot_state.final_selected_linear_T_mask,linear_idx,p);
            final_mask_phi = expand_support_mask( ...
                pilot_state.final_selected_linear_C_mask,linear_idx,p);
            final_mask_add_T = expand_support_mask( ...
                pilot_state.final_selected_additive_T_mask,additive_idx,q);
            final_mask_add_C = expand_support_mask( ...
                pilot_state.final_selected_additive_C_mask,additive_idx,q);
            if pilot_state.enabled && pilot_state.completed
                rep(iter).mpl_pilot_final_agree_linear_T = double(isequal( ...
                    abs(pilot_state.beta(:)) > 0, ...
                    pilot_state.final_selected_linear_T_mask(:)));
                rep(iter).mpl_pilot_final_agree_linear_C = double(isequal( ...
                    abs(pilot_state.phi(:)) > 0, ...
                    pilot_state.final_selected_linear_C_mask(:)));
                rep(iter).mpl_pilot_final_agree_add_T = double(isequal( ...
                    spline_group_support(pilot_state.bspline_beta,num_basis,q_fit), ...
                    pilot_state.final_selected_additive_T_mask(:)));
                rep(iter).mpl_pilot_final_agree_add_C = double(isequal( ...
                    spline_group_support(pilot_state.bspline_phi,num_basis,q_fit), ...
                    pilot_state.final_selected_additive_C_mask(:)));
            end
            rep(iter).linear_penalty_tuned = tune_linear_penalty;
            if tune_linear_penalty
                rep(iter).linear_penalty_tuning_criterion = ...
                    linear_penalty_tuning_criterion;
            end
            rep(iter).linear_penalty_candidate_count = height(tuning_path);
            rep(iter).linear_penalty_valid_candidate_count = ...
                sum(tuning_path.numerically_valid);
            rep(iter).linear_penalty_selected_multiplier = best_fit.multiplier;
            rep(iter).linear_penalty_selected_lambda_T = ...
                mpl_diag.mpl_linear_mic_lambda_T*best_fit.multiplier;
            rep(iter).linear_penalty_selected_lambda_C = ...
                mpl_diag.mpl_linear_mic_lambda_C*best_fit.multiplier;
            rep(iter).linear_penalty_selected_criterion = best_value;
            rep(iter).linear_penalty_selected_bic = best_fit.bic;
            rep(iter).linear_penalty_selected_ebic = best_fit.ebic;
            rep(iter).linear_penalty_selected_df = best_fit.df;
            rep(iter).linear_penalty_selected_at_grid_boundary = ...
                tune_linear_penalty && ...
                (best_index == 1 || best_index == height(tuning_path));
            if pilot_state.enabled
                pilot_beta = expand_linear_coefficients(pilot_state.beta, linear_idx, p);
                pilot_phi = expand_linear_coefficients(pilot_state.phi, linear_idx, p);
                pilot_spline_beta = expand_spline_coefficients( ...
                    pilot_state.bspline_beta, additive_idx, num_basis, q);
                pilot_spline_phi = expand_spline_coefficients( ...
                    pilot_state.bspline_phi, additive_idx, num_basis, q);
                rep(iter).mpl_pilot_completed = pilot_state.completed;
                rep(iter).mpl_pilot_iterations = pilot_state.iterations;
                rep(iter).mpl_pilot_beta_fp = count_fp( ...
                    pilot_beta,active_beta,threshold_linear_beta);
                rep(iter).mpl_pilot_beta_fn = count_fn( ...
                    pilot_beta,active_beta,threshold_linear_beta);
                rep(iter).mpl_pilot_phi_fp = count_fp( ...
                    pilot_phi,active_phi,threshold_linear_phi);
                rep(iter).mpl_pilot_phi_fn = count_fn( ...
                    pilot_phi,active_phi,threshold_linear_phi);
                rep(iter).mpl_pilot_beta_min_signal = min(abs(pilot_beta(active_beta)));
                rep(iter).mpl_pilot_phi_min_signal = min(abs(pilot_phi(active_phi)));
                rep(iter).mpl_pilot_beta_max_noise = max_abs_noise(pilot_beta,active_beta);
                rep(iter).mpl_pilot_phi_max_noise = max_abs_noise(pilot_phi,active_phi);
                rep(iter).mpl_pilot_add_T_fp = count_group_fp_metric( ...
                    B1,pilot_spline_beta,num_basis,q,active_additive_T,0,"effect");
                rep(iter).mpl_pilot_add_T_fn = count_group_fn_metric( ...
                    B1,pilot_spline_beta,num_basis,q,active_additive_T,0,"effect");
                rep(iter).mpl_pilot_add_C_fp = count_group_fp_metric( ...
                    B1,pilot_spline_phi,num_basis,q,active_additive_C,0,"effect");
                rep(iter).mpl_pilot_add_C_fn = count_group_fn_metric( ...
                    B1,pilot_spline_phi,num_basis,q,active_additive_C,0,"effect");
                rep(iter).mpl_pilot_eff_T_active_min = min_active_effect_norm( ...
                    B1,pilot_spline_beta,num_basis,q,active_additive_T);
                rep(iter).mpl_pilot_eff_T_noise_max = max_noise_effect_norm( ...
                    B1,pilot_spline_beta,num_basis,q,active_additive_T);
                rep(iter).mpl_pilot_eff_C_active_min = min_active_effect_norm( ...
                    B1,pilot_spline_phi,num_basis,q,active_additive_C);
                rep(iter).mpl_pilot_eff_C_noise_max = max_noise_effect_norm( ...
                    B1,pilot_spline_phi,num_basis,q,active_additive_C);
                rep(iter).mpl_pilot_stationarity_beta = pilot_state.stationarity_beta;
                rep(iter).mpl_pilot_stationarity_phi = pilot_state.stationarity_phi;
                rep(iter).mpl_pilot_stationarity_group_T = pilot_state.stationarity_group_T;
                rep(iter).mpl_pilot_stationarity_group_C = pilot_state.stationarity_group_C;

                target_beta_candidate = expand_linear_coefficients( ...
                    pilot_state.first_target_beta_candidate,linear_idx,p);
                target_beta_threshold = expand_linear_thresholds( ...
                    pilot_state.first_target_beta_threshold,linear_idx,p);
                target_phi_candidate = expand_linear_coefficients( ...
                    pilot_state.first_target_phi_candidate,linear_idx,p);
                target_phi_threshold = expand_linear_thresholds( ...
                    pilot_state.first_target_phi_threshold,linear_idx,p);
                beta_margin = abs(target_beta_candidate)-target_beta_threshold;
                phi_margin = abs(target_phi_candidate)-target_phi_threshold;
                rep(iter).mpl_target_beta_candidate_active_min = ...
                    min(abs(target_beta_candidate(active_beta)));
                rep(iter).mpl_target_beta_candidate_noise_max = ...
                    max_abs_noise(target_beta_candidate,active_beta);
                rep(iter).mpl_target_beta_active_margin_min = min(beta_margin(active_beta));
                rep(iter).mpl_target_beta_noise_margin_max = ...
                    max_noise_value(beta_margin,active_beta);
                rep(iter).mpl_target_beta_strict_separation = double( ...
                    rep(iter).mpl_target_beta_active_margin_min > 0 && ...
                    rep(iter).mpl_target_beta_noise_margin_max <= 0);
                rep(iter).mpl_target_beta_threshold_min = ...
                    finite_min(target_beta_threshold);
                rep(iter).mpl_target_beta_threshold_max = ...
                    finite_max(target_beta_threshold);
                rep(iter).mpl_target_phi_candidate_active_min = ...
                    min(abs(target_phi_candidate(active_phi)));
                rep(iter).mpl_target_phi_candidate_noise_max = ...
                    max_abs_noise(target_phi_candidate,active_phi);
                rep(iter).mpl_target_phi_active_margin_min = min(phi_margin(active_phi));
                rep(iter).mpl_target_phi_noise_margin_max = ...
                    max_noise_value(phi_margin,active_phi);
                rep(iter).mpl_target_phi_strict_separation = double( ...
                    rep(iter).mpl_target_phi_active_margin_min > 0 && ...
                    rep(iter).mpl_target_phi_noise_margin_max <= 0);
                rep(iter).mpl_target_phi_threshold_min = ...
                    finite_min(target_phi_threshold);
                rep(iter).mpl_target_phi_threshold_max = ...
                    finite_max(target_phi_threshold);

                target_group_T_candidate = expand_group_values( ...
                    pilot_state.first_target_group_T_candidate_norm,additive_idx,q);
                target_group_C_candidate = expand_group_values( ...
                    pilot_state.first_target_group_C_candidate_norm,additive_idx,q);
                group_T_margin = target_group_T_candidate- ...
                    pilot_state.first_target_group_T_threshold;
                group_C_margin = target_group_C_candidate- ...
                    pilot_state.first_target_group_C_threshold;
                rep(iter).mpl_target_group_T_candidate_active_min = ...
                    min(target_group_T_candidate(active_additive_T));
                rep(iter).mpl_target_group_T_candidate_noise_max = ...
                    max_noise_value(target_group_T_candidate,active_additive_T);
                rep(iter).mpl_target_group_T_active_margin_min = ...
                    min(group_T_margin(active_additive_T));
                rep(iter).mpl_target_group_T_noise_margin_max = ...
                    max_noise_value(group_T_margin,active_additive_T);
                rep(iter).mpl_target_group_T_threshold = ...
                    pilot_state.first_target_group_T_threshold;
                rep(iter).mpl_target_group_T_strict_separation = double( ...
                    rep(iter).mpl_target_group_T_active_margin_min > 0 && ...
                    rep(iter).mpl_target_group_T_noise_margin_max <= 0);
                rep(iter).mpl_target_group_C_candidate_active_min = ...
                    min(target_group_C_candidate(active_additive_C));
                rep(iter).mpl_target_group_C_candidate_noise_max = ...
                    max_noise_value(target_group_C_candidate,active_additive_C);
                rep(iter).mpl_target_group_C_active_margin_min = ...
                    min(group_C_margin(active_additive_C));
                rep(iter).mpl_target_group_C_noise_margin_max = ...
                    max_noise_value(group_C_margin,active_additive_C);
                rep(iter).mpl_target_group_C_threshold = ...
                    pilot_state.first_target_group_C_threshold;
                rep(iter).mpl_target_group_C_strict_separation = double( ...
                    rep(iter).mpl_target_group_C_active_margin_min > 0 && ...
                    rep(iter).mpl_target_group_C_noise_margin_max <= 0);
            end
            mpl_beta_se = expand_linear_standard_errors( ...
                mpl_inference.linear_beta_se,linear_idx,p);
            mpl_phi_se = expand_linear_standard_errors( ...
                mpl_inference.linear_phi_se,linear_idx,p);
            mpl_beta = expand_linear_coefficients(mpl_beta, linear_idx, p);
            mpl_phi = expand_linear_coefficients(mpl_phi, linear_idx, p);
            bspline_beta = expand_spline_coefficients(bspline_beta, additive_idx, num_basis, q);
            bspline_phi = expand_spline_coefficients(bspline_phi, additive_idx, num_basis, q);
            rep(iter).mse_mpl_beta = (mpl_beta-Beta)' * sigma * (mpl_beta-Beta);
            rep(iter).mse_mpl_phi = (mpl_phi-PHI)' * sigma * (mpl_phi-PHI);
            rep(iter) = assign_active_coefficient_fields( ...
                rep(iter),"mpl_beta",mpl_beta,active_beta);
            rep(iter) = assign_active_coefficient_fields( ...
                rep(iter),"mpl_phi",mpl_phi,active_phi);
            if ~isempty(inactive_beta_reference)
                rep(iter).mpl_beta_inactive1 = ...
                    mpl_beta(inactive_beta_reference);
                rep(iter).mpl_beta_inactive1_index = ...
                    inactive_beta_reference;
            end
            rep(iter) = assign_active_inference_fields( ...
                rep(iter),"mpl_beta",mpl_beta_se,mpl_beta,Beta,active_beta, ...
                mpl_inference.valid);
            rep(iter) = assign_active_inference_fields( ...
                rep(iter),"mpl_phi",mpl_phi_se,mpl_phi,PHI,active_phi, ...
                mpl_inference.valid);
            rep(iter).tau_hat = tau_hat;
            rep(iter).alpha_hat = alph_hat;
            rep(iter).mpl_baseline_armijo_failed = mpl_diag.mpl_baseline_armijo_failed;
            rep(iter).mpl_armijo_rejected = mpl_diag.mpl_armijo_rejected;
            % v2: supports are the frozen pre-refit masks (no second threshold).
            rep(iter).beta_fp = support_false_positives(final_mask_beta,active_beta);
            rep(iter).beta_fn = support_false_negatives(final_mask_beta,active_beta);
            rep(iter).phi_fp = support_false_positives(final_mask_phi,active_phi);
            rep(iter).phi_fn = support_false_negatives(final_mask_phi,active_phi);
            rep(iter).mpl_beta_selected = sum(final_mask_beta);
            rep(iter).mpl_phi_selected = sum(final_mask_phi);
            rep(iter).mpl_beta_min_signal = min(abs(mpl_beta(active_beta)));
            rep(iter).mpl_phi_min_signal = min(abs(mpl_phi(active_phi)));
            rep(iter).mpl_beta_max_noise = max_abs_noise(mpl_beta, active_beta);
            rep(iter).mpl_phi_max_noise = max_abs_noise(mpl_phi, active_phi);
            rep(iter).mpl_beta_noise_q50 = abs_noise_quantile(mpl_beta, active_beta, 0.50);
            rep(iter).mpl_beta_noise_q75 = abs_noise_quantile(mpl_beta, active_beta, 0.75);
            rep(iter).mpl_beta_noise_q90 = abs_noise_quantile(mpl_beta, active_beta, 0.90);
            rep(iter).mpl_beta_noise_q95 = abs_noise_quantile(mpl_beta, active_beta, 0.95);
            rep(iter).mpl_phi_noise_q50 = abs_noise_quantile(mpl_phi, active_phi, 0.50);
            rep(iter).mpl_phi_noise_q75 = abs_noise_quantile(mpl_phi, active_phi, 0.75);
            rep(iter).mpl_phi_noise_q90 = abs_noise_quantile(mpl_phi, active_phi, 0.90);
            rep(iter).mpl_phi_noise_q95 = abs_noise_quantile(mpl_phi, active_phi, 0.95);
            rep(iter).mpl_selected_additive_T = sum(final_mask_add_T);
            rep(iter).mpl_selected_additive_C = sum(final_mask_add_C);
            rep(iter).mpl_add_T_fp = support_false_positives(final_mask_add_T,active_additive_T);
            rep(iter).mpl_add_T_fn = support_false_negatives(final_mask_add_T,active_additive_T);
            rep(iter).mpl_add_C_fp = support_false_positives(final_mask_add_C,active_additive_C);
            rep(iter).mpl_add_C_fn = support_false_negatives(final_mask_add_C,active_additive_C);
            selected_supports{iter} = struct( ...
                'linear_T',find(final_mask_beta), ...
                'linear_C',find(final_mask_phi), ...
                'additive_T',find(final_mask_add_T), ...
                'additive_C',find(final_mask_add_C));
            rep(iter).mpl_add_T_active_min = min_active_group_norm(bspline_beta, num_basis, q, active_additive_T);
            rep(iter).mpl_add_T_noise_max = max_noise_group_norm(bspline_beta, num_basis, q, active_additive_T);
            rep(iter).mpl_add_C_active_min = min_active_group_norm(bspline_phi, num_basis, q, active_additive_C);
            rep(iter).mpl_add_C_noise_max = max_noise_group_norm(bspline_phi, num_basis, q, active_additive_C);
            rep(iter).mpl_eff_T_active_min = min_active_effect_norm(B1, bspline_beta, num_basis, q, active_additive_T);
            rep(iter).mpl_eff_T_noise_max = max_noise_effect_norm(B1, bspline_beta, num_basis, q, active_additive_T);
            rep(iter).mpl_eff_C_active_min = min_active_effect_norm(B1, bspline_phi, num_basis, q, active_additive_C);
            rep(iter).mpl_eff_C_noise_max = max_noise_effect_norm(B1, bspline_phi, num_basis, q, active_additive_C);
            rep(iter).mpl_eff_T_noise_q50 = noise_effect_quantile(B1, bspline_beta, num_basis, q, active_additive_T, 0.50);
            rep(iter).mpl_eff_T_noise_q75 = noise_effect_quantile(B1, bspline_beta, num_basis, q, active_additive_T, 0.75);
            rep(iter).mpl_eff_T_noise_q90 = noise_effect_quantile(B1, bspline_beta, num_basis, q, active_additive_T, 0.90);
            rep(iter).mpl_eff_T_noise_q95 = noise_effect_quantile(B1, bspline_beta, num_basis, q, active_additive_T, 0.95);
            rep(iter).mpl_eff_C_noise_q50 = noise_effect_quantile(B1, bspline_phi, num_basis, q, active_additive_C, 0.50);
            rep(iter).mpl_eff_C_noise_q75 = noise_effect_quantile(B1, bspline_phi, num_basis, q, active_additive_C, 0.75);
            rep(iter).mpl_eff_C_noise_q90 = noise_effect_quantile(B1, bspline_phi, num_basis, q, active_additive_C, 0.90);
            rep(iter).mpl_eff_C_noise_q95 = noise_effect_quantile(B1, bspline_phi, num_basis, q, active_additive_C, 0.95);
            rep(iter) = assign_active_effect_fields(rep(iter), "mpl_eff_T", ...
                active_effect_norms(B1,bspline_beta,num_basis,q,active_additive_T));
            rep(iter) = assign_active_effect_fields(rep(iter), "mpl_eff_C", ...
                active_effect_norms(B1,bspline_phi,num_basis,q,active_additive_C));
            rep(iter).imse_mpl_add_T = additive_imse(B1, X1, bspline_beta, ...
                num_basis,q,active_additive_T,additive_linear_coefficient, ...
                additive_signal_scale,dgp_type,"T",false,censoring_function_design, ...
                opts.FunctionFitGridSize);
            rep(iter).imse_mpl_add_C = additive_imse(B1, X1, bspline_phi, ...
                num_basis,q,active_additive_C,additive_linear_coefficient, ...
                additive_signal_scale,dgp_type,"C",false,censoring_function_design, ...
                opts.FunctionFitGridSize);
            if collect_function_fits
                function_fits = store_function_fit_curves( ...
                    function_fits,"mpl_T",iter,B1,X1,bspline_beta,num_basis, ...
                    active_additive_union);
                function_fits = store_function_fit_curves( ...
                    function_fits,"mpl_C",iter,B1,X1,bspline_phi,num_basis, ...
                    active_additive_union);
            end
            if collect_diagnostics
                diagnostics(iter) = merge_mpl_diagnostics(diagnostics(iter), mpl_diag);
            end
        end

        if any(methods == "PL")
            pl_linear_idx = screen.linear_T(:);
            pl_additive_idx = screen.additive_T(:);
            [linear_present,pl_linear_pos] = ismember(pl_linear_idx,linear_idx);
            [additive_present,pl_additive_pos] = ismember( ...
                pl_additive_idx,additive_idx);
            if ~all(linear_present) || ~all(additive_present)
                error('PLAPH:PLScreenEmbedding', ...
                    'Failure-margin PL screen is not contained in the union design.');
            end
            pl_spline_pos = spline_group_columns(pl_additive_pos,num_basis);
            pl_spline_cols = spline_group_columns(pl_additive_idx,num_basis);
            if lower(string(opts.PLVersion)) == "v2"
                [pl_beta,pl_bspline_beta,pl_refit_diag] = run_pl_v2_tuned( ...
                    n,numel(pl_additive_idx),ini_beta(pl_linear_pos), ...
                    ini_spline_beta(pl_spline_pos),Z(:,pl_linear_idx), ...
                    B(:,pl_spline_cols),Delta,a,threshold_linear_beta, ...
                    linear_kink_coefficient,mic_penalty_scale_T, ...
                    spline_group_iht_threshold_T,mpl_step_scale, ...
                    decomposable_pilot_options,max_outer_iter, ...
                    opts.PostSelectionRefitMaxIter,opts,num_basis, ...
                    linear_penalty_ebic_gamma);
                rep(iter).pl_selection_gate_pass = pl_refit_diag.selection_gate_pass;
                rep(iter).pl_refit_stationarity_pass = pl_refit_diag.stationarity_pass;
                rep(iter).pl_refit_stationarity_inf = pl_refit_diag.stationarity_inf;
                rep(iter).pl_accepted = pl_refit_diag.accepted;
                rep(iter).pl_tuning_cL = pl_refit_diag.tuning_cL;
                rep(iter).pl_tuning_kappa = pl_refit_diag.tuning_kappa;
                rep(iter).pl_step_doublings = pl_refit_diag.step_doublings;
            else
            [pl_beta, pl_bspline_beta, pl_refit_diag] = PL( ...
                n,numel(pl_additive_idx),ini_beta(pl_linear_pos), ...
                ini_spline_beta(pl_spline_pos),Z(:,pl_linear_idx), ...
                X(:,pl_additive_idx),B(:,pl_spline_cols),T,1e-5,Delta, ...
                a,block_size,q_size, ...
                threshold_linear_beta,threshold_spline,mic_penalty_scale_T,mic_penalty_normalization, ...
                spline_threshold_stage,spline_selection_metric,spline_group_penalty_scale_T, ...
                spline_group_penalty_type,spline_group_penalty_normalization,spline_group_mic_sharpness, ...
                spline_group_mic_continuation,spline_group_mic_sharpness_path, ...
                spline_group_mic_form,spline_group_kink_coefficient,spline_group_mcp_concavity, ...
                spline_group_iht_threshold_T,spline_group_iht_warmup_iterations, ...
                linear_mic_form,linear_kink_coefficient, ...
                opts.PostSelectionRefitMaxIter);
            end
            pl_beta = expand_linear_coefficients(pl_beta, pl_linear_idx, p);
            pl_bspline_beta = expand_spline_coefficients( ...
                pl_bspline_beta,pl_additive_idx,num_basis,q);
            pl_linear_support = expand_support_mask( ...
                pl_refit_diag.linear_support,pl_linear_idx,p);
            pl_additive_support = expand_support_mask( ...
                pl_refit_diag.additive_support,pl_additive_idx,q);
            rep(iter).mse_pl = (pl_beta-Beta)' * sigma * (pl_beta-Beta);
            rep(iter) = assign_active_coefficient_fields( ...
                rep(iter),"pl_beta",pl_beta,active_beta);
            rep(iter).pl_fp = support_false_positives( ...
                pl_linear_support,active_beta);
            rep(iter).pl_fn = support_false_negatives( ...
                pl_linear_support,active_beta);
            rep(iter).pl_beta_selected = sum(pl_linear_support);
            rep(iter).pl_beta_min_signal = min(abs(pl_beta(active_beta)));
            rep(iter).pl_beta_max_noise = max_abs_noise(pl_beta, active_beta);
            rep(iter).pl_refit_iterations = pl_refit_diag.iterations;
            rep(iter).pl_refit_converged = pl_refit_diag.converged;
            rep(iter).pl_refit_rejections = pl_refit_diag.rejections;
            rep(iter).pl_refit_final_loglik = pl_refit_diag.final_loglik;
            rep(iter).pl_selected_additive = sum(pl_additive_support);
            rep(iter).pl_add_fp = support_false_positives( ...
                pl_additive_support,active_additive_T);
            rep(iter).pl_add_fn = support_false_negatives( ...
                pl_additive_support,active_additive_T);
            rep(iter).pl_add_active_min = min_active_group_norm(pl_bspline_beta,num_basis,q,active_additive_T);
            rep(iter).pl_add_noise_max = max_noise_group_norm(pl_bspline_beta,num_basis,q,active_additive_T);
            rep(iter).pl_eff_active_min = min_active_effect_norm(B,pl_bspline_beta,num_basis,q,active_additive_T);
            rep(iter).pl_eff_noise_max = max_noise_effect_norm(B,pl_bspline_beta,num_basis,q,active_additive_T);
            rep(iter).pl_eff_noise_q50 = noise_effect_quantile(B,pl_bspline_beta,num_basis,q,active_additive_T,0.50);
            rep(iter).pl_eff_noise_q75 = noise_effect_quantile(B,pl_bspline_beta,num_basis,q,active_additive_T,0.75);
            rep(iter).pl_eff_noise_q90 = noise_effect_quantile(B,pl_bspline_beta,num_basis,q,active_additive_T,0.90);
            rep(iter).pl_eff_noise_q95 = noise_effect_quantile(B,pl_bspline_beta,num_basis,q,active_additive_T,0.95);
            rep(iter) = assign_active_effect_fields(rep(iter), "pl_eff", ...
                active_effect_norms(B,pl_bspline_beta,num_basis,q,active_additive_T));
            rep(iter).imse_pl_add = additive_imse(B, X, pl_bspline_beta, ...
                num_basis,q,active_additive_T,additive_linear_coefficient, ...
                additive_signal_scale,dgp_type,"T",false, ...
                censoring_function_design, ...
                opts.FunctionFitGridSize);
            if collect_function_fits
                function_fits = store_function_fit_curves( ...
                    function_fits,"pl",iter,B,X,pl_bspline_beta,num_basis, ...
                    active_additive_union);
            end
        end

        validate_requested_metrics(methods, rep(iter));
        if rep(iter).mpl_armijo_rejected
            rep(iter).ok = false;
            rep(iter).message = "MPL Armijo outer update rejected";
        else
            rep(iter).ok = true;
        end
    catch ME
        rep(iter).ok = false;
        rep(iter).message = string(ME.message);
        warning('Scenario %s replication %d failed: %s', scenario.scenario_id, iter, ME.message);
    end
end

result = struct();
result.scenario = scenario;
result.replicates = struct2table(rep);
result.summary = summarize_replicates(result.replicates);
if collect_diagnostics
    result.diagnostics = struct2table(diagnostics);
end
if collect_function_fits
    result.function_fits = function_fits;
end
if pre_screening_type ~= "none"
    result.screening = screening;
end
result.linear_penalty_tuning = linear_penalty_tuning;
result.selected_supports = selected_supports;
result.elapsed_seconds = toc;
end


function scenario = normalize_scenario(scenario)
fields = ["scenario_id" "experiment" "generate_copula" "fit_copula" "dgp_type"];
for f = fields
    if isfield(scenario, f)
        scenario.(f) = string(scenario.(f));
    end
end
scenario.methods = parse_method_list(scenario.methods);
if ~isfield(scenario, "additive_linear_coefficient")
    scenario.additive_linear_coefficient = 0.5;
end
if ~isfield(scenario, "linear_signal_scale")
    scenario.linear_signal_scale = 1;
end
if ~isfield(scenario, "additive_signal_scale")
    scenario.additive_signal_scale = 1;
end
if ~isfield(scenario,"dgp_type")
    scenario.dgp_type = "legacy_shared";
end
scenario.dgp_type = lower(string(scenario.dgp_type));
if ~isfield(scenario,"censoring_function_design")
    scenario.censoring_function_design = "current";
end
scenario.censoring_function_design = ...
    lower(string(scenario.censoring_function_design));
if ~isfield(scenario,"failure_scale")
    scenario.failure_scale = 1;
end
if ~isfield(scenario,"censoring_scale")
    scenario.censoring_scale = 5;
end
end


function [time_signature,linear_signature,additive_signature, ...
        basis_signature,event_signature] = paired_data_signatures(T,Z,X,B,Delta)
% Deterministic signatures verify paired data reuse without storing matrices.
time_signature = weighted_numeric_signature(T);
linear_signature = weighted_numeric_signature(Z);
additive_signature = weighted_numeric_signature(X);
basis_signature = weighted_numeric_signature(B);
event_signature = weighted_numeric_signature(double(Delta));
end


function value = weighted_numeric_signature(values)
values = double(values(:));
if isempty(values)
    value = 0;
    return;
end
weights = 1+mod((1:numel(values))',104729);
value = sum(values.*weights,'double');
end


function methods = parse_method_list(raw_methods)
methods = string(raw_methods);
if isscalar(methods) && contains(methods, "|")
    methods = split(methods, "|")';
end
methods = methods(strlength(methods) > 0);
end


function n0 = choose_bin_count(n,target_bin_count)
if nargin >= 2 && ~isempty(target_bin_count)
    n0 = ceil(n/target_bin_count);
    return;
end
n0 = min(100, max(30, round(n/5)));
end


function step_options = mpl_step_options_from_opts(opts)
step_options = default_mpl_step_options(opts.MPLStepScale);
names = {'MPLStepStrategy','MPLBacktrackingFactor','MPLMaxBacktracks', ...
    'MPLAutoExpansion','MPLStepMinScale','MPLStepMaxScale', ...
    'MPLOracleRefitFminconPolish','MPLOracleRefitFminconAlgorithm', ...
    'MPLOracleRefitFminconScaleProblem', ...
    'MPLOracleRefitGradientAudit', ...
    'MPLCollectModelBasedInference'};
targets = {'strategy','backtracking_factor','max_backtracks', ...
    'auto_expansion','min_scale','max_scale','fmincon_polish', ...
    'fmincon_algorithm','fmincon_scale_problem', ...
    'fmincon_gradient_audit', ...
    'collect_model_based_inference'};
for j = 1:numel(names)
    if isfield(opts, names{j}) && ~isempty(opts.(names{j}))
        step_options.(targets{j}) = opts.(names{j});
    end
end
step_options.strategy = lower(string(step_options.strategy));
end


function step_options = default_mpl_step_options(step_scale)
step_options = struct( ...
    'strategy', "fixed", ...
    'backtracking_factor', 2, ...
    'max_backtracks', 12, ...
    'auto_expansion', 1.25, ...
    'min_scale', 0.5, ...
    'max_scale', 16, ...
    'fmincon_polish', false, ...
    'fmincon_algorithm', "sqp", ...
    'fmincon_scale_problem', false, ...
    'fmincon_gradient_audit', false, ...
    'collect_model_based_inference', false, ...
    'step_scale', step_scale);
end


function inference = empty_mpl_inference(p)
inference = struct( ...
    'attempted', false, 'valid', false, 'reason', "not_requested", ...
    'dimension', NaN, 'rank', NaN, 'condition_number', NaN, ...
    'min_eigenvalue', NaN, 'linear_beta_se', NaN(p,1), ...
    'linear_phi_se', NaN(p,1));
end


function pilot_options = decomposable_pilot_options_from_opts(opts)
pilot_options = struct( ...
    'iterations', opts.MPLDecomposablePilotIterations, ...
    'linear_scale', opts.MPLDecomposablePilotLinearScale, ...
    'group_scale', opts.MPLDecomposablePilotGroupScale, ...
    'start', opts.MPLDecomposablePilotStart, ...
    'weighting', opts.MPLDecomposablePilotWeighting, ...
    'weight_gamma', opts.AdaptiveWeightGamma, ...
    'group_weight_gamma', opts.MPLDecomposablePilotGroupWeightGamma, ...
    'weight_epsilon', opts.AdaptiveWeightEpsilon, ...
    'weight_cap', opts.AdaptiveWeightCap);
end


function state = empty_decomposable_pilot_state(p, qs, q, enabled)
state = struct( ...
    'enabled', logical(enabled), 'completed', false, 'iterations', 0, ...
    'beta', zeros(p,1), 'phi', zeros(p,1), ...
    'bspline_beta', zeros(qs,1), 'bspline_phi', zeros(qs,1), ...
    'first_target_beta_candidate', NaN(p,1), ...
    'first_target_beta_threshold', NaN(p,1), ...
    'first_target_phi_candidate', NaN(p,1), ...
    'first_target_phi_threshold', NaN(p,1), ...
    'first_target_group_T_candidate_norm', NaN(q,1), ...
    'first_target_group_C_candidate_norm', NaN(q,1), ...
    'first_target_group_T_threshold', NaN, ...
    'first_target_group_C_threshold', NaN, ...
    'stationarity_beta', NaN, 'stationarity_phi', NaN, ...
    'stationarity_group_T', NaN, 'stationarity_group_C', NaN, ...
    'linear_weight_beta', ones(p,1), 'linear_weight_phi', ones(p,1), ...
    'group_weight_T', ones(q,1), 'group_weight_C', ones(q,1));
end


function lamt = lambda_guess(cr_target, tau)
% Coarse starting values based on the pilot script; Experiment 0 will refine them.
if cr_target <= 0.30
    if tau >= 0.45
        lamt = 0.755;
    else
        lamt = 0.65;
    end
else
    if tau >= 0.45
        lamt = 2.55;
    else
        lamt = 2.32;
    end
end
end


function screen = build_pre_screening(screening_type, T, Z, X, B, status, num_basis, linear_size, additive_size, additive_rank)
% Rank failure and censoring candidates separately, then use their unions.
screening_type = lower(string(screening_type));
p = size(Z, 2);
q = size(B, 2) / num_basis;
if screening_type == "none"
    screen = struct( ...
        'type', screening_type, ...
        'linear_T', (1:p)', 'linear_C', (1:p)', 'linear_union', (1:p)', ...
        'additive_T', (1:q)', 'additive_C', (1:q)', 'additive_union', (1:q)', ...
        'linear_score_T', NaN(p, 1), 'linear_score_C', NaN(p, 1), ...
        'additive_score_T', NaN(q, 1), 'additive_score_C', NaN(q, 1));
    return;
end
if screening_type == "joint_distance_correlation"
    linear_scores = joint_outcome_distance_scores(Z, T, status);
    additive_scores = joint_outcome_distance_scores(X, T, status);
    linear_common = top_score_indices(linear_scores, min(linear_size, p));
    additive_common = top_score_indices(additive_scores, min(additive_size, q));
    screen = struct( ...
        'type', screening_type, ...
        'linear_T', linear_common, 'linear_C', linear_common, 'linear_union', linear_common, ...
        'additive_T', additive_common, 'additive_C', additive_common, 'additive_union', additive_common, ...
        'linear_score_T', linear_scores, 'linear_score_C', linear_scores, ...
        'additive_score_T', additive_scores, 'additive_score_C', additive_scores);
    return;
end
if screening_type == "joint_fused_observed"
    linear_dcor = joint_outcome_distance_scores(Z, T, status);
    additive_dcor = joint_outcome_distance_scores(X, T, status);
    linear_time = observed_log_time_linear_scores(Z, T);
    additive_time = observed_log_time_group_scores(B, T, num_basis, q, additive_rank);
    linear_scores = fuse_screening_scores(linear_dcor, linear_time);
    additive_scores = fuse_screening_scores(additive_dcor, additive_time);
    linear_common = top_score_indices(linear_scores, min(linear_size, p));
    additive_common = top_score_indices(additive_scores, min(additive_size, q));
    screen = struct( ...
        'type', screening_type, ...
        'linear_T', linear_common, 'linear_C', linear_common, 'linear_union', linear_common, ...
        'additive_T', additive_common, 'additive_C', additive_common, 'additive_union', additive_common, ...
        'linear_score_T', linear_scores, 'linear_score_C', linear_scores, ...
        'additive_score_T', additive_scores, 'additive_score_C', additive_scores, ...
        'linear_score_dcor', linear_dcor, 'linear_score_log_time', linear_time, ...
        'additive_score_dcor', additive_dcor, 'additive_score_log_time', additive_time);
    return;
end
if screening_type == "dual_margin_rank_fused_information"
    linear_dcor = joint_outcome_distance_scores(Z, T, status);
    additive_dcor = joint_outcome_distance_scores(X, T, status);
    linear_time = observed_log_time_linear_scores(Z, T);
    additive_time = observed_log_time_group_scores(B, T, num_basis, q, additive_rank);
    linear_observed = fuse_screening_scores(linear_dcor, linear_time);
    additive_observed = fuse_screening_scores(additive_dcor, additive_time);

    linear_cox_T = cox_null_score_statistics(Z, status);
    linear_cox_C = cox_null_score_statistics(Z, ~status);
    additive_cox_T = cox_null_group_information_statistics( ...
        B, status, num_basis, q, additive_rank);
    additive_cox_C = cox_null_group_information_statistics( ...
        B, ~status, num_basis, q, additive_rank);

    linear_score_T = rank_fuse_screening_scores(linear_cox_T, linear_observed);
    linear_score_C = rank_fuse_screening_scores(linear_cox_C, linear_observed);
    additive_score_T = rank_fuse_screening_scores(additive_cox_T, additive_observed);
    additive_score_C = rank_fuse_screening_scores(additive_cox_C, additive_observed);
    linear_T = top_score_indices(linear_score_T, min(linear_size, p));
    linear_C = top_score_indices(linear_score_C, min(linear_size, p));
    additive_T = top_score_indices(additive_score_T, min(additive_size, q));
    additive_C = top_score_indices(additive_score_C, min(additive_size, q));
    screen = struct( ...
        'type', screening_type, ...
        'linear_T', linear_T, 'linear_C', linear_C, ...
        'linear_union', unique([linear_T; linear_C], 'stable'), ...
        'additive_T', additive_T, 'additive_C', additive_C, ...
        'additive_union', unique([additive_T; additive_C], 'stable'), ...
        'linear_score_T', linear_score_T, 'linear_score_C', linear_score_C, ...
        'additive_score_T', additive_score_T, 'additive_score_C', additive_score_C, ...
        'linear_score_cox_T', linear_cox_T, 'linear_score_cox_C', linear_cox_C, ...
        'additive_score_cox_T', additive_cox_T, 'additive_score_cox_C', additive_cox_C, ...
        'linear_score_observed', linear_observed, ...
        'additive_score_observed', additive_observed);
    return;
end
if ~any(screening_type == ["dual_margin_cox_score" "dual_margin_cox_information"])
    error('PLAPH:UnknownScreeningType', 'Unsupported PreScreeningType: %s', screening_type);
end

linear_size = min(linear_size, p);
additive_size = min(additive_size, q);
linear_score_T = cox_null_score_statistics(Z, status);
linear_score_C = cox_null_score_statistics(Z, ~status);
if screening_type == "dual_margin_cox_information"
    additive_score_T = cox_null_group_information_statistics( ...
        B, status, num_basis, q, additive_rank);
    additive_score_C = cox_null_group_information_statistics( ...
        B, ~status, num_basis, q, additive_rank);
else
    additive_score_T = cox_null_group_score_statistics( ...
        B, status, num_basis, q, additive_rank);
    additive_score_C = cox_null_group_score_statistics( ...
        B, ~status, num_basis, q, additive_rank);
end
linear_T = top_score_indices(linear_score_T, linear_size);
linear_C = top_score_indices(linear_score_C, linear_size);
additive_T = top_score_indices(additive_score_T, additive_size);
additive_C = top_score_indices(additive_score_C, additive_size);
screen = struct( ...
    'type', screening_type, ...
    'linear_T', linear_T, 'linear_C', linear_C, ...
    'linear_union', unique([linear_T; linear_C], 'stable'), ...
    'additive_T', additive_T, 'additive_C', additive_C, ...
    'additive_union', unique([additive_T; additive_C], 'stable'), ...
    'linear_score_T', linear_score_T, 'linear_score_C', linear_score_C, ...
    'additive_score_T', additive_score_T, 'additive_score_C', additive_score_C);
end


function scores = joint_outcome_distance_scores(design, T, status)
% Distance correlation with the fully observed pair (log Y, Delta). This is
% model-free and does not treat informative censoring as independent.
log_time = log(max(T(:), realmin));
log_time = (log_time - mean(log_time)) / (std(log_time) + eps);
event = double(status(:));
outcome_distance = sqrt((log_time-log_time').^2 + (event-event').^2);
outcome_centered = double_center_distance(outcome_distance);
outcome_variance = mean(outcome_centered(:).^2);
scores = zeros(size(design, 2), 1);
if outcome_variance <= eps
    return;
end
for j = 1:size(design, 2)
    predictor_distance = abs(design(:, j)-design(:, j)');
    predictor_centered = double_center_distance(predictor_distance);
    predictor_variance = mean(predictor_centered(:).^2);
    if predictor_variance > eps
        covariance_squared = max(mean(predictor_centered(:).*outcome_centered(:)), 0);
        scores(j) = sqrt(covariance_squared / sqrt(predictor_variance*outcome_variance));
    end
end
scores(~isfinite(scores)) = 0;
end


function centered = double_center_distance(distance)
row_mean = mean(distance, 2);
centered = distance - row_mean - row_mean' + mean(row_mean);
end


function scores = observed_log_time_linear_scores(design, T)
y = log(max(T(:), realmin));
y = y - mean(y);
design_centered = design - mean(design, 1);
numerator = abs(design_centered' * y);
denominator = sqrt(sum(design_centered.^2, 1)' * sum(y.^2) + eps);
scores = numerator ./ denominator;
scores(~isfinite(scores)) = 0;
end


function scores = observed_log_time_group_scores(B, T, num_basis, q, screening_rank)
y = log(max(T(:), realmin));
y = y - mean(y);
n = size(B, 1);
scores = zeros(q, 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block = B(:, idx);
    gram = (block' * block) / n;
    gram = (gram + gram') / 2;
    [vectors, values] = eig(gram, 'vector');
    [values, order] = sort(real(values), 'descend');
    tolerance = max(1e-10, max([values; 0])*num_basis*eps);
    keep = find(values > tolerance, min(screening_rank, sum(values > tolerance)), 'first');
    if ~isempty(keep)
        projection = vectors(:, order(keep));
        inverse = projection * diag(1 ./ values(keep)) * projection';
        block_score = (block' * y) / n;
        scores(j) = sqrt(max(real(block_score' * inverse * block_score), 0));
    end
end
scores(~isfinite(scores)) = 0;
end


function fused = fuse_screening_scores(varargin)
dimension = numel(varargin{1});
best_rank = inf(dimension, 1);
for k = 1:nargin
    [~, order] = sort(varargin{k}, 'descend');
    ranks = zeros(dimension, 1);
    ranks(order) = (1:dimension)';
    best_rank = min(best_rank, ranks);
end
fused = dimension - best_rank + 1;
end


function fused = rank_fuse_screening_scores(score_a, score_b)
% Prefer a variable ranked highly by either valid screening statistic.
% The worse of the two ranks provides a deterministic secondary ordering.
dimension = numel(score_a);
[~, order_a] = sort(score_a, 'descend');
[~, order_b] = sort(score_b, 'descend');
rank_a = zeros(dimension, 1);
rank_b = zeros(dimension, 1);
rank_a(order_a) = (1:dimension)';
rank_b(order_b) = (1:dimension)';
best_rank = min(rank_a, rank_b);
second_rank = max(rank_a, rank_b);
fused = (dimension + 1 - best_rank) + ...
    (dimension + 1 - second_rank) / (dimension + 1);
end


function scores = cox_null_score_statistics(design, event)
% Diagonally standardized Cox score at the zero-coefficient model.
n = size(design, 1);
risk_size = (1:n)';
risk_mean = cumsum(design, 1) ./ risk_size;
risk_second = cumsum(design.^2, 1) ./ risk_size;
event = double(event(:));
score = sum(event .* (design - risk_mean), 1);
information = sum(event .* max(risk_second - risk_mean.^2, 0), 1);
scores = (abs(score) ./ sqrt(information + eps))';
scores(~isfinite(scores)) = 0;
end


function scores = cox_null_group_score_statistics(B, event, num_basis, q, screening_rank)
% Low-rank Gram-standardized block score for screening only. The estimator
% retains the full spline basis after candidates have been screened.
n = size(B, 1);
risk_size = (1:n)';
risk_mean = cumsum(B, 1) ./ risk_size;
event = double(event(:));
score = sum(event .* (B - risk_mean), 1)';
scores = zeros(q, 1);
event_count = max(sum(event), 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block_score = score(idx);
    block_design = B(:, idx);
    gram = (block_design' * block_design) / n;
    gram = (gram + gram') / 2;
    [vectors, values] = eig(gram, 'vector');
    [values, order] = sort(real(values), 'descend');
    tolerance = max(1e-10, max([values; 0])*num_basis*eps);
    keep = find(values > tolerance, min(screening_rank, sum(values > tolerance)), 'first');
    if ~isempty(keep)
        projection = vectors(:, order(keep));
        inverse = projection * diag(1 ./ values(keep)) * projection';
        scores(j) = sqrt(max(real(block_score' * inverse * block_score), 0)) / event_count;
    end
end
scores(~isfinite(scores)) = 0;
end


function scores = cox_null_group_information_statistics(B, event, num_basis, q, screening_rank)
% Risk-set information-standardized Cox block score at the null model.
% The data are ordered so rows 1:i form the risk set at row i. Restricting
% the inverse information to its leading estimable eigenspace removes the
% centered spline null direction without adding arbitrary ridge constants.
n = size(B, 1);
risk_size = (1:n)';
event = double(event(:));
scores = zeros(q, 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block = B(:, idx);
    risk_mean = cumsum(block, 1) ./ risk_size;
    block_score = sum(event .* (block-risk_mean), 1)';
    information = zeros(num_basis);
    for a = 1:num_basis
        for b = a:num_basis
            risk_cross = cumsum(block(:,a).*block(:,b)) ./ risk_size;
            value = sum(event .* ...
                (risk_cross-risk_mean(:,a).*risk_mean(:,b)));
            information(a,b) = value;
            information(b,a) = value;
        end
    end
    information = (information+information')/2;
    [vectors, values] = eig(information, 'vector');
    [values, order] = sort(real(values), 'descend');
    tolerance = max(1e-10, max([values; 0])*num_basis*eps);
    keep = find(values > tolerance, ...
        min(screening_rank, sum(values > tolerance)), 'first');
    if ~isempty(keep)
        projection = vectors(:,order(keep));
        inverse = projection*diag(1./values(keep))*projection';
        scores(j) = sqrt(max(real(block_score'*inverse*block_score),0));
    end
end
scores(~isfinite(scores)) = 0;
end


function idx = top_score_indices(scores, count)
candidate = [(1:numel(scores))', -scores(:)];
candidate = sortrows(candidate, [2 1]);
idx = candidate(1:min(count, size(candidate, 1)), 1);
end


function cols = spline_group_columns(groups, num_basis)
cols = zeros(numel(groups)*num_basis, 1);
for j = 1:numel(groups)
    target = ((j-1)*num_basis+1):(j*num_basis);
    source = ((groups(j)-1)*num_basis+1):(groups(j)*num_basis);
    cols(target) = source;
end
end


function full = expand_linear_coefficients(coef, indices, p)
full = zeros(p, 1);
full(indices) = coef;
end


function full = expand_linear_standard_errors(standard_error,indices,p)
full = NaN(p,1);
if ~isempty(standard_error)
    full(indices) = standard_error;
end
end


function full = expand_linear_thresholds(threshold, indices, p)
full = Inf(p,1);
full(indices) = threshold;
end


function full = expand_group_values(values, groups, q)
full = zeros(q,1);
full(groups) = values;
end


function full = expand_support_mask(mask,indices,total_dimension)
full = false(total_dimension,1);
full(indices) = logical(mask(:));
end


function full = expand_spline_coefficients(coef, groups, num_basis, q)
full = zeros(q*num_basis, 1);
for j = 1:numel(groups)
    source = ((j-1)*num_basis+1):(j*num_basis);
    target = ((groups(j)-1)*num_basis+1):(groups(j)*num_basis);
    full(target) = coef(source);
end
end


function rep = assign_screening_fields(rep,screen,active_beta,active_phi, ...
    active_additive_T,active_additive_C)
rep.screen_linear_T_size = numel(screen.linear_T);
rep.screen_linear_C_size = numel(screen.linear_C);
rep.screen_linear_union_size = numel(screen.linear_union);
rep.screen_additive_T_size = numel(screen.additive_T);
rep.screen_additive_C_size = numel(screen.additive_C);
rep.screen_additive_union_size = numel(screen.additive_union);
rep.screen_beta_sure = double(all(ismember(active_beta, screen.linear_T)));
rep.screen_phi_sure = double(all(ismember(active_phi, screen.linear_C)));
rep.screen_add_T_sure = double(all(ismember(active_additive_T,screen.additive_T)));
rep.screen_add_C_sure = double(all(ismember(active_additive_C,screen.additive_C)));
rep.screen_all_sure = double(rep.screen_beta_sure && rep.screen_phi_sure && ...
    rep.screen_add_T_sure && rep.screen_add_C_sure);
end


function rep = empty_replication_result()
rep = struct( ...
    'replication', NaN, 'seed', NaN, 'ok', false, 'message', "", ...
    'mpl_selection_gate_pass', false, 'mpl_selection_status_code', NaN, ...
    'mpl_selection_gate_maximum', NaN, 'mpl_refit_stationarity_pass', false, ...
    'mpl_refit_stationarity_inf', NaN, 'mpl_accepted', false, ...
    'mpl_step_doublings', NaN, 'mpl_target_iterations', NaN, ...
    'mpl_df_rank', NaN, 'tuning_cL', NaN, 'tuning_kappa', NaN, ...
    'tuning_ebic', NaN, 'tuning_valid_count', NaN, 'tuning_at_boundary', false, ...
    'tuning_tie_count', NaN, ...
    'mpl_pilot_final_agree_linear_T', NaN, 'mpl_pilot_final_agree_linear_C', NaN, ...
    'mpl_pilot_final_agree_add_T', NaN, 'mpl_pilot_final_agree_add_C', NaN, ...
    'pl_selection_gate_pass', false, 'pl_refit_stationarity_pass', false, ...
    'pl_refit_stationarity_inf', NaN, 'pl_accepted', false, ...
    'pl_tuning_cL', NaN, 'pl_tuning_kappa', NaN, 'pl_step_doublings', NaN, ...
    'spline_order', NaN, ...
    'generate_tau', NaN, 'fit_tau', NaN, ...
    'data_time_signature', NaN, 'data_linear_signature', NaN, ...
    'data_additive_signature', NaN, 'data_basis_signature', NaN, ...
    'data_event_signature', NaN, ...
    'linear_penalty_tuned', false, ...
    'linear_penalty_tuning_criterion', "fixed", ...
    'linear_penalty_candidate_count', NaN, ...
    'linear_penalty_valid_candidate_count', NaN, ...
    'linear_penalty_selected_multiplier', NaN, ...
    'linear_penalty_selected_lambda_T', NaN, ...
    'linear_penalty_selected_lambda_C', NaN, ...
    'linear_penalty_selected_criterion', NaN, ...
    'linear_penalty_selected_bic', NaN, ...
    'linear_penalty_selected_ebic', NaN, ...
    'linear_penalty_selected_df', NaN, ...
    'linear_penalty_selected_at_grid_boundary', false, ...
    'censor_rate', NaN, 'mse_mpl_beta', NaN, 'mse_mpl_phi', NaN, ...
    'n_failure_events', NaN, 'n_censoring_events', NaN, ...
    'censor_rate_x2_center', NaN, 'censor_rate_x2_tail', NaN, ...
    'x2_tail_to_center_censor_ratio', NaN, ...
    'screen_linear_T_size', NaN, 'screen_linear_C_size', NaN, ...
    'screen_linear_union_size', NaN, 'screen_additive_T_size', NaN, ...
    'screen_additive_C_size', NaN, 'screen_additive_union_size', NaN, ...
    'screen_beta_sure', NaN, 'screen_phi_sure', NaN, ...
    'screen_add_T_sure', NaN, 'screen_add_C_sure', NaN, 'screen_all_sure', NaN, ...
    'mse_pl', NaN, 'mse_oracle', NaN, 'tau_hat', NaN, 'alpha_hat', NaN, ...
    'mpl_baseline_armijo_failed', false, 'mpl_armijo_rejected', false, ...
    'beta_fp', NaN, 'beta_fn', NaN, 'phi_fp', NaN, 'phi_fn', NaN, ...
    'pl_fp', NaN, 'pl_fn', NaN, 'ini_beta_selected', NaN, ...
    'pl_refit_iterations', NaN, 'pl_refit_converged', false, ...
    'pl_refit_rejections', NaN, 'pl_refit_final_loglik', NaN, ...
    'ini_phi_selected', NaN, 'ini_beta_min_signal', NaN, ...
    'ini_phi_min_signal', NaN, 'ini_beta_max_noise', NaN, ...
    'ini_phi_max_noise', NaN, ...
    'mpl_pilot_completed', false, 'mpl_pilot_iterations', 0, ...
    'mpl_pilot_beta_fp', NaN, 'mpl_pilot_beta_fn', NaN, ...
    'mpl_pilot_phi_fp', NaN, 'mpl_pilot_phi_fn', NaN, ...
    'mpl_pilot_beta_min_signal', NaN, 'mpl_pilot_phi_min_signal', NaN, ...
    'mpl_pilot_beta_max_noise', NaN, 'mpl_pilot_phi_max_noise', NaN, ...
    'mpl_pilot_add_T_fp', NaN, 'mpl_pilot_add_T_fn', NaN, ...
    'mpl_pilot_add_C_fp', NaN, 'mpl_pilot_add_C_fn', NaN, ...
    'mpl_pilot_eff_T_active_min', NaN, 'mpl_pilot_eff_T_noise_max', NaN, ...
    'mpl_pilot_eff_C_active_min', NaN, 'mpl_pilot_eff_C_noise_max', NaN, ...
    'mpl_pilot_stationarity_beta', NaN, 'mpl_pilot_stationarity_phi', NaN, ...
    'mpl_pilot_stationarity_group_T', NaN, ...
    'mpl_pilot_stationarity_group_C', NaN, ...
    'mpl_target_beta_candidate_active_min', NaN, ...
    'mpl_target_beta_candidate_noise_max', NaN, ...
    'mpl_target_beta_active_margin_min', NaN, ...
    'mpl_target_beta_noise_margin_max', NaN, ...
    'mpl_target_beta_threshold_min', NaN, 'mpl_target_beta_threshold_max', NaN, ...
    'mpl_target_beta_strict_separation', NaN, ...
    'mpl_target_phi_candidate_active_min', NaN, ...
    'mpl_target_phi_candidate_noise_max', NaN, ...
    'mpl_target_phi_active_margin_min', NaN, ...
    'mpl_target_phi_noise_margin_max', NaN, ...
    'mpl_target_phi_threshold_min', NaN, 'mpl_target_phi_threshold_max', NaN, ...
    'mpl_target_phi_strict_separation', NaN, ...
    'mpl_target_group_T_candidate_active_min', NaN, ...
    'mpl_target_group_T_candidate_noise_max', NaN, ...
    'mpl_target_group_T_active_margin_min', NaN, ...
    'mpl_target_group_T_noise_margin_max', NaN, ...
    'mpl_target_group_T_threshold', NaN, ...
    'mpl_target_group_T_strict_separation', NaN, ...
    'mpl_target_group_C_candidate_active_min', NaN, ...
    'mpl_target_group_C_candidate_noise_max', NaN, ...
    'mpl_target_group_C_active_margin_min', NaN, ...
    'mpl_target_group_C_noise_margin_max', NaN, ...
    'mpl_target_group_C_threshold', NaN, ...
    'mpl_target_group_C_strict_separation', NaN, ...
    'mpl_beta_selected', NaN, ...
    'mpl_phi_selected', NaN, 'mpl_beta_min_signal', NaN, ...
    'mpl_phi_min_signal', NaN, 'mpl_beta_max_noise', NaN, ...
    'mpl_phi_max_noise', NaN, 'pl_beta_selected', NaN, ...
    'pl_beta_min_signal', NaN, 'pl_beta_max_noise', NaN, ...
    'mpl_beta_noise_q50', NaN, 'mpl_beta_noise_q75', NaN, ...
    'mpl_beta_noise_q90', NaN, 'mpl_beta_noise_q95', NaN, ...
    'mpl_phi_noise_q50', NaN, 'mpl_phi_noise_q75', NaN, ...
    'mpl_phi_noise_q90', NaN, 'mpl_phi_noise_q95', NaN, ...
    'mpl_selected_additive_T', NaN, ...
    'mpl_selected_additive_C', NaN, 'pl_selected_additive', NaN, ...
    'oracle_selected_additive', NaN, ...
    'mpl_add_T_fp', NaN, 'mpl_add_T_fn', NaN, ...
    'mpl_add_C_fp', NaN, 'mpl_add_C_fn', NaN, ...
    'pl_add_fp', NaN, 'pl_add_fn', NaN, ...
    'mpl_add_T_active_min', NaN, 'mpl_add_T_noise_max', NaN, ...
    'mpl_add_C_active_min', NaN, 'mpl_add_C_noise_max', NaN, ...
    'pl_add_active_min', NaN, 'pl_add_noise_max', NaN, ...
    'mpl_eff_T_active_min', NaN, 'mpl_eff_T_noise_max', NaN, ...
    'mpl_eff_C_active_min', NaN, 'mpl_eff_C_noise_max', NaN, ...
    'pl_eff_active_min', NaN, 'pl_eff_noise_max', NaN, ...
    'mpl_eff_T_noise_q50', NaN, 'mpl_eff_T_noise_q75', NaN, ...
    'mpl_eff_T_noise_q90', NaN, 'mpl_eff_T_noise_q95', NaN, ...
    'mpl_eff_C_noise_q50', NaN, 'mpl_eff_C_noise_q75', NaN, ...
    'mpl_eff_C_noise_q90', NaN, 'mpl_eff_C_noise_q95', NaN, ...
    'pl_eff_noise_q50', NaN, 'pl_eff_noise_q75', NaN, ...
    'pl_eff_noise_q90', NaN, 'pl_eff_noise_q95', NaN, ...
    'mpl_eff_T_g1', NaN, 'mpl_eff_T_g2', NaN, ...
    'mpl_eff_T_g3', NaN, 'mpl_eff_T_g4', NaN, ...
    'mpl_eff_C_g1', NaN, 'mpl_eff_C_g2', NaN, ...
    'mpl_eff_C_g3', NaN, 'mpl_eff_C_g4', NaN, ...
    'pl_eff_g1', NaN, 'pl_eff_g2', NaN, ...
    'pl_eff_g3', NaN, 'pl_eff_g4', NaN, ...
    'imse_mpl_add_T', NaN, 'imse_mpl_add_C', NaN, ...
    'imse_pl_add', NaN, 'imse_oracle_add', NaN, ...
    'mpl_beta_active1', NaN, 'mpl_beta_active2', NaN, ...
    'mpl_beta_active3', NaN, 'mpl_phi_active1', NaN, ...
    'mpl_phi_active2', NaN, 'mpl_phi_active3', NaN, ...
    'mpl_beta_inactive1', NaN, 'mpl_beta_inactive1_index', NaN, ...
    'mpl_beta_ase_active1', NaN, 'mpl_beta_ase_active2', NaN, ...
    'mpl_beta_ase_active3', NaN, 'mpl_phi_ase_active1', NaN, ...
    'mpl_phi_ase_active2', NaN, 'mpl_phi_ase_active3', NaN, ...
    'mpl_beta_cp_active1', NaN, 'mpl_beta_cp_active2', NaN, ...
    'mpl_beta_cp_active3', NaN, 'mpl_phi_cp_active1', NaN, ...
    'mpl_phi_cp_active2', NaN, 'mpl_phi_cp_active3', NaN, ...
    'pl_beta_active1', NaN, 'pl_beta_active2', NaN, ...
    'pl_beta_active3', NaN, 'oracle_beta_active1', NaN, ...
    'oracle_beta_active2', NaN, 'oracle_beta_active3', NaN);
end


function rep = assign_active_coefficient_fields(rep,prefix,coef,active)
for j = 1:min(3,numel(active))
    rep.(sprintf('%s_active%d',prefix,j)) = coef(active(j));
end
end


function rep = assign_active_inference_fields( ...
        rep,prefix,standard_error,estimate,truth,active,inference_valid)
for j = 1:min(3,numel(active))
    index = active(j);
    se = standard_error(index);
    rep.(sprintf('%s_ase_active%d',prefix,j)) = se;
    if ~inference_valid
        coverage = NaN;
    elseif isfinite(se) && se >= 0
        coverage = double(abs(estimate(index)-truth(index)) <= 1.96*se);
    else
        % A true signal omitted by selection has no post-selection interval.
        coverage = 0;
    end
    rep.(sprintf('%s_cp_active%d',prefix,j)) = coverage;
end
end


function diag = empty_diagnostic_result()
diag = struct( ...
    'replication', NaN, 'seed', NaN, 'n', NaN, 'p', NaN, 'q', NaN, ...
    'num_basis', NaN, 'event_count', NaN, 'censor_rate', NaN, ...
    'T_min', NaN, 'T_max', NaN, 'Z_all_finite', false, ...
    'X_all_finite', false, 'B_all_finite', false, 'Z_max_abs', NaN, ...
    'X_max_abs', NaN, 'B_max_abs', NaN, 'B_col_mean_max_abs', NaN, ...
    'ini_beta_selected', NaN, 'ini_phi_selected', NaN, ...
    'ini_beta_min_signal', NaN, 'ini_phi_min_signal', NaN, ...
    'ini_beta_max_noise', NaN, 'ini_phi_max_noise', NaN, ...
    'ini_spline_beta_max_abs', NaN, 'ini_spline_phi_max_abs', NaN, ...
    'ini_spline_beta_eff_g1', NaN, 'ini_spline_beta_eff_g2', NaN, ...
    'ini_spline_beta_eff_g3', NaN, 'ini_spline_beta_eff_g4', NaN, ...
    'ini_spline_phi_eff_g1', NaN, 'ini_spline_phi_eff_g2', NaN, ...
    'ini_spline_phi_eff_g3', NaN, 'ini_spline_phi_eff_g4', NaN, ...
    'mpl_baseline_armijo_failed', false, 'mpl_armijo_rejected', false, ...
    'mpl_diag_collected', false, 'mpl_first_L_before_coef', NaN, ...
    'mpl_first_L_after_coef', NaN, 'mpl_first_coef_delta', NaN, ...
    'mpl_start_etaT_min', NaN, 'mpl_start_etaT_max', NaN, ...
    'mpl_start_etaC_min', NaN, 'mpl_start_etaC_max', NaN, ...
    'mpl_start_exp_etaT_max', NaN, 'mpl_start_exp_etaC_max', NaN, ...
    'mpl_start_linear_etaT_max_abs', NaN, 'mpl_start_spline_etaT_max_abs', NaN, ...
    'mpl_start_linear_etaC_max_abs', NaN, 'mpl_start_spline_etaC_max_abs', NaN, ...
    'mpl_first_DL_seta_norm', NaN, 'mpl_first_DL_gama_norm', NaN, ...
    'mpl_first_Ps1_min', NaN, 'mpl_first_Ps2_min', NaN, ...
    'mpl_first_Ps1_nonpos', NaN, 'mpl_first_Ps2_nonpos', NaN, ...
    'mpl_first_seta_min', NaN, 'mpl_first_gama_min', NaN, ...
    'mpl_first_HT_max', NaN, 'mpl_first_HC_max', NaN, ...
    'mpl_first_ST_min', NaN, 'mpl_first_SC_min', NaN, ...
    'mpl_first_etaT_min', NaN, 'mpl_first_etaT_max', NaN, ...
    'mpl_first_etaC_min', NaN, 'mpl_first_etaC_max', NaN, ...
    'mpl_first_exp_etaT_max', NaN, 'mpl_first_exp_etaC_max', NaN, ...
    'mpl_first_K1_min', NaN, 'mpl_first_K2_min', NaN, ...
    'mpl_first_LAM1_abs_max', NaN, 'mpl_first_LAM2_abs_max', NaN, ...
    'mpl_first_LAM11_abs_max', NaN, 'mpl_first_LAM22_abs_max', NaN, ...
    'mpl_first_precond_plus_1e4', NaN, 'mpl_first_precond_minus_1e4', NaN, ...
    'mpl_first_raw_plus_1e4', NaN, 'mpl_first_raw_minus_1e4', NaN, ...
    'mpl_first_precond_plus_1e2', NaN, 'mpl_first_precond_minus_1e2', NaN, ...
    'mpl_first_raw_plus_1e2', NaN, 'mpl_first_raw_minus_1e2', NaN, ...
    'mpl_final_mean_loglik', NaN, 'mpl_final_total_loglik', NaN, ...
    'mpl_linear_mic_lambda_T', NaN, 'mpl_linear_mic_lambda_C', NaN, ...
    'mpl_profile_tau_enabled', false, ...
    'mpl_profile_tau_attempts', 0, 'mpl_profile_tau_accepted', 0, ...
    'mpl_profile_tau_rejected', 0, 'mpl_profile_tau_failures', 0, ...
    'mpl_profile_tau_boundary_hit', false, ...
    'mpl_profile_tau_min_loglik_delta', NaN, ...
    'mpl_joint_refit_enabled', false, ...
    'mpl_joint_refit_initial_loglik', NaN, 'mpl_joint_refit_final_loglik', NaN, ...
    'mpl_joint_refit_iterations', 0, 'mpl_joint_refit_rejected', false, ...
    'mpl_oracle_refit_enabled', false, 'mpl_oracle_refit_converged', false, ...
    'mpl_oracle_refit_coefficient_rejections', 0, ...
    'mpl_oracle_refit_baseline_rejections', 0, ...
    'mpl_oracle_refit_alpha_updates', 0, ...
    'mpl_oracle_refit_alpha_rejections', 0, ...
    'mpl_oracle_refit_alpha_boundary_hit', false, ...
    'mpl_oracle_refit_coefficient_backtracks', 0, ...
    'mpl_oracle_refit_baseline_backtracks', 0, ...
    'mpl_oracle_refit_final_coefficient_step', NaN, ...
    'mpl_oracle_refit_final_baseline_step', NaN, ...
    'mpl_oracle_refit_min_loglik_delta', NaN, ...
    'mpl_oracle_refit_final_relative_change', NaN, ...
    'mpl_oracle_refit_final_relative_loglik_change', NaN, ...
    'mpl_oracle_refit_final_coefficient_score_inf', NaN, ...
    'mpl_oracle_refit_final_baseline_kkt_inf', NaN, ...
    'mpl_oracle_refit_final_stationarity_inf', NaN, ...
    'mpl_oracle_refit_stationarity_pass', false, ...
    'mpl_oracle_refit_stagnation_handoff', false, ...
    'mpl_oracle_refit_fmincon_attempted', false, ...
    'mpl_oracle_refit_fmincon_accepted', false, ...
    'mpl_oracle_refit_fmincon_exitflag', NaN, ...
    'mpl_oracle_refit_fmincon_iterations', NaN, ...
    'mpl_oracle_refit_fmincon_firstorderopt', NaN, ...
    'mpl_oracle_refit_gradient_check_max_abs', NaN, ...
    'mpl_oracle_refit_gradient_check_max_rel', NaN, ...
    'mpl_oracle_refit_gradient_check_max_index', NaN, ...
    'mpl_oracle_refit_gradient_check_analytic', NaN, ...
    'mpl_oracle_refit_gradient_check_numeric', NaN, ...
    'mpl_oracle_refit_gradient_check_block', "not_audited", ...
    'mpl_inference_attempted', false, 'mpl_inference_valid', false, ...
    'mpl_inference_reason', "not_requested", ...
    'mpl_inference_dimension', NaN, 'mpl_inference_rank', NaN, ...
    'mpl_inference_condition_number', NaN, ...
    'mpl_inference_min_eigenvalue', NaN, ...
    'mpl_decomposable_pilot_enabled', false, ...
    'mpl_decomposable_pilot_iterations', 0, ...
    'mpl_decomposable_pilot_completed', false, ...
    'mpl_decomposable_pilot_linear_lambda_T', NaN, ...
    'mpl_decomposable_pilot_linear_lambda_C', NaN, ...
    'mpl_decomposable_pilot_group_lambda_T', NaN, ...
    'mpl_decomposable_pilot_group_lambda_C', NaN, ...
    'mpl_decomposable_pilot_objective_initial', NaN, ...
    'mpl_decomposable_pilot_objective_final', NaN, ...
    'mpl_decomposable_pilot_objective_min_delta', NaN, ...
    'mpl_decomposable_pilot_objective_rejected', false, ...
    'mpl_decomposable_pilot_stationarity_beta', NaN, ...
    'mpl_decomposable_pilot_stationarity_phi', NaN, ...
    'mpl_decomposable_pilot_stationarity_group_T', NaN, ...
    'mpl_decomposable_pilot_stationarity_group_C', NaN, ...
    'mpl_group_iht_enabled', false, ...
    'mpl_group_iht_warmup_iterations', NaN, ...
    'mpl_group_iht_continuation_iterations', NaN, ...
    'mpl_group_iht_start_fraction', NaN, ...
    'mpl_group_iht_threshold_T', NaN, 'mpl_group_iht_threshold_C', NaN, ...
    'mpl_group_iht_lambda_T', NaN, 'mpl_group_iht_lambda_C', NaN, ...
    'mpl_group_iht_objective_initial', NaN, 'mpl_group_iht_objective_final', NaN, ...
    'mpl_group_iht_objective_min_delta', NaN, ...
    'mpl_group_iht_objective_rejected', false, ...
    'mpl_group_iht_support_changes', 0, ...
    'mpl_group_iht_support_stable_iterations', 0, ...
    'mpl_step_strategy', "fixed", ...
    'mpl_step_backtracks_beta', 0, 'mpl_step_backtracks_phi', 0, ...
    'mpl_step_backtracks_spline_T', 0, 'mpl_step_backtracks_spline_C', 0, ...
    'mpl_step_fallbacks_beta', 0, 'mpl_step_fallbacks_phi', 0, ...
    'mpl_step_fallbacks_spline_T', 0, 'mpl_step_fallbacks_spline_C', 0, ...
    'mpl_step_final_L_beta', NaN, 'mpl_step_final_L_phi', NaN, ...
    'mpl_step_final_L_spline_T', NaN, 'mpl_step_final_L_spline_C', NaN, ...
    'mpl_outer_iterations', 0, ...
    'mpl_selection_gate_pass', false, 'mpl_selection_gate_maximum', NaN, ...
    'mpl_selection_status_code', NaN, 'mpl_step_doublings', 0, ...
    'mpl_target_iterations', NaN, ...
    'mpl_oracle_refit_coefficient_score_inf_raw', NaN);
end


function diag = data_diagnostics(iter, seed, n, p, q, Z, X, B, Delta, T, num_basis)
diag = empty_diagnostic_result();
diag.replication = iter;
diag.seed = seed;
diag.n = n;
diag.p = p;
diag.q = q;
diag.num_basis = num_basis;
diag.event_count = sum(Delta);
diag.censor_rate = 1 - mean(Delta);
diag.T_min = min(T);
diag.T_max = max(T);
diag.Z_all_finite = all(isfinite(Z(:)));
diag.X_all_finite = all(isfinite(X(:)));
diag.B_all_finite = all(isfinite(B(:)));
diag.Z_max_abs = max(abs(Z(:)));
diag.X_max_abs = max(abs(X(:)));
diag.B_max_abs = max(abs(B(:)));
diag.B_col_mean_max_abs = max(abs(mean(B,1)));
end


function diag = merge_mpl_diagnostics(diag, mpl_diag)
names = fieldnames(mpl_diag);
for j = 1:numel(names)
    diag.(names{j}) = mpl_diag.(names{j});
end
end


function nsel = count_selected_linear(est, threshold_linear)
nsel = sum(abs(est) > threshold_linear);
end


function val = max_abs_noise(est, active)
noise = setdiff((1:numel(est))', active(:));
if isempty(noise)
    val = NaN;
else
    val = max(abs(est(noise)));
end
end


function val = max_noise_value(values, active)
noise = setdiff((1:numel(values))', active(:));
if isempty(noise)
    val = NaN;
else
    val = max(values(noise));
end
end


function val = finite_min(values)
values = values(isfinite(values));
if isempty(values)
    val = NaN;
else
    val = min(values);
end
end


function val = finite_max(values)
values = values(isfinite(values));
if isempty(values)
    val = NaN;
else
    val = max(values);
end
end


function val = abs_noise_quantile(est, active, prob)
noise = setdiff((1:numel(est))', active(:));
if isempty(noise)
    val = NaN;
else
    vals = sort(abs(est(noise)));
    idx = max(1, min(numel(vals), ceil(prob*numel(vals))));
    val = vals(idx);
end
end


function [bic,ebic,df,s_linear_T,s_linear_C,s_additive_T,s_additive_C] = ...
    linear_penalty_information_criteria(total_loglik,n,n0,num_basis, ...
    pilot_state,linear_mask_T,linear_mask_C,group_mask_T,group_mask_C, ...
    ebic_gamma,Z,B)
% Compare SCAD/MCP fits using the refitted observed-data likelihood. The
% baseline contribution is constant over the candidate path but is retained
% in the reported degrees of freedom. Centering removes one spline direction.
required = {'final_selected_linear_T_mask', ...
    'final_selected_linear_C_mask','final_selected_additive_T_mask', ...
    'final_selected_additive_C_mask'};
if ~all(isfield(pilot_state,required))
    error('PLAPH:LinearPenaltyTuningSupport', ...
        'The fixed-support masks required for BIC/EBIC are unavailable.');
end
s_linear_T = sum(pilot_state.final_selected_linear_T_mask);
s_linear_C = sum(pilot_state.final_selected_linear_C_mask);
s_additive_T = sum(pilot_state.final_selected_additive_T_mask);
s_additive_C = sum(pilot_state.final_selected_additive_C_mask);
baseline_df = 2*ceil(n/n0);
if nargin >= 12 && ~isempty(Z) && ~isempty(B)
    % v2: the regression degrees of freedom in each margin are the numerical
    % rank of the selected linear-plus-spline design, so duplicated covariate
    % profiles and the centering null direction are not double counted.
    df = joint_selected_rank(Z,B,num_basis, ...
        pilot_state.final_selected_linear_T_mask, ...
        pilot_state.final_selected_additive_T_mask) + ...
        joint_selected_rank(Z,B,num_basis, ...
        pilot_state.final_selected_linear_C_mask, ...
        pilot_state.final_selected_additive_C_mask) + baseline_df;
else
    spline_rank = max(num_basis-1,1);
    df = s_linear_T+s_linear_C + ...
        spline_rank*(s_additive_T+s_additive_C) + baseline_df;
end
bic = -2*total_loglik + log(n)*df;
model_space = log_binomial(sum(linear_mask_T),s_linear_T) + ...
    log_binomial(sum(linear_mask_C),s_linear_C) + ...
    log_binomial(sum(group_mask_T),s_additive_T) + ...
    log_binomial(sum(group_mask_C),s_additive_C);
ebic = bic + 2*ebic_gamma*model_space;
end


function value = log_binomial(n,k)
if k < 0 || k > n || n < 0
    value = Inf;
else
    value = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1);
end
end


function validate_requested_metrics(methods, rep)
if (any(methods == "MPL_fixed_tau") || any(methods == "MPL_profile_tau")) && ...
        (~isfinite(rep.mse_mpl_beta) || ~isfinite(rep.mse_mpl_phi))
    error('PLAPH:NonfiniteMPL', 'MPL produced non-finite MSE values.');
end

if any(methods == "PL") && ~isfinite(rep.mse_pl)
    error('PLAPH:NonfinitePL', 'PL produced a non-finite MSE value.');
end

if any(methods == "Oracle") && ~isfinite(rep.mse_oracle)
    error('PLAPH:NonfiniteOracle', 'Oracle produced a non-finite MSE value.');
end
end


function diag = empty_mpl_diagnostic()
diag = struct( ...
    'mpl_diag_collected', false, 'mpl_first_L_before_coef', NaN, ...
    'mpl_baseline_armijo_failed', false, 'mpl_armijo_rejected', false, ...
    'mpl_first_L_after_coef', NaN, 'mpl_first_coef_delta', NaN, ...
    'mpl_start_etaT_min', NaN, 'mpl_start_etaT_max', NaN, ...
    'mpl_start_etaC_min', NaN, 'mpl_start_etaC_max', NaN, ...
    'mpl_start_exp_etaT_max', NaN, 'mpl_start_exp_etaC_max', NaN, ...
    'mpl_start_linear_etaT_max_abs', NaN, 'mpl_start_spline_etaT_max_abs', NaN, ...
    'mpl_start_linear_etaC_max_abs', NaN, 'mpl_start_spline_etaC_max_abs', NaN, ...
    'mpl_first_DL_seta_norm', NaN, 'mpl_first_DL_gama_norm', NaN, ...
    'mpl_first_Ps1_min', NaN, 'mpl_first_Ps2_min', NaN, ...
    'mpl_first_Ps1_nonpos', NaN, 'mpl_first_Ps2_nonpos', NaN, ...
    'mpl_first_seta_min', NaN, 'mpl_first_gama_min', NaN, ...
    'mpl_first_HT_max', NaN, 'mpl_first_HC_max', NaN, ...
    'mpl_first_ST_min', NaN, 'mpl_first_SC_min', NaN, ...
    'mpl_first_etaT_min', NaN, 'mpl_first_etaT_max', NaN, ...
    'mpl_first_etaC_min', NaN, 'mpl_first_etaC_max', NaN, ...
    'mpl_first_exp_etaT_max', NaN, 'mpl_first_exp_etaC_max', NaN, ...
    'mpl_first_K1_min', NaN, 'mpl_first_K2_min', NaN, ...
    'mpl_first_LAM1_abs_max', NaN, 'mpl_first_LAM2_abs_max', NaN, ...
    'mpl_first_LAM11_abs_max', NaN, 'mpl_first_LAM22_abs_max', NaN, ...
    'mpl_first_precond_plus_1e4', NaN, 'mpl_first_precond_minus_1e4', NaN, ...
    'mpl_first_raw_plus_1e4', NaN, 'mpl_first_raw_minus_1e4', NaN, ...
    'mpl_first_precond_plus_1e2', NaN, 'mpl_first_precond_minus_1e2', NaN, ...
    'mpl_first_raw_plus_1e2', NaN, 'mpl_first_raw_minus_1e2', NaN, ...
    'mpl_final_mean_loglik', NaN, 'mpl_final_total_loglik', NaN, ...
    'mpl_linear_mic_lambda_T', NaN, 'mpl_linear_mic_lambda_C', NaN, ...
    'mpl_profile_tau_enabled', false, ...
    'mpl_profile_tau_attempts', 0, 'mpl_profile_tau_accepted', 0, ...
    'mpl_profile_tau_rejected', 0, 'mpl_profile_tau_failures', 0, ...
    'mpl_profile_tau_boundary_hit', false, ...
    'mpl_profile_tau_min_loglik_delta', NaN, ...
    'mpl_joint_refit_enabled', false, ...
    'mpl_joint_refit_initial_loglik', NaN, 'mpl_joint_refit_final_loglik', NaN, ...
    'mpl_joint_refit_iterations', 0, 'mpl_joint_refit_rejected', false, ...
    'mpl_oracle_refit_enabled', false, 'mpl_oracle_refit_converged', false, ...
    'mpl_oracle_refit_coefficient_rejections', 0, ...
    'mpl_oracle_refit_baseline_rejections', 0, ...
    'mpl_oracle_refit_alpha_updates', 0, ...
    'mpl_oracle_refit_alpha_rejections', 0, ...
    'mpl_oracle_refit_alpha_boundary_hit', false, ...
    'mpl_oracle_refit_coefficient_backtracks', 0, ...
    'mpl_oracle_refit_baseline_backtracks', 0, ...
    'mpl_oracle_refit_final_coefficient_step', NaN, ...
    'mpl_oracle_refit_final_baseline_step', NaN, ...
    'mpl_oracle_refit_min_loglik_delta', NaN, ...
    'mpl_oracle_refit_final_relative_change', NaN, ...
    'mpl_oracle_refit_final_relative_loglik_change', NaN, ...
    'mpl_oracle_refit_final_coefficient_score_inf', NaN, ...
    'mpl_oracle_refit_final_baseline_kkt_inf', NaN, ...
    'mpl_oracle_refit_final_stationarity_inf', NaN, ...
    'mpl_oracle_refit_stationarity_pass', false, ...
    'mpl_oracle_refit_stagnation_handoff', false, ...
    'mpl_oracle_refit_fmincon_attempted', false, ...
    'mpl_oracle_refit_fmincon_accepted', false, ...
    'mpl_oracle_refit_fmincon_exitflag', NaN, ...
    'mpl_oracle_refit_fmincon_iterations', NaN, ...
    'mpl_oracle_refit_fmincon_firstorderopt', NaN, ...
    'mpl_oracle_refit_gradient_check_max_abs', NaN, ...
    'mpl_oracle_refit_gradient_check_max_rel', NaN, ...
    'mpl_oracle_refit_gradient_check_max_index', NaN, ...
    'mpl_oracle_refit_gradient_check_analytic', NaN, ...
    'mpl_oracle_refit_gradient_check_numeric', NaN, ...
    'mpl_oracle_refit_gradient_check_block', "not_audited", ...
    'mpl_inference_attempted', false, 'mpl_inference_valid', false, ...
    'mpl_inference_reason', "not_requested", ...
    'mpl_inference_dimension', NaN, 'mpl_inference_rank', NaN, ...
    'mpl_inference_condition_number', NaN, ...
    'mpl_inference_min_eigenvalue', NaN, ...
    'mpl_decomposable_pilot_enabled', false, ...
    'mpl_decomposable_pilot_iterations', 0, ...
    'mpl_decomposable_pilot_completed', false, ...
    'mpl_decomposable_pilot_linear_lambda_T', NaN, ...
    'mpl_decomposable_pilot_linear_lambda_C', NaN, ...
    'mpl_decomposable_pilot_group_lambda_T', NaN, ...
    'mpl_decomposable_pilot_group_lambda_C', NaN, ...
    'mpl_decomposable_pilot_objective_initial', NaN, ...
    'mpl_decomposable_pilot_objective_final', NaN, ...
    'mpl_decomposable_pilot_objective_min_delta', NaN, ...
    'mpl_decomposable_pilot_objective_rejected', false, ...
    'mpl_decomposable_pilot_stationarity_beta', NaN, ...
    'mpl_decomposable_pilot_stationarity_phi', NaN, ...
    'mpl_decomposable_pilot_stationarity_group_T', NaN, ...
    'mpl_decomposable_pilot_stationarity_group_C', NaN, ...
    'mpl_group_iht_enabled', false, ...
    'mpl_group_iht_warmup_iterations', NaN, ...
    'mpl_group_iht_continuation_iterations', NaN, ...
    'mpl_group_iht_start_fraction', NaN, ...
    'mpl_group_iht_threshold_T', NaN, 'mpl_group_iht_threshold_C', NaN, ...
    'mpl_group_iht_lambda_T', NaN, 'mpl_group_iht_lambda_C', NaN, ...
    'mpl_group_iht_objective_initial', NaN, 'mpl_group_iht_objective_final', NaN, ...
    'mpl_group_iht_objective_min_delta', NaN, ...
    'mpl_group_iht_objective_rejected', false, ...
    'mpl_group_iht_support_changes', 0, ...
    'mpl_group_iht_support_stable_iterations', 0, ...
    'mpl_step_strategy', "fixed", ...
    'mpl_step_backtracks_beta', 0, 'mpl_step_backtracks_phi', 0, ...
    'mpl_step_backtracks_spline_T', 0, 'mpl_step_backtracks_spline_C', 0, ...
    'mpl_step_fallbacks_beta', 0, 'mpl_step_fallbacks_phi', 0, ...
    'mpl_step_fallbacks_spline_T', 0, 'mpl_step_fallbacks_spline_C', 0, ...
    'mpl_step_final_L_beta', NaN, 'mpl_step_final_L_phi', NaN, ...
    'mpl_step_final_L_spline_T', NaN, 'mpl_step_final_L_spline_C', NaN, ...
    'mpl_outer_iterations', 0, ...
    'mpl_selection_gate_pass', false, 'mpl_selection_gate_maximum', NaN, ...
    'mpl_selection_status_code', NaN, 'mpl_step_doublings', 0, ...
    'mpl_target_iterations', NaN, ...
    'mpl_oracle_refit_coefficient_score_inf_raw', NaN);
end


function lambda0 = mic_penalty_lambda(n_events, n, mic_penalty_scale, normalization)
% MIC/SELO complexity weight. The default "sum" scale reproduces the original
% simulation code; "mean" is retained for likelihood-scale diagnostics.
n_events = max(n_events, 2);
base_lambda = log(n_events);
if strcmpi(string(normalization), "mean")
    base_lambda = base_lambda / n;
end
lambda0 = mic_penalty_scale * base_lambda;
end


function [w_beta,w_phi] = linear_l1_weights(beta0,phi0,penalty_type,gamma,epsilon,weight_cap)
penalty_type = lower(string(penalty_type));
switch penalty_type
    case "none"
        w_beta = zeros(size(beta0));
        w_phi = zeros(size(phi0));
    case "lasso"
        w_beta = ones(size(beta0));
        w_phi = ones(size(phi0));
    case "adaptive"
        w_beta = min((abs(beta0) + epsilon).^(-gamma), weight_cap);
        w_phi = min((abs(phi0) + epsilon).^(-gamma), weight_cap);
    otherwise
        error('Unsupported LinearL1Type: %s', penalty_type);
end
end


function [mic_weight,l1_coefficient] = linear_penalty_local_terms(coef,a,lambda0,mic_form,kink_coefficient,external_l1_scale,external_l1_weight)
% MM/LLA coefficients for the supported linear penalties. SCAD and MCP use
% the same derivative at zero, lambda0*kink_coefficient, as kinked MIC.
mic_form = lower(string(mic_form));
switch mic_form
    case "smooth"
        mic_weight = smooth_mic_curvature_weight(coef, a, lambda0, 1e-6);
        internal_l1 = zeros(size(coef));
    case "kinked"
        inner = a*coef.^2 + kink_coefficient*abs(coef);
        exp_neg = exp(-2*min(inner,350));
        sech_sq = 4*exp_neg./(1+exp_neg).^2;
        outer_weight = lambda0.*sech_sq;
        mic_weight = a*outer_weight + 1e-6;
        internal_l1 = kink_coefficient*outer_weight;
    case "scad"
        if kink_coefficient <= 0
            error('PLAPH:LinearPenaltySlope', ...
                'SCAD requires a positive local-slope multiplier.');
        end
        penalty_lambda = lambda0*kink_coefficient;
        scad_shape = 3.7;
        abs_coef = abs(coef);
        internal_l1 = penalty_lambda*ones(size(coef));
        middle = abs_coef > penalty_lambda & ...
            abs_coef <= scad_shape*penalty_lambda;
        internal_l1(middle) = ...
            (scad_shape*penalty_lambda-abs_coef(middle))/(scad_shape-1);
        internal_l1(abs_coef > scad_shape*penalty_lambda) = 0;
        mic_weight = zeros(size(coef));
    case "mcp"
        if kink_coefficient <= 0
            error('PLAPH:LinearPenaltySlope', ...
                'MCP requires a positive local-slope multiplier.');
        end
        penalty_lambda = lambda0*kink_coefficient;
        mcp_concavity = 3;
        internal_l1 = max(penalty_lambda-abs(coef)/mcp_concavity,0);
        mic_weight = zeros(size(coef));
    otherwise
        error('Unsupported LinearMICForm: %s', mic_form);
end
l1_coefficient = internal_l1 + external_l1_scale.*external_l1_weight;
end


function [mic_weight,l1_coefficient,ridge_scale] = ...
    linear_penalty_terms_for_iteration(coef,a,lambda0,mic_form, ...
    kink_coefficient,external_l1_scale,external_l1_weight, ...
    target_ridge_scale,pilot_active,pilot_lambda)
if pilot_active
    mic_weight = zeros(size(coef));
    if isscalar(pilot_lambda)
        l1_coefficient = pilot_lambda*ones(size(coef));
    else
        l1_coefficient = pilot_lambda(:);
        if numel(l1_coefficient) ~= numel(coef)
            error('PLAPH:PilotLinearWeightSize', ...
                'Pilot linear penalty weights must match the coefficient block.');
        end
    end
    ridge_scale = 0;
else
    [mic_weight,l1_coefficient] = linear_penalty_local_terms( ...
        coef,a,lambda0,mic_form,kink_coefficient, ...
        external_l1_scale,external_l1_weight);
    ridge_scale = target_ridge_scale;
end
end


function weight = smooth_mic_curvature_weight(coef, a, lambda0, epsilon)
% Evaluate 4*a*lambda*exp(x)/(exp(x)+1)^2 without creating Inf/Inf.
% The original arithmetic is retained for x<=350 so ordinary runs are unchanged.
exp_arg = 2*a*coef.^2;
weight = zeros(size(coef));
regular = isfinite(exp_arg) & exp_arg <= 350;
exp_term = exp(exp_arg(regular));
weight(regular) = (4*a*lambda0*exp_term + epsilon)./ ...
    ((exp_term+1).^2 + epsilon);
weight(~isfinite(weight)) = 0;
end


function coef_new = hybrid_linear_prox(z,mic_weight,step_size,l1_coefficient,ridge_lambda)
% Exact proximal update for all nonsmooth absolute-value terms.
threshold = l1_coefficient/step_size;
shrunk = sign(z).*max(abs(z)-threshold,0);
denom = 1 + (2*mic_weight + ridge_lambda)/step_size;
coef_new = shrunk./denom;
end


function [coef_new,L_new,n_backtracks,used_fallback] = ...
    adaptive_linear_prox_update(coef,gradient,mic_weight,l1_coefficient, ...
    ridge_lambda,base_L,current_L,step_options,loss_function)
% Safeguarded proximal-gradient update using the smooth-loss majorizer.
strategy = lower(string(step_options.strategy));
if strategy == "auto"
    L_trial = current_L;
else
    L_trial = base_L;
end
L_min = base_L*step_options.min_scale;
L_max = base_L*step_options.max_scale;
L_trial = min(max(L_trial,L_min),L_max);
base_loss = loss_function(coef);
n_backtracks = 0;
used_fallback = false;
accepted = false;
coef_new = coef;
for attempt = 0:step_options.max_backtracks
    candidate = hybrid_linear_prox(coef-gradient/L_trial,mic_weight, ...
        L_trial,l1_coefficient,ridge_lambda);
    candidate(~isfinite(candidate)) = 0;
    direction = candidate-coef;
    trial_loss = loss_function(candidate);
    upper_bound = base_loss + gradient'*direction + ...
        0.5*L_trial*sum(direction.^2);
    tolerance = 1e-10*(1+abs(base_loss));
    if isfinite(trial_loss) && trial_loss <= upper_bound+tolerance
        coef_new = candidate;
        accepted = true;
        break;
    end
    n_backtracks = n_backtracks+1;
    L_trial = min(L_trial*step_options.backtracking_factor,L_max);
end
if ~accepted
    used_fallback = true;
    L_trial = base_L;
    coef_new = hybrid_linear_prox(coef-gradient/base_L,mic_weight, ...
        base_L,l1_coefficient,ridge_lambda);
    coef_new(~isfinite(coef_new)) = 0;
end
L_new = L_trial;
end


function [coef_new,L_new,n_backtracks,used_fallback] = ...
    adaptive_group_iht_update(metric,coef,gradient,num_basis,q, ...
    base_L,current_L,stage_lambda,step_options,loss_function)
% Backtracking IHT with a fixed statistical penalty and adaptive threshold.
strategy = lower(string(step_options.strategy));
% Group-IHT support is discontinuous. Auto mode carries its last accepted
% curvature but does not proactively enlarge the step; it only backtracks.
if strategy == "auto"
    L_trial = current_L;
else
    L_trial = base_L;
end
L_min = base_L*step_options.min_scale;
L_max = base_L*step_options.max_scale;
L_trial = min(max(L_trial,L_min),L_max);
base_loss = loss_function(coef);
n_backtracks = 0;
used_fallback = false;
accepted = false;
coef_new = coef;
for attempt = 0:step_options.max_backtracks
    threshold = sqrt(max(2*stage_lambda/L_trial,0));
    candidate = apply_spline_group_iht_metric( ...
        metric,coef,gradient,num_basis,q,L_trial,threshold);
    direction = candidate-coef;
    metric_norm_squared = spline_metric_norm_squared( ...
        metric,direction,num_basis,q);
    trial_loss = loss_function(candidate);
    upper_bound = base_loss + gradient'*direction + ...
        0.5*L_trial*metric_norm_squared;
    tolerance = 1e-10*(1+abs(base_loss));
    if isfinite(trial_loss) && trial_loss <= upper_bound+tolerance
        coef_new = candidate;
        accepted = true;
        break;
    end
    n_backtracks = n_backtracks+1;
    L_trial = min(L_trial*step_options.backtracking_factor,L_max);
end
if ~accepted
    used_fallback = true;
    L_trial = base_L;
    threshold = sqrt(max(2*stage_lambda/base_L,0));
    coef_new = apply_spline_group_iht_metric( ...
        metric,coef,gradient,num_basis,q,base_L,threshold);
end
L_new = L_trial;
end


function value = spline_metric_norm_squared(metric,direction,num_basis,q)
value = 0;
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block = direction(idx);
    value = value + max(real(block'*metric.gram{j}*block),0);
end
end


function value = max_spline_metric_effect_norm(metric,direction,num_basis,q)
value = 0;
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block = direction(idx);
    value = max(value,sqrt(max(real(block'*metric.gram{j}*block),0)));
end
end


function nfp = count_fp(est, active, threshold_linear)
selected = find(abs(est) > threshold_linear);
nfp = numel(setdiff(selected(:), active(:)));
end


function nfn = count_fn(est, active, threshold_linear)
selected = find(abs(est) > threshold_linear);
nfn = numel(setdiff(active(:), selected(:)));
end


function nfp = support_false_positives(selected,active)
truth = false(numel(selected),1);
truth(active) = true;
nfp = sum(logical(selected(:)) & ~truth);
end


function nfn = support_false_negatives(selected,active)
truth = false(numel(selected),1);
truth(active) = true;
nfn = sum(~logical(selected(:)) & truth);
end


function nsel = count_selected_groups(coef, num_basis, q, threshold)
nsel = 0;
strict_zero = threshold <= 0;
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    group_norm = norm(coef(idx), 2);
    if (strict_zero && group_norm > 0) || (~strict_zero && group_norm >= threshold)
        nsel = nsel + 1;
    end
end
end


function support = spline_group_support(coef, num_basis, q)
support = false(q, 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    support(j) = norm(coef(idx), 2) > 1e-10;
end
end


function coef = apply_spline_group_threshold(coef, num_basis, q, threshold)
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    if norm(coef(idx), 2) < threshold
        coef(idx) = 0;
    end
end
end


function nfp = count_group_fp(coef, num_basis, q, active_groups, threshold)
norms = spline_group_norms(coef, num_basis, q);
selected = find(group_selected(norms, threshold));
active_groups = active_groups(active_groups <= q);
nfp = numel(setdiff(selected(:), active_groups(:)));
end


function nfn = count_group_fn(coef, num_basis, q, active_groups, threshold)
norms = spline_group_norms(coef, num_basis, q);
selected = find(group_selected(norms, threshold));
active_groups = active_groups(active_groups <= q);
nfn = numel(setdiff(active_groups(:), selected(:)));
end


function nsel = count_selected_groups_metric(B, coef, num_basis, q, threshold, metric)
scores = spline_group_scores(B, coef, num_basis, q, metric);
nsel = sum(group_selected(scores, threshold));
end


function coef = apply_spline_group_threshold_metric(B, coef, num_basis, q, threshold, metric)
if threshold <= 0
    return;
end
scores = spline_group_scores(B, coef, num_basis, q, metric);
for j = 1:q
    if scores(j) < threshold
        idx = ((j-1)*num_basis+1):(j*num_basis);
        coef(idx) = 0;
    end
end
end


function nfp = count_group_fp_metric(B, coef, num_basis, q, active_groups, threshold, metric)
scores = spline_group_scores(B, coef, num_basis, q, metric);
selected = find(group_selected(scores, threshold));
active_groups = active_groups(active_groups <= q);
nfp = numel(setdiff(selected(:), active_groups(:)));
end


function nfn = count_group_fn_metric(B, coef, num_basis, q, active_groups, threshold, metric)
scores = spline_group_scores(B, coef, num_basis, q, metric);
selected = find(group_selected(scores, threshold));
active_groups = active_groups(active_groups <= q);
nfn = numel(setdiff(active_groups(:), selected(:)));
end


function selected = group_selected(scores, threshold)
if threshold <= 0
    selected = scores > 0;
else
    selected = scores >= threshold;
end
end


function scores = spline_group_scores(B, coef, num_basis, q, metric)
metric = lower(string(metric));
switch metric
    case "coef"
        scores = spline_group_norms(coef, num_basis, q);
    case "effect"
        scores = spline_effect_norms(B, coef, num_basis, q);
    otherwise
        error('Unsupported spline selection metric: %s', metric);
end
end


function coef = apply_spline_group_lasso_prox(coef, num_basis, q, penalty_over_step)
if penalty_over_step <= 0
    return;
end
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    group_norm = norm(coef(idx), 2);
    if group_norm <= penalty_over_step
        coef(idx) = 0;
    else
        coef(idx) = (1 - penalty_over_step/group_norm).*coef(idx);
    end
end
end


function coef = apply_spline_group_penalty_step(coef, num_basis, q, lambda, step_size, penalty_type, sharpness, use_continuation, sharpness_path, mic_form, kink_coefficient, mcp_concavity, iht_threshold)
if nargin < 8 || isempty(use_continuation)
    use_continuation = false;
end
if nargin < 9
    sharpness_path = [];
end
if nargin < 10 || isempty(mic_form)
    mic_form = "smooth";
end
if nargin < 11 || isempty(kink_coefficient)
    kink_coefficient = 0;
end
if nargin < 12 || isempty(mcp_concavity)
    mcp_concavity = 3;
end
if nargin < 13 || isempty(iht_threshold)
    iht_threshold = 0.1;
end
penalty_type = lower(string(penalty_type));
mic_form = lower(string(mic_form));
if penalty_type == "group_iht"
    coef = apply_spline_group_iht(coef, num_basis, q, iht_threshold);
    return;
end
if lambda <= 0
    return;
end
switch penalty_type
    case "group_lasso"
        coef = apply_spline_group_lasso_prox(coef, num_basis, q, lambda/step_size);
    case "group_mic"
        if use_continuation
            coef = apply_spline_group_mic_continuation(coef, num_basis, q, lambda, sharpness, step_size, sharpness_path, mic_form, kink_coefficient);
        else
            coef = apply_spline_group_mic_update(coef, num_basis, q, lambda, sharpness, step_size, mic_form, kink_coefficient);
        end
    case "group_mcp"
        coef = apply_spline_group_mcp_prox(coef, num_basis, q, lambda, step_size, mcp_concavity);
    otherwise
        error('Unsupported spline group penalty type: %s', penalty_type);
end
end


function coef = apply_spline_group_iht(coef, num_basis, q, threshold)
% Exact proximal map for lambda_A*||c||_{2,0}, with
% threshold=sqrt(2*lambda_A/step_size) supplied by the caller.
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    if norm(coef(idx), 2) <= threshold
        coef(idx) = 0;
    end
end
end


function coef_new = apply_spline_group_iht_metric(metric, coef, gradient, num_basis, q, step_size, threshold)
% Equivalent implementation of IHT in c_j=G_j^(1/2)b_j coordinates without
% explicitly whitening the design matrix. The pseudoinverse removes the
% centered B-spline null direction and the threshold is the fitted-effect RMS.
coef_new = zeros(size(coef));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    gram = metric.gram{j};
    candidate = coef(idx) - metric.gram_pinv{j}*gradient(idx)/step_size;
    candidate(~isfinite(candidate)) = 0;
    effect_norm = sqrt(max(real(candidate' * gram * candidate), 0));
    if effect_norm > threshold
        coef_new(idx) = candidate;
    end
end
end


function effect_norms = spline_group_candidate_effect_norms( ...
    metric,coef,gradient,num_basis,q,step_size)
effect_norms = zeros(q,1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    candidate = coef(idx)-metric.gram_pinv{j}*gradient(idx)/step_size;
    candidate(~isfinite(candidate)) = 0;
    effect_norms(j) = sqrt(max(real( ...
        candidate'*metric.gram{j}*candidate),0));
end
end


function coef_new = apply_spline_group_lasso_metric(metric, coef, gradient, ...
    num_basis, q, step_size, lambda)
% Proximal map for lambda*sqrt(b_j'*G_j*b_j) in the fitted-effect metric.
if ~(isscalar(lambda) || numel(lambda) == q)
    error('PLAPH:PilotGroupWeightSize', ...
        'Group penalty weights must be scalar or have one value per group.');
end
coef_new = zeros(size(coef));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    gram = metric.gram{j};
    candidate = coef(idx) - metric.gram_pinv{j}*gradient(idx)/step_size;
    candidate(~isfinite(candidate)) = 0;
    effect_norm = sqrt(max(real(candidate' * gram * candidate), 0));
    lambda_j = lambda(min(j,numel(lambda)));
    if effect_norm > lambda_j/step_size
        coef_new(idx) = (1-lambda_j/(step_size*effect_norm))*candidate;
    end
end
end


function metric = prepare_spline_group_metric(B, num_basis, q)
% The spline design is fixed within a fit, so cache every Gram matrix and
% pseudoinverse instead of repeating the eigendecomposition in each IHT step.
metric.gram = cell(q, 1);
metric.gram_pinv = cell(q, 1);
n = size(B, 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block_design = B(:, idx);
    gram = (block_design' * block_design) / n;
    gram = (gram + gram') / 2;
    [vectors, values] = eig(gram);
    eigenvalues = real(diag(values));
    largest = max([eigenvalues; 0]);
    tolerance = max(1e-10, largest * num_basis * eps);
    inverse_values = zeros(size(eigenvalues));
    estimable = eigenvalues > tolerance;
    inverse_values(estimable) = 1 ./ eigenvalues(estimable);
    metric.gram{j} = gram;
    metric.gram_pinv{j} = real(vectors * diag(inverse_values) * vectors');
end
end


function coef = apply_spline_group_mcp_prox(coef, num_basis, q, lambda, step_size, concavity)
if lambda <= 0
    return;
end
if concavity <= 1
    error('The group-MCP concavity parameter must exceed one.');
end

% The block MCP proximal map is unique when the quadratic curvature exceeds
% 1/concavity. The fitted algorithms use much larger step curvatures, while
% this guard keeps the update well defined for diagnostic settings.
step_size = max(step_size, 1/concavity + 1e-8);
zero_cutoff = lambda/step_size;
flat_cutoff = concavity*lambda;
middle_denominator = 1 - 1/(concavity*step_size);

for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    block = coef(idx);
    block_norm = norm(block, 2);
    if block_norm <= zero_cutoff
        coef(idx) = 0;
    elseif block_norm <= flat_cutoff
        shrinkage = (1 - lambda/(step_size*block_norm))/middle_denominator;
        coef(idx) = shrinkage.*block;
    end
end
end


function coef = apply_spline_group_mic_continuation(coef, num_basis, q, lambda, final_sharpness, step_size, sharpness_path, mic_form, kink_coefficient)
if isempty(sharpness_path)
    sharpness_path = unique([2 4 8 16 32 final_sharpness], "stable");
end
sharpness_path = sharpness_path(:)';
sharpness_path = sharpness_path(isfinite(sharpness_path) & sharpness_path > 0);
if isempty(sharpness_path)
    sharpness_path = final_sharpness;
end
lambda_stage = lambda/numel(sharpness_path);
for sharpness = sharpness_path
    coef = apply_spline_group_mic_update(coef, num_basis, q, lambda_stage, sharpness, step_size, mic_form, kink_coefficient);
end
end


function coef = apply_spline_group_mic_update(coef, num_basis, q, lambda, sharpness, step_size, mic_form, kink_coefficient)
if lambda <= 0
    return;
end
mic_form = lower(string(mic_form));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    group_coef = coef(idx);
    group_norm = norm(group_coef, 2);
    if mic_form == "kinked"
        penalty_argument = sharpness*group_norm^2 + kink_coefficient*group_norm;
    else
        penalty_argument = sharpness*group_norm^2;
    end
    exp_arg = min(2*penalty_argument, 700);
    exp_term = exp(exp_arg);
    if mic_form == "kinked"
        outer_weight = (4*lambda*exp_term + 1e-8)/((exp_term + 1)^2 + 1e-8);
        ridge_weight = sharpness*outer_weight;
        threshold = kink_coefficient*outer_weight/step_size;
        if group_norm <= threshold
            coef(idx) = 0;
        else
            coef(idx) = ((1-threshold/group_norm)/(1+2*ridge_weight/step_size)).*group_coef;
        end
    elseif mic_form == "smooth"
        % Preserve the archived smooth group-MIC update exactly.
        ridge_weight = (4*sharpness*lambda*exp_term + 1e-8)/((exp_term + 1)^2 + 1e-8);
        coef(idx) = group_coef/(1+2*ridge_weight/step_size);
    else
        error('Unsupported spline group MIC form: %s', mic_form);
    end
end
end


function lambda = spline_group_penalty_lambda(n_events, n, scale, normalization)
if scale <= 0
    lambda = 0;
    return;
end
n_events = max(n_events, 2);
if lower(string(normalization)) == "mean"
    lambda = scale * log(n_events) / n;
else
    lambda = scale * log(n_events);
end
end


function val = min_active_group_norm(coef, num_basis, q, active_groups)
norms = spline_group_norms(coef, num_basis, q);
active_groups = active_groups(active_groups <= q);
if isempty(active_groups)
    val = NaN;
else
    val = min(norms(active_groups));
end
end


function val = max_noise_group_norm(coef, num_basis, q, active_groups)
norms = spline_group_norms(coef, num_basis, q);
noise_groups = setdiff(1:q, active_groups(active_groups <= q));
if isempty(noise_groups)
    val = NaN;
else
    val = max(norms(noise_groups));
end
end


function norms = spline_group_norms(coef, num_basis, q)
norms = zeros(q,1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    norms(j) = norm(coef(idx), 2);
end
end


function val = min_active_effect_norm(B, coef, num_basis, q, active_groups)
norms = spline_effect_norms(B, coef, num_basis, q);
active_groups = active_groups(active_groups <= q);
if isempty(active_groups)
    val = NaN;
else
    val = min(norms(active_groups));
end
end


function val = max_noise_effect_norm(B, coef, num_basis, q, active_groups)
norms = spline_effect_norms(B, coef, num_basis, q);
noise_groups = setdiff(1:q, active_groups(active_groups <= q));
if isempty(noise_groups)
    val = NaN;
else
    val = max(norms(noise_groups));
end
end


function val = noise_effect_quantile(B, coef, num_basis, q, active_groups, prob)
norms = spline_effect_norms(B, coef, num_basis, q);
noise_groups = setdiff(1:q, active_groups(active_groups <= q));
if isempty(noise_groups)
    val = NaN;
else
    val = quantile(norms(noise_groups), prob);
end
end


function norms = spline_effect_norms(B, coef, num_basis, q)
norms = zeros(q,1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    fitted_effect = B(:,idx)*coef(idx);
    norms(j) = sqrt(mean(fitted_effect.^2));
end
end


function norms = active_effect_norms(B, coef, num_basis, q, active_groups)
all_norms = spline_effect_norms(B, coef, num_basis, q);
active_groups = active_groups(active_groups <= q);
norms = all_norms(active_groups);
end


function rep = assign_active_effect_fields(rep, prefix, norms)
for j = 1:min(4, numel(norms))
    rep.(sprintf('%s_g%d', prefix, j)) = norms(j);
end
end


function fits = initialize_function_fit_storage(grid_size,active_additive,R, ...
    additive_linear_coefficient,additive_signal_scale,dgp_type, ...
    censoring_function_design)
n_active = numel(active_additive);
grid_size = round(grid_size);
fits = struct();
fits.components = active_additive(:)';
fits.grid = NaN(grid_size, n_active);
fits.true = NaN(grid_size, n_active);
fits.true_T = NaN(grid_size,n_active);
fits.true_C = NaN(grid_size,n_active);
fits.mpl_T = NaN(grid_size, n_active, R);
fits.mpl_C = NaN(grid_size, n_active, R);
fits.pl = NaN(grid_size, n_active, R);
fits.oracle = NaN(grid_size, n_active, R);
fits.note = "Curves are centered over the plotting grid to compare identifiable shapes.";

for j = 1:n_active
    component = active_additive(j);
    fits.grid(:,j) = function_plot_grid(component,grid_size,dgp_type, ...
        censoring_function_design);
    fits.true_T(:,j) = center_curve(true_margin_component_values( ...
        dgp_type,"T",component,fits.grid(:,j), ...
        additive_linear_coefficient,additive_signal_scale, ...
        censoring_function_design));
    fits.true_C(:,j) = center_curve(true_margin_component_values( ...
        dgp_type,"C",component,fits.grid(:,j), ...
        additive_linear_coefficient,additive_signal_scale, ...
        censoring_function_design));
end
fits.true = fits.true_T;
end


function grid = function_plot_grid(component,grid_size,dgp_type, ...
        censoring_function_design)
if nargin < 4
    censoring_function_design = "current";
end
if lower(string(dgp_type)) == "discriminating"
    bounded_quadratic = any(lower(string(censoring_function_design)) == ...
        ["shared_opposite_bounded" "shared_balanced_bounded"]) && ...
        component == 2;
    gaussian_component = component <= 2 && ~bounded_quadratic;
else
    gaussian_component = component <= 2;
end
if gaussian_component
    grid = linspace(-2, 2, grid_size)';
else
    grid = linspace(0, 1, grid_size)';
end
end


function fits = store_function_fit_curves(fits, field_name, iter, B, X, coef, num_basis, active_additive)
if isempty(coef)
    return;
end

for j = 1:numel(active_additive)
    component = active_additive(j);
    idx = ((component-1)*num_basis+1):(component*num_basis);
    if idx(end) > numel(coef) || idx(end) > size(B,2)
        continue;
    end
    observed_fit = B(:,idx)*coef(idx);
    curve = curve_from_observed_values(X(:,component), observed_fit, fits.grid(:,j));
    fits.(field_name)(:,j,iter) = center_curve(curve);
end
end


function curve = curve_from_observed_values(x_observed, y_observed, grid)
valid = isfinite(x_observed) & isfinite(y_observed);
x_observed = x_observed(valid);
y_observed = y_observed(valid);
if numel(x_observed) < 2
    curve = NaN(size(grid));
    return;
end

[x_sorted, order] = sort(x_observed);
y_sorted = y_observed(order);
[x_unique, ~, groups] = unique(x_sorted);
y_unique = accumarray(groups, y_sorted, [], @mean);
if numel(x_unique) < 2
    curve = NaN(size(grid));
else
    curve = interp1(x_unique, y_unique, grid, "linear", "extrap");
end
end


function centered = center_curve(values)
centered = values - mean(values, 'omitnan');
end


function values = true_margin_component_values(dgp_type,margin,component,x, ...
    additive_linear_coefficient,additive_signal_scale,censoring_function_design)
if nargin < 7
    censoring_function_design = "current";
end
dgp_type = lower(string(dgp_type));
margin = upper(string(margin));
bounded_design = dgp_type == "discriminating" && ...
    any(lower(string(censoring_function_design)) == ...
    ["shared_opposite_bounded" "shared_balanced_bounded"]);
if bounded_design && component == 2
    values = sqrt(5/2).*(6*x.^2-6*x+1);
    if margin == "C" && lower(string(censoring_function_design)) == ...
            "shared_opposite_bounded"
        values = -values;
    end
elseif dgp_type == "discriminating" && margin == "C" && ...
        any(lower(string(censoring_function_design)) == ...
        ["shared_opposite" "shared_opposite_bounded" ...
        "shared_balanced_bounded"])
    switch component
        case 1
            values = -additive_linear_coefficient*x;
        case 2
            values = -(0.5*x.^2+0.5*x-0.5);
        case 3
            values = -sin(2*pi*x);
        case 4
            values = -cos(2*pi*x);
        otherwise
            values = zeros(size(x));
    end
elseif dgp_type == "discriminating" && margin == "C"
    switch component
        case 1
            values = -additive_linear_coefficient*x;
        case 3
            values = cos(2*pi*x);
        case 5
            values = sin(4*pi*x);
        case 6
            % Centered monotone quadratic: nonlinear, but with nonvanishing
            % marginal screening utility under the prespecified null score.
            values = 2*x.^2-2/3;
        otherwise
            values = zeros(size(x));
    end
else
    switch component
        case 1
            values = additive_linear_coefficient*x;
        case 2
            values = 0.5*x.^2+0.5*x-0.5;
        case 3
            values = sin(2*pi*x);
        case 4
            values = cos(2*pi*x);
        otherwise
            values = zeros(size(x));
    end
end
values = additive_signal_scale.*values;
end


function val = additive_imse(B,X,coef,num_basis,q,active_groups, ...
    additive_linear_coefficient,additive_signal_scale,dgp_type,margin,compact_coef, ...
    censoring_function_design,grid_size)
if nargin < 12
    censoring_function_design = "current";
end
if nargin < 13 || isempty(grid_size)
    grid_size = 101;
end
active_groups = active_groups(active_groups <= q);
n_active = numel(active_groups);
if n_active == 0
    val = NaN;
    return;
end
fit = NaN(round(grid_size),n_active);
truth = NaN(round(grid_size),n_active);
for j = 1:n_active
    component = active_groups(j);
    design_idx = ((component-1)*num_basis+1):(component*num_basis);
    if compact_coef
        coef_idx = ((j-1)*num_basis+1):(j*num_basis);
    else
        coef_idx = design_idx;
    end
    observed_fit = B(:,design_idx)*coef(coef_idx);
    grid = function_plot_grid(component,grid_size,dgp_type, ...
        censoring_function_design);
    fit(:,j) = curve_from_observed_values( ...
        X(:,component),observed_fit,grid);
    truth(:,j) = true_margin_component_values(dgp_type,margin,component, ...
        grid,additive_linear_coefficient,additive_signal_scale, ...
        censoring_function_design);
    % Additive PH components are identifiable only up to a constant.  The
    % spline blocks are sample-centered, so compare the same centered shapes.
    % v2: normalized trapezoidal weights on the (equally spaced) grid are used
    % both for centering and for integration, as stated in the manuscript.
    w = trapezoid_weights(size(fit,1));
    fit(:,j) = fit(:,j)-sum(w.*fit(:,j),'omitnan');
    truth(:,j) = truth(:,j)-sum(w.*truth(:,j),'omitnan');
end
w = trapezoid_weights(size(fit,1));
val = mean(sum(w.*(fit-truth).^2,1,'omitnan'));
end


function summary = summarize_replicates(T)
ok = T.ok;
summary = table();
summary.R = height(T);
summary.n_ok = sum(ok);
summary.fail_rate = 1 - mean(ok);
summary.mpl_baseline_armijo_fail_rate = mean(T.mpl_baseline_armijo_failed, 'omitnan');
summary.mpl_armijo_reject_rate = mean(T.mpl_armijo_rejected, 'omitnan');
summary.mean_censor_rate = mean(T.censor_rate(ok), 'omitnan');
summary.mean_n_failure_events = mean(T.n_failure_events(ok), 'omitnan');
summary.mean_n_censoring_events = mean(T.n_censoring_events(ok), 'omitnan');
summary.mean_screen_linear_union_size = mean(T.screen_linear_union_size, 'omitnan');
summary.mean_screen_additive_union_size = mean(T.screen_additive_union_size, 'omitnan');
summary.screen_beta_sure_rate = mean(T.screen_beta_sure, 'omitnan');
summary.screen_phi_sure_rate = mean(T.screen_phi_sure, 'omitnan');
summary.screen_add_T_sure_rate = mean(T.screen_add_T_sure, 'omitnan');
summary.screen_add_C_sure_rate = mean(T.screen_add_C_sure, 'omitnan');
summary.screen_all_sure_rate = mean(T.screen_all_sure, 'omitnan');
summary.mean_mse_mpl_beta = mean(T.mse_mpl_beta(ok), 'omitnan');
summary.mean_mse_mpl_phi = mean(T.mse_mpl_phi(ok), 'omitnan');
summary.mean_mse_pl = mean(T.mse_pl(ok), 'omitnan');
summary.mean_mse_oracle = mean(T.mse_oracle(ok), 'omitnan');
summary.mean_imse_mpl_add_T = mean(T.imse_mpl_add_T(ok), 'omitnan');
summary.mean_imse_mpl_add_C = mean(T.imse_mpl_add_C(ok), 'omitnan');
summary.mean_imse_pl_add = mean(T.imse_pl_add(ok), 'omitnan');
summary.mean_imse_oracle_add = mean(T.imse_oracle_add(ok), 'omitnan');
summary.mean_tau_hat = mean(T.tau_hat(ok), 'omitnan');
summary.mean_beta_fp = mean(T.beta_fp(ok), 'omitnan');
summary.mean_beta_fn = mean(T.beta_fn(ok), 'omitnan');
summary.mean_phi_fp = mean(T.phi_fp(ok), 'omitnan');
summary.mean_phi_fn = mean(T.phi_fn(ok), 'omitnan');
summary.mean_mpl_beta_selected = mean(T.mpl_beta_selected(ok), 'omitnan');
summary.mean_mpl_phi_selected = mean(T.mpl_phi_selected(ok), 'omitnan');
summary.mean_mpl_add_T_fp = mean(T.mpl_add_T_fp(ok), 'omitnan');
summary.mean_mpl_add_T_fn = mean(T.mpl_add_T_fn(ok), 'omitnan');
summary.mean_mpl_add_C_fp = mean(T.mpl_add_C_fp(ok), 'omitnan');
summary.mean_mpl_add_C_fn = mean(T.mpl_add_C_fn(ok), 'omitnan');
summary.mean_pl_add_fp = mean(T.pl_add_fp(ok), 'omitnan');
summary.mean_pl_add_fn = mean(T.pl_add_fn(ok), 'omitnan');
summary.mean_mpl_eff_T_active_min = mean(T.mpl_eff_T_active_min(ok), 'omitnan');
summary.mean_mpl_eff_T_noise_q90 = mean(T.mpl_eff_T_noise_q90(ok), 'omitnan');
summary.mean_mpl_eff_T_noise_q95 = mean(T.mpl_eff_T_noise_q95(ok), 'omitnan');
summary.mean_mpl_eff_T_noise_max = mean(T.mpl_eff_T_noise_max(ok), 'omitnan');
summary.mean_mpl_eff_C_active_min = mean(T.mpl_eff_C_active_min(ok), 'omitnan');
summary.mean_mpl_eff_C_noise_q90 = mean(T.mpl_eff_C_noise_q90(ok), 'omitnan');
summary.mean_mpl_eff_C_noise_q95 = mean(T.mpl_eff_C_noise_q95(ok), 'omitnan');
summary.mean_mpl_eff_C_noise_max = mean(T.mpl_eff_C_noise_max(ok), 'omitnan');
summary.mean_pl_eff_active_min = mean(T.pl_eff_active_min(ok), 'omitnan');
summary.mean_pl_eff_noise_q90 = mean(T.pl_eff_noise_q90(ok), 'omitnan');
summary.mean_pl_eff_noise_q95 = mean(T.pl_eff_noise_q95(ok), 'omitnan');
summary.mean_pl_eff_noise_max = mean(T.pl_eff_noise_max(ok), 'omitnan');
end

function  [T,Z,X,B,status,T1,Z1,X1,B1,status1,num_basis] = survival_data(n,q,Beta,Phi,mu,sigma,copula_family,tau,lamt,iter,num_internal_knots,spline_basis_scaling,additive_linear_coefficient,additive_signal_scale,dgp_type,additive_rho,censoring_scale,censoring_function_design,build_spline_basis)
if nargin < 12 || isempty(spline_basis_scaling)
    spline_basis_scaling = "center";
end
if nargin < 19
    build_spline_basis = true;
end
spline_basis_scaling = lower(string(spline_basis_scaling));
dgp_type = lower(string(dgp_type));
%% Generating the survival data
% p = length(Beta);
ZZ = mvnrnd(mu,sigma,n);  % n-by-p covariate matrix

if dgp_type == "discriminating"
    additive_index = 1:q;
    additive_sigma = additive_rho.^abs(additive_index-additive_index');
    latent_X = mvnrnd(zeros(q,1),additive_sigma,n);
    XX = latent_X;
    uniform_components = 3:min(6,q);
    XX(:,uniform_components) = normcdf(latent_X(:,uniform_components));
    if any(lower(string(censoring_function_design)) == ...
            ["shared_opposite_bounded" "shared_balanced_bounded"])
        XX(:,2) = normcdf(latent_X(:,2));
    end
    item_T_X = zeros(n,1);
    item_C_X = zeros(n,1);
    for component = [1 2 3 4]
        item_T_X = item_T_X+true_margin_component_values( ...
            dgp_type,"T",component,XX(:,component), ...
            additive_linear_coefficient,additive_signal_scale, ...
            censoring_function_design);
    end
    if any(lower(string(censoring_function_design)) == ...
            ["shared_opposite" "shared_opposite_bounded" ...
            "shared_balanced_bounded"])
        censoring_components = [1 2 3 4];
    else
        censoring_components = [1 3 5 6];
    end
    for component = censoring_components
        item_C_X = item_C_X+true_margin_component_values( ...
            dgp_type,"C",component,XX(:,component), ...
            additive_linear_coefficient,additive_signal_scale, ...
            censoring_function_design);
    end
else
    XX = randn(n,q);  % Gaussian covariates
    XX(:,3:4) = rand(n,2);  % Uniform(0,1) covariates
    item_T_X = additive_linear_coefficient*XX(:,1) + ...
        (0.5*XX(:,2).^2+0.5*XX(:,2)-0.5) +...
        sin(2*pi*XX(:,3)) + cos(2*pi*XX(:,4));
    item_T_X = additive_signal_scale.*item_T_X;
    item_C_X = item_T_X;
end

theta_data = copula_param_from_tau(copula_family,tau);  % copula dependence parameter for data generation
U = sample_copula(copula_family,theta_data,n);
V = plaph_survival_uniforms(copula_family,U);

%% failure time & censoring time
Failure_time = sqrt(-lamt^2*log(V(:,1))./safe_exp(ZZ*Beta+item_T_X));  % baseline failure hazard h0T = 2t/lamt^2
Censoring_time = -censoring_scale*log(V(:,2))./safe_exp(ZZ*Phi+item_C_X);

status_raw = (Failure_time <= Censoring_time);
TT = min(Failure_time,Censoring_time);     % observed survival time

if ~build_spline_basis
    [T,I] = sort(TT,'descend');
    Z = ZZ(I,:);
    X = XX(I,:);
    status = status_raw(I);
    B = zeros(n,0);
    [T1,I1] = sort(TT,'ascend');
    Z1 = ZZ(I1,:);
    X1 = XX(I1,:);
    status1 = status_raw(I1);
    B1 = zeros(n,0);
    num_basis = 0;
    return;
end



%% Computing the B-spline basis matrix

% v2: bspline_basis(j,n,...) takes the ORDER n (degree n-1). The legacy
% engine passed K=3 (quadratic) with order-4 knot multiplicity. The order is
% now explicit and consistent: 4 = cubic (default), 3 = quadratic.
order = plaph_spline_order();
K = order - 1; %#ok<NASGU> % degree

% num_internal_knots is supplied by the pilot-tuning step or by the default setting.
total_knots = num_internal_knots+2; % total number of knots
num_basis = num_internal_knots + order; % clamped B-splines: interior knots + order

BB = zeros(n, num_basis*q);
B_eval = zeros(n, num_basis);

%%


for i = 1:q
    % Construct B-spline basis functions for every nonparametric covariate.
    % This is important in high-dimensional simulations: if the basis matrix is
    % constructed only for the truly active covariates, false-positive selection
    % for inactive nonlinear covariates cannot be evaluated properly.
    knot_grid = quantile(XX(:,i), linspace(0,1,total_knots));
    knot_grid = unique(knot_grid, 'stable');

    if numel(knot_grid) < total_knots
        % Fall back to equally spaced knots if repeated quantiles occur.
        knot_grid = linspace(min(XX(:,i)), max(XX(:,i)), total_knots);
    end

    internal_knots = knot_grid(2:(end-1));
    left_boundary = knot_grid(1);
    right_boundary = knot_grid(end);

    if right_boundary <= left_boundary
        right_boundary = left_boundary + 1e-6;
    end

    knots = [left_boundary*ones(1, order), internal_knots, right_boundary*ones(1, order)];

    B_eval(:,:) = 0;
    for j = 1:num_basis
        [B_eval(:,j),~] = bspline_basis(j-1, order, knots, XX(:,i));
    end

    BB(:,((i-1)*num_basis+1):(i*num_basis)) = B_eval;
end

% for i = 1:q
%     for j = 1:num_basis
%         internal_knots = linspace(min(XX(:,i)),max(XX(:,i)), total_knots);
%         internal_knots = internal_knots(2:(total_knots-1));
%         knots = [min(XX(:,i))*ones(1, order), internal_knots, max(XX(:,i))*ones(1, order)];
%         % internal_knots = quantile(XX(:,i),linspace(0.001,0.99,num_internal_knots));
%         % knots = augknt([min(XX(:,i)), internal_knots, max(XX(:,i))], order);
%         [B_eval(:,j),~] = bspline_basis(j-1, K, knots, XX(:,i));
% 
%     end
% 
%     BB(:,((i-1)*num_basis+1):(i*num_basis) ) = B_eval;
% end

% Column-wise centering and scaling improves numerical stability of spline updates.
if spline_basis_scaling == "zscore"
    BB = BB - mean(BB,1);
    BB = BB ./ (std(BB,0,1) + 1e-8);
elseif spline_basis_scaling == "center"
    BB = BB - mean(BB,1);
elseif spline_basis_scaling == "whiten"
    BB = whiten_spline_basis_groups(BB, num_basis, q);
elseif spline_basis_scaling == "none"
    BB = BB;
else
    error('PLAPH:UnknownSplineBasisScaling', 'Unknown SplineBasisScaling: %s', spline_basis_scaling);
end
BB(~isfinite(BB)) = 0;

%% sorting the time and other covariates
[T,I] = sort(TT,'descend');
Z = ZZ(I,:);
X = XX(I,:);
status = status_raw(I);
B = BB(I,:);

[T1,I1] = sort(TT,'ascend');  % % sorting the time
Z1 = ZZ(I1,:);
X1 = XX(I1,:);
status1 = status_raw(I1);
B1 = BB(I1,:);

end



function result = plaph_run_real_data_impl(scenario, opts)
scenario = normalize_real_data_scenario(scenario);
estimate_tau = isfield(scenario,"profile_tau") && ...
    logical(scenario.profile_tau);
rng(opts.BaseSeed);

[clinical, patient_ids, Zraw, Xraw, z_features, x_features] = read_real_data_working_set(scenario);
if scenario.p > 0
    Zraw = Zraw(:, 1:min(scenario.p, size(Zraw, 2)));
    z_features = z_features(1:size(Zraw, 2));
end
if scenario.q > 0
    Xraw = Xraw(:, 1:min(scenario.q, size(Xraw, 2)));
    x_features = x_features(1:size(Xraw, 2));
end

TT = clinical.os_months;
Delta = logical(clinical.os_event);
% v2: optional unpenalized clinical adjustment (age, male, stage II/III/IV)
% entering both margins with penalty factor zero; complete cases only.
use_clinical = isfield(opts,"UnpenalizedClinical") && logical(opts.UnpenalizedClinical);
clinical_names = strings(0,1);
Cmat = zeros(numel(TT),0);
if use_clinical
    [Cmat,clinical_names] = clinical_adjustment_design(clinical);
end
valid = isfinite(TT) & TT > 0 & all(isfinite(Zraw), 2) & all(isfinite(Xraw), 2) & ...
    all(isfinite(Cmat), 2);
TT = TT(valid);
Delta = Delta(valid);
Zraw = Zraw(valid, :);
Xraw = Xraw(valid, :);
Cmat = Cmat(valid, :);
patient_ids = patient_ids(valid);
if isfield(opts,"SubsampleSeed") && opts.SubsampleSeed > 0
    % v2: stratified (by death) subsample of the analysis cohort for the
    % conditional stability analysis; screening is not repeated.
    stream = RandStream("twister","Seed",opts.SubsampleSeed);
    keep_rows = false(numel(TT),1);
    for stratum = [true false]
        idx = find(Delta == stratum);
        m = round(opts.SubsampleFraction*numel(idx));
        keep_rows(idx(randperm(stream,numel(idx),m))) = true;
    end
    TT = TT(keep_rows); Delta = Delta(keep_rows);
    Zraw = Zraw(keep_rows,:); Xraw = Xraw(keep_rows,:);
    Cmat = Cmat(keep_rows,:); patient_ids = patient_ids(keep_rows);
end
duplicate_map = [];
if isfield(opts,"MergeDuplicateAdditive") && logical(opts.MergeDuplicateAdditive)
    % Identical observed CNA profiles are not separately identifiable within
    % a margin; keep the first label of each class and record the classes.
    [~,keep,class_id] = unique(Xraw','rows','stable');
    duplicate_map = table(x_features(:),x_features(keep(class_id)), ...
        'VariableNames',{'feature','representative'});
    Xraw = Xraw(:,keep);
    x_features = x_features(keep);
end
if use_clinical
    Cmat(:,1) = (Cmat(:,1)-mean(Cmat(:,1)))/std(Cmat(:,1));
    Zraw = [Cmat, Zraw];
    z_features = [clinical_names; z_features(:)];
end

[n, p] = size(Zraw);
q = size(Xraw, 2);
linear_penalty_factor = ones(p,1);
linear_penalty_factor(1:numel(clinical_names)) = 0;
n0 = choose_bin_count(n,opts.BaselineBinCount);
block_size = 10;
q_size = 5;
a = opts.MICSharpness;
mpl_step_options = mpl_step_options_from_opts(opts);

[T, Z, X, B, status, T1, Z1, X1, B1, status1, num_basis, spline_knots] = ...
    real_data_design(TT, Delta, Zraw, Xraw, opts.NumInternalKnots, opts.SplineBasisScaling);

fprintf('Real-data PLAPH: n=%d, p=%d, q=%d, num_basis=%d, events=%d, censored=%d\n', ...
    n, p, q, num_basis, sum(status), sum(~status));

[ini_beta, ini_phi, ini_spline_beta, ini_spline_phi] = initial_value( ...
    n, Z, X, B, 1e-5, status, a, block_size, q_size, ...
    opts.MICPenaltyScale, opts.MICPenaltyNormalization, opts.InitialLinearStrategy, ...
    opts.InitialLinearCap,opts.InitialSplineWeighting,opts.InitialLinearPilotScale, ...
    opts.InitialSplineStrategy,opts.InitialSplinePilotScale, ...
    opts.MICPenaltyScaleFailure,opts.MICPenaltyScaleCensoring, ...
    opts.MICPenaltyMarginSpecific);
decomposable_pilot_options = decomposable_pilot_options_from_opts(opts);
preserve_pilot_initial_spline = decomposable_pilot_options.iterations > 0 && ...
    lower(string(decomposable_pilot_options.start)) == "initial";
if opts.SplineGroupPenaltyType == "group_iht" && ...
        opts.SplineGroupIHTWarmupIterations == 0 && ...
        ~preserve_pilot_initial_spline
    ini_spline_beta(:) = 0;
    ini_spline_phi(:) = 0;
end

[mpl_beta, mpl_phi, bspline_beta, bspline_phi, seta, gama, H_T, H_C, tau_hat, alph_hat, mpl_diag, support_state] = MPL( ...
    n, q, n0, ini_beta, ini_phi, ini_spline_beta, ini_spline_phi, ...
    Z1, X1, B1, T1, 1e-5, status1, a, char(scenario.fit_copula), scenario.tau, block_size, q_size, ...
    opts.ThresholdLinearFailure, opts.ThresholdLinearCensoring, opts.ThresholdSpline, opts.MaxOuterIter, ...
    opts.MICPenaltyScale, opts.MICPenaltyNormalization, opts.MPLStepScale, opts.CollectDiagnostics, ...
    opts.SplineThresholdStage, opts.SplineSelectionMetric, opts.SplineGroupPenaltyScale, ...
    opts.SplineGroupPenaltyType, opts.SplineGroupPenaltyNormalization, opts.SplineGroupMICSharpness, ...
    opts.SplineGroupMICContinuation, opts.SplineGroupMICSharpnessPath, ...
    estimate_tau, ...
    opts.LinearMICForm, opts.LinearKinkCoefficient, ...
    opts.LinearL1Type, opts.LinearL1Scale, opts.LinearRidgeScale, opts.AdaptiveWeightGamma, ...
    opts.AdaptiveWeightEpsilon, opts.AdaptiveWeightCap, ...
    opts.SplineGroupMICForm, opts.SplineGroupKinkCoefficient, ...
    opts.SplineGroupPenaltyScaleFailure, opts.SplineGroupPenaltyScaleCensoring, ...
    opts.SplineGroupMCPConcavity, opts.SplineGroupIHTThresholdFailure, ...
    opts.SplineGroupIHTThresholdCensoring, opts.SplineGroupIHTWarmupIterations, ...
    opts.PostSelectionRefitType, opts.SplineGroupIHTContinuationIterations, ...
    opts.SplineGroupIHTStartFraction, mpl_step_options, ...
    decomposable_pilot_options,[],[],[],[], ...
    opts.PostSelectionRefitMaxIter,opts.BaselineRoughnessScale, ...
    opts.MICPenaltyScaleFailure,opts.MICPenaltyScaleCensoring, ...
    opts.MICPenaltyMarginSpecific,linear_penalty_factor);

result = struct();
result.scenario = scenario;
result.df_rank = joint_selected_rank(Z1,B1,num_basis, ...
    support_state.final_selected_linear_T_mask,support_state.final_selected_additive_T_mask) + ...
    joint_selected_rank(Z1,B1,num_basis, ...
    support_state.final_selected_linear_C_mask,support_state.final_selected_additive_C_mask) + ...
    numel(seta) + numel(gama);
result.selection_gate_pass = mpl_diag.mpl_selection_gate_pass;
result.selection_status_code = mpl_diag.mpl_selection_status_code;
result.selection_gate_maximum = mpl_diag.mpl_selection_gate_maximum;
result.refit_stationarity_pass = mpl_diag.mpl_oracle_refit_stationarity_pass;
result.refit_stationarity_inf = mpl_diag.mpl_oracle_refit_final_stationarity_inf;
result.accepted = result.selection_gate_pass && result.refit_stationarity_pass;
result.unpenalized_clinical = clinical_names;
result.duplicate_additive_map = duplicate_map;
result.n = n;
result.p = p;
result.q = q;
result.num_basis = num_basis;
result.n_failure_events = sum(status);
result.n_censoring_events = sum(~status);
result.censor_rate = mean(~status);
result.tau_hat = tau_hat;
result.alpha_hat = alph_hat;
result.estimate_tau = estimate_tau;
result.initial_beta = ini_beta;
result.initial_phi = ini_phi;
result.initial_bspline_beta = ini_spline_beta;
result.initial_bspline_phi = ini_spline_phi;
result.preserved_initial_spline_for_pilot = preserve_pilot_initial_spline;
result.beta = mpl_beta;
result.phi = mpl_phi;
result.bspline_beta = bspline_beta;
result.bspline_phi = bspline_phi;
result.seta = seta;
result.gama = gama;
result.H_T = H_T;
result.H_C = H_C;
result.diagnostics = mpl_diag;
result.pilot_state = support_state;
result.selection_support = struct( ...
    'linear_failure',logical(support_state.final_selected_linear_T_mask(:)), ...
    'linear_censoring',logical(support_state.final_selected_linear_C_mask(:)), ...
    'additive_failure',logical(support_state.final_selected_additive_T_mask(:)), ...
    'additive_censoring',logical(support_state.final_selected_additive_C_mask(:)), ...
    'stage',"frozen_before_postselection_refit");
result.patient_ids = patient_ids;
result.z_features = z_features;
result.x_features = x_features;
result.spline_knots = spline_knots;
result.linear_failure = selected_linear_table(z_features,mpl_beta, ...
    opts.ThresholdLinearFailure,result.selection_support.linear_failure);
result.linear_censoring = selected_linear_table(z_features,mpl_phi, ...
    opts.ThresholdLinearCensoring,result.selection_support.linear_censoring);
result.additive_failure = selected_additive_table(x_features,B1, ...
    bspline_beta,num_basis,q,opts.ThresholdSpline, ...
    opts.SplineSelectionMetric,result.selection_support.additive_failure);
result.additive_censoring = selected_additive_table(x_features,B1, ...
    bspline_phi,num_basis,q,opts.ThresholdSpline, ...
    opts.SplineSelectionMetric,result.selection_support.additive_censoring);
result.function_curves_failure = real_data_function_curves(x_features, Xraw, spline_knots, bspline_beta, num_basis, opts.SplineBasisScaling);
result.function_curves_censoring = real_data_function_curves(x_features, Xraw, spline_knots, bspline_phi, num_basis, opts.SplineBasisScaling);
result.elapsed_seconds = toc;
end


function scenario = normalize_real_data_scenario(scenario)
if ~isfield(scenario, "scenario_id")
    scenario.scenario_id = "real_data";
end
if ~isfield(scenario, "fit_copula")
    scenario.fit_copula = "Frank";
end
if ~isfield(scenario, "tau")
    scenario.tau = 0.2;
end
if ~isfield(scenario, "p")
    scenario.p = 0;
end
if ~isfield(scenario, "q")
    scenario.q = 0;
end
scenario.scenario_id = string(scenario.scenario_id);
scenario.fit_copula = string(scenario.fit_copula);
end


function [clinical, patient_ids, Z, X, z_features, x_features] = read_real_data_working_set(scenario)
clinical = readtable(char(scenario.clinical_file), "VariableNamingRule", "preserve");
z_table = readtable(char(scenario.z_file), "VariableNamingRule", "preserve");
x_table = readtable(char(scenario.x_file), "VariableNamingRule", "preserve");
z_feature_table = readtable(char(scenario.z_feature_file), "VariableNamingRule", "preserve");
x_feature_table = readtable(char(scenario.x_feature_file), "VariableNamingRule", "preserve");

z_ids = strtrim(string(z_table{:, 1}));
x_ids = strtrim(string(x_table{:, 1}));
clinical_ids = strtrim(string(clinical.patient_id));
[common_zc, iz0, ic0] = intersect(z_ids, clinical_ids, 'stable');
[patient_ids, ia, ix] = intersect(common_zc, x_ids, 'stable');
iz = iz0(ia);
ic = ic0(ia);
if numel(patient_ids) < 100
    error('PLAPH:RealDataPatientMismatch', 'Too few patient IDs are shared across clinical, Z, and X files.');
end
clinical = clinical(ic, :);
z_table = z_table(iz, :);
x_table = x_table(ix, :);
Z = table2array(z_table(:, 2:end));
X = table2array(x_table(:, 2:end));
z_features = string(z_feature_table.feature);
x_features = string(x_feature_table.feature);
end


function [T, Z, X, B, status, T1, Z1, X1, B1, status1, num_basis, spline_knots] = ...
    real_data_design(TT, Delta, ZZ, XX, num_internal_knots, spline_basis_scaling)
spline_basis_scaling = lower(string(spline_basis_scaling));
[n, q] = size(XX);

order = plaph_spline_order();
total_knots = num_internal_knots + 2;
num_basis = num_internal_knots + order;
BB = zeros(n, num_basis*q);
spline_knots = cell(q, 1);

for i = 1:q
    knot_grid = quantile(XX(:, i), linspace(0, 1, total_knots));
    knot_grid = unique(knot_grid, 'stable');
    if numel(knot_grid) < total_knots
        knot_grid = linspace(min(XX(:, i)), max(XX(:, i)), total_knots);
    end
    internal_knots = knot_grid(2:(end-1));
    left_boundary = knot_grid(1);
    right_boundary = knot_grid(end);
    if right_boundary <= left_boundary
        right_boundary = left_boundary + 1e-6;
    end
    knots = [left_boundary*ones(1, order), internal_knots, right_boundary*ones(1, order)];
    spline_knots{i} = knots;
    B_eval = zeros(n, num_basis);
    for j = 1:num_basis
        [B_eval(:, j), ~] = bspline_basis(j-1, order, knots, XX(:, i));
    end
    BB(:, ((i-1)*num_basis+1):(i*num_basis)) = B_eval;
end

if spline_basis_scaling == "zscore"
    BB = BB - mean(BB, 1);
    BB = BB ./ (std(BB, 0, 1) + 1e-8);
elseif spline_basis_scaling == "center"
    BB = BB - mean(BB, 1);
elseif spline_basis_scaling == "whiten"
    BB = whiten_spline_basis_groups(BB, num_basis, q);
elseif spline_basis_scaling == "none"
    BB = BB;
else
    error('PLAPH:UnknownSplineBasisScaling', 'Unknown SplineBasisScaling: %s', spline_basis_scaling);
end
BB(~isfinite(BB)) = 0;

[T, I] = sort(TT, 'descend');
Z = ZZ(I, :);
X = XX(I, :);
status = Delta(I);
B = BB(I, :);

[T1, I1] = sort(TT, 'ascend');
Z1 = ZZ(I1, :);
X1 = XX(I1, :);
status1 = Delta(I1);
B1 = BB(I1, :);
end


function tab = selected_linear_table(feature_names, coef, threshold, selected)
if nargin < 4 || isempty(selected)
    selected = abs(coef) > threshold;
end
selected = logical(selected(:));
if numel(selected) ~= numel(coef)
    error('PLAPH:RealDataLinearSupportDimension', ...
        'The frozen linear support has the wrong dimension.');
end
tab = table(feature_names(:), coef(:), abs(coef(:)), selected(:), ...
    'VariableNames', {'feature', 'coefficient', 'abs_coefficient', 'selected'});
tab = sortrows(tab, "abs_coefficient", "descend");
end


function tab = selected_additive_table(feature_names, B, coef, num_basis, q, threshold, metric, selected)
metric = lower(string(metric));
coef_norm = zeros(q, 1);
effect_norm = zeros(q, 1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    coef_norm(j) = norm(coef(idx), 2);
    effect_norm(j) = sqrt(mean((B(:, idx)*coef(idx)).^2));
end
if metric == "effect"
    score = effect_norm;
else
    score = coef_norm;
end
if nargin < 8 || isempty(selected)
    selected = score > threshold;
end
selected = logical(selected(:));
if numel(selected) ~= q
    error('PLAPH:RealDataAdditiveSupportDimension', ...
        'The frozen additive support has the wrong dimension.');
end
tab = table(feature_names(:), coef_norm, effect_norm, score, selected, ...
    'VariableNames', {'feature', 'coef_group_norm', 'effect_rms', 'selection_score', 'selected'});
tab = sortrows(tab, "selection_score", "descend");
end


function curves = real_data_function_curves(feature_names, X, spline_knots, coef, num_basis, spline_basis_scaling)
q = numel(feature_names);
grid_n = 101;
rows = cell(q, 1);
order = plaph_spline_order();
for j = 1:q
    xj = X(:, j);
    grid = linspace(quantile(xj, 0.02), quantile(xj, 0.98), grid_n)';
    B_grid = zeros(grid_n, num_basis);
    for k = 1:num_basis
        [B_grid(:, k), ~] = bspline_basis(k-1, order, spline_knots{j}, grid);
    end
    train_B = zeros(size(X, 1), num_basis);
    for k = 1:num_basis
        [train_B(:, k), ~] = bspline_basis(k-1, order, spline_knots{j}, xj);
    end
    if lower(string(spline_basis_scaling)) == "zscore"
        mu = mean(train_B, 1);
        sd = std(train_B, 0, 1) + 1e-8;
        B_grid = (B_grid - mu) ./ sd;
    elseif lower(string(spline_basis_scaling)) == "center"
        B_grid = B_grid - mean(train_B, 1);
    elseif lower(string(spline_basis_scaling)) == "whiten"
        [~, whitening_transform, basis_center] = whiten_spline_block(train_B);
        B_grid = (B_grid - basis_center) * whitening_transform;
    end
    idx = ((j-1)*num_basis+1):(j*num_basis);
    fit = B_grid * coef(idx);
    fit = fit - mean(fit);
    rows{j} = table(repmat(feature_names(j), grid_n, 1), grid, fit, ...
        'VariableNames', {'feature', 'x', 'centered_fit'});
end
curves = vertcat(rows{:});
end


function BB = whiten_spline_basis_groups(BB, num_basis, q)
% Standardize each centered block so its empirical fitted-effect norm equals
% the Euclidean coefficient norm on the estimable spline subspace.
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    BB(:, idx) = whiten_spline_block(BB(:, idx));
end
end


function [B_white, transform, center] = whiten_spline_block(B_raw)
center = mean(B_raw, 1);
B_centered = B_raw - center;
gram = (B_centered' * B_centered) / size(B_centered, 1);
gram = (gram + gram') / 2;
[vectors, values] = eig(gram);
eigenvalues = real(diag(values));
[eigenvalues, order] = sort(eigenvalues, 'descend');
vectors = real(vectors(:, order));

% Fix eigenvector signs so repeated construction gives the same basis.
for k = 1:size(vectors, 2)
    [~, pivot] = max(abs(vectors(:, k)));
    if vectors(pivot, k) < 0
        vectors(:, k) = -vectors(:, k);
    end
end

largest = max([eigenvalues; 0]);
tolerance = max(1e-12, largest * size(gram, 1) * eps);
inverse_sqrt = zeros(size(eigenvalues));
estimable = eigenvalues > tolerance;
inverse_sqrt(estimable) = 1 ./ sqrt(eigenvalues(estimable));
transform = vectors * diag(inverse_sqrt);
B_white = B_centered * transform;
B_white(~isfinite(B_white)) = 0;
end


%% ===================================================
%                 initial_value()
% ============================================================
function [initial_beta,initial_phi,initial_spline_beta,initial_spline_phi]= initial_value(n,Z,X,B, r,status,a,block_size,q_size,mic_penalty_scale,mic_penalty_normalization,initial_linear_strategy,initial_linear_cap,initial_spline_weighting,initial_linear_pilot_scale,initial_spline_strategy,initial_spline_pilot_scale,mic_penalty_scale_T,mic_penalty_scale_C,mic_penalty_margin_specific)
if nargin < 10 || isempty(mic_penalty_scale)
    mic_penalty_scale = 1;
end
if nargin < 11 || isempty(mic_penalty_normalization)
    mic_penalty_normalization = "sum";
end
if nargin < 12 || isempty(initial_linear_strategy)
    initial_linear_strategy = "legacy";
end
if nargin < 13 || isempty(initial_linear_cap)
    initial_linear_cap = Inf;
end
if nargin < 14 || isempty(initial_spline_weighting)
    initial_spline_weighting = "legacy_tail";
end
if nargin < 15 || isempty(initial_linear_pilot_scale)
    initial_linear_pilot_scale = 1;
end
if nargin < 16 || isempty(initial_spline_strategy)
    initial_spline_strategy = "legacy";
end
if nargin < 17 || isempty(initial_spline_pilot_scale)
    initial_spline_pilot_scale = 1;
end
if nargin < 18 || isempty(mic_penalty_scale_T)
    mic_penalty_scale_T = mic_penalty_scale;
end
if nargin < 19 || isempty(mic_penalty_scale_C)
    mic_penalty_scale_C = mic_penalty_scale;
end
if nargin < 20 || isempty(mic_penalty_margin_specific)
    mic_penalty_margin_specific = false;
end
initial_linear_strategy = lower(string(initial_linear_strategy));
initial_spline_weighting = lower(string(initial_spline_weighting));
initial_spline_strategy = lower(string(initial_spline_strategy));
initial_linear_warmup = initial_linear_strategy == "legacy";
initial_spline_warmup = initial_spline_strategy == "legacy";
[~,p] = size(Z);
[~,q] = size(X);
[~,qs] = size(B);
beta = zeros(p,1);
phi = zeros(p,1);
beta1 = zeros(p,1);
phi1 = zeros(p,1);
bspline_beta = zeros(qs,1);
bspline_phi = zeros(qs,1);

n0 = sum(status);
lambda0_T = mic_penalty_lambda(n0, n, mic_penalty_scale_T, ...
    mic_penalty_normalization);
if mic_penalty_margin_specific
    lambda0_C = mic_penalty_lambda(n-n0, n, mic_penalty_scale_C, ...
        mic_penalty_normalization);
else
    lambda0_C = lambda0_T;
end

%% Set up block coordinate descent parameters
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% Construct block indices
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


%% Set up block coordinate descent parameters--spline
num_basis = qs/q;
block_size_spline = block_size*num_basis;  % num_basis  % number of spline basis functions
num_blocks_spline = ceil(q / block_size);
blocks_spline = cell(num_blocks_spline, 1);

%% Construct block indices
for b = 1:num_blocks_spline
    start_idx1 = (b-1)*block_size_spline + 1;
    end_idx1 = min(b*block_size_spline, q*num_basis);
    blocks_spline{b} = start_idx1:end_idx1;
end


switch initial_spline_weighting
    case "none"
        weights = ones(n,1);
    case "legacy_tail"
        % Historical simulation-only warm start. These weights use the known
        % identities of X1 and X2 and are not part of the general estimator.
        weights1 = 1 + 2.0*(X(:,1) > 1.8) + 0.35*(X(:,1) < -1.8);
        weights2 = 1 + 3.0*(X(:,2) > 2.0) + 3.0*(X(:,2) < -2.0);
        weights = weights1+weights2;
    otherwise
        error('PLAPH:UnknownInitialSplineWeighting', ...
            'Unsupported InitialSplineWeighting: %s', initial_spline_weighting);
end
 

k = 1; err1 = ~initial_spline_warmup; err2 = ~initial_spline_warmup;
err3 = ~initial_linear_warmup; err4 = ~initial_linear_warmup;
tk1 = 4; tk2 = 4;

while k<=1000 && (err1==0 || err2==0 || err3==0 || err4==0)
    %k

    if initial_spline_warmup
    L_prime2 = -(sum((status.*weights).*(B-cumsum((safe_exp(Z*beta+B*bspline_beta).*B))./cumsum(safe_exp(Z*beta+B*bspline_beta)))))'/n;
    % L_prime2 = -(sum(status.*(B-cumsum((safe_exp(Z*beta+B*bspline_beta).*B))./cumsum(safe_exp(Z*beta+B*bspline_beta)))))'/n;
    %bspline_beta1 =  bspline_beta- L_prime2/tk1;

    if q<=q_size
        bspline_beta1 =  bspline_beta- L_prime2/tk1;
    else
        block_order_spline = randperm(num_blocks_spline);
        for b = block_order_spline
            current_block_spline = blocks_spline{b};
            bspline_beta_block = bspline_beta(current_block_spline);
            bspline_beta_new1 = bspline_beta_block - ...
                L_prime2(current_block_spline)/tk1;
            bspline_beta1(current_block_spline,1) = bspline_beta_new1;
        end
        delta = norm(bspline_beta - bspline_beta1,2)^2 / ...
            (norm(bspline_beta,2)^2 + 1e-8);
        if delta >= 1e-6
            for iter = 2:50
                randperm(num_blocks_spline);
            end
        end
    end

    w1 = bspline_beta1-bspline_beta;
    err1 = norm(w1,2)^2 <= r*norm(bspline_beta,2)^2;
    bspline_beta = bspline_beta1;

    risk_phi_init = safe_exp(Z*phi+B*bspline_phi);
    L_prime4 = -(sum(((1-status).*weights).*(B-cumsum((risk_phi_init.*B))./cumsum(risk_phi_init))))'/n;
    % bspline_phi1 =  bspline_phi- L_prime4/tk2;

    if q<=q_size
        bspline_phi1 =  bspline_phi- L_prime4/tk2;
    else
        block_order_spline = randperm(num_blocks_spline);
        for b = block_order_spline
            current_block_spline = blocks_spline{b};
            bspline_phi_block = bspline_phi(current_block_spline);
            bspline_phi_new1 = bspline_phi_block - ...
                L_prime4(current_block_spline)/tk2;
            bspline_phi1(current_block_spline,1) = bspline_phi_new1;
        end
        delta = norm(bspline_phi - bspline_phi1,2)^2 / ...
            (norm(bspline_phi,2)^2 + 1e-8);
        if delta >= 1e-6
            for iter = 2:50
                randperm(num_blocks_spline);
            end
        end
    end

    w2 = bspline_phi1-bspline_phi;
    err2 = norm(w2,2)^2 <= r*norm(bspline_phi,2)^2;
    bspline_phi = bspline_phi1;
    else
        err1 = 1;
        err2 = 1;
    end

    if initial_linear_warmup
    L_prime1 = -(sum(status.*(Z-cumsum((safe_exp(Z*beta+B*bspline_beta).*Z))./...
        cumsum(safe_exp(Z*beta+B*bspline_beta)))))'/n;

    if p<=30
        beta1 =  beta - L_prime1/tk1;
    else
        block_order = randperm(num_blocks);
        for b = block_order
            current_block = blocks{b};
            beta_block = beta(current_block);
            beta_new1 = beta_block - L_prime1(current_block)/(tk1*10);
            beta1(current_block,1) = beta_new1;
        end
        delta = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
        if delta >= 1e-6
            for iter = 2:50
                randperm(num_blocks);
            end
        end
    end
    
    w3 = beta1-beta;
    err3 = norm(w3,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;

    risk_phi = safe_exp(Z*phi+B*bspline_phi);
    denom_phi = max(cumsum(risk_phi), 1e-8);
    L_prime3 = -(sum((1-status).*(Z-cumsum(risk_phi.*Z)./denom_phi)))'/n;
    L_prime3(~isfinite(L_prime3)) = 0;

    if p <= 30
        phi1 =  phi - L_prime3/tk2;
    else
        block_order2 = randperm(num_blocks);
        for b = block_order2
            current_block = blocks{b};
            phi_block = phi(current_block);
            phi_new1 = phi_block - L_prime3(current_block)/(tk2*10);
            phi1(current_block,1) = phi_new1;
        end
        delta2 = norm(phi - phi1,2)^2 / (norm(phi,2)^2 + 1e-8);
        if delta2 >= 1e-6
            for iter = 2:50
                randperm(num_blocks);
            end
        end
    end

    w4 = phi1-phi;
    err4 = norm(w4,2)^2 <= r*norm(phi,2)^2;
    phi = phi1;
    else
        err3 = 1;
        err4 = 1;
    end

    k = k+1;

end


%% ============================================================
if initial_linear_strategy == "lasso"
    [beta, phi] = initial_linear_lasso_pilot( ...
        Z, B*bspline_beta, B*bspline_phi, status, initial_linear_pilot_scale);
end

if initial_spline_strategy == "group_lasso"
    [bspline_beta,bspline_phi] = initial_spline_group_lasso_pilot( ...
        B,Z*beta,Z*phi,status,q,initial_spline_pilot_scale);
    if initial_linear_strategy == "lasso"
        % One additional block sweep reduces confounding between the scalar
        % and additive pilot components while retaining convex subproblems.
        [beta,phi] = initial_linear_lasso_pilot( ...
            Z,B*bspline_beta,B*bspline_phi,status,initial_linear_pilot_scale);
        [bspline_beta,bspline_phi] = initial_spline_group_lasso_pilot( ...
            B,Z*beta,Z*phi,status,q,initial_spline_pilot_scale);
    end
end

if isfinite(initial_linear_cap)
    beta = sign(beta).*min(abs(beta), initial_linear_cap);
    phi = sign(phi).*min(abs(phi), initial_linear_cap);
end

if initial_linear_strategy == "continuation"
    lambda_mean_T = mic_penalty_lambda(n0, n, mic_penalty_scale_T, "mean");
    lambda_sum_T = mic_penalty_lambda(n0, n, mic_penalty_scale_T, "sum");
    lambda_path_T = unique([lambda_mean_T, ...
        sqrt(lambda_mean_T*lambda_sum_T),lambda_sum_T],"stable");
    if mic_penalty_margin_specific
        lambda_mean_C = mic_penalty_lambda(n-n0,n,mic_penalty_scale_C,"mean");
        lambda_sum_C = mic_penalty_lambda(n-n0,n,mic_penalty_scale_C,"sum");
        lambda_path_C = unique([lambda_mean_C, ...
            sqrt(lambda_mean_C*lambda_sum_C),lambda_sum_C],"stable");
    else
        lambda_path_C = lambda_path_T;
    end
    for lambda_iter = 1:numel(lambda_path_T)
        [beta, phi] = initial_linear_mic_refine( ...
            beta, phi, bspline_beta, bspline_phi, Z, B, status, a, ...
            lambda_path_T(lambda_iter),lambda_path_C(lambda_iter), ...
            blocks,p,block_size,r,4,4,200);
    end
end

if initial_linear_strategy == "zero" || initial_linear_strategy == "lasso"
    % Preserve the requested sparse linear start. The smooth-MIC refinement
    % below is part of the legacy initialization and would otherwise replace
    % zero or L1-pilot coefficients before kinked-MIC iterations begin.
    initial_beta = beta;
    initial_phi = phi;
    initial_spline_beta = bspline_beta;
    initial_spline_phi = bspline_phi;
    return;
end

%% ============================================================
k = 1; err5 = 0; tk3 = 4;

while k<=1000 && err5==0
    %k

    risk_beta = safe_exp(Z*beta+B*bspline_beta);
    denom_beta = max(cumsum(risk_beta), 1e-8);
    L_prime3 = -(sum(status.*(Z-cumsum(risk_beta.*Z)./denom_beta)))'/n;
    L_prime3(~isfinite(L_prime3)) = 0;
    W2 = diag(smooth_mic_curvature_weight(beta, a, lambda0_T, 1e-8));
    u1 = eye(p) + 2*W2/tk3 + 1e-6*eye(p);
    beta_tilde = beta - L_prime3/tk3;
    beta_tilde(~isfinite(beta_tilde)) = 0;
    beta1 = u1\beta_tilde;
    beta1(~isfinite(beta1)) = 0;

    w5 = beta1-beta;
    err5 = norm(w5,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

end

beta2 = beta.*(abs(beta)>=2.5*1e-4);
opt_beta = beta2;


%% ============================================================
k = 1; err6 = 0; tk4 = 4;

while k<=1000 && err6==0
    % k
    risk_phi = safe_exp(Z*phi+B*bspline_phi);
    denom_phi = max(cumsum(risk_phi), 1e-8);
    L_prime6 = -(sum((1-status).*(Z-cumsum(risk_phi.*Z)./denom_phi)))'/n;
    L_prime6(~isfinite(L_prime6)) = 0;

    W5 = diag(smooth_mic_curvature_weight(phi, a, lambda0_C, 1e-8));
    u2 = eye(p) + 2*W5/tk4 + 1e-6*eye(p);
    phi_tilde = phi - L_prime6/tk4;
    phi_tilde(~isfinite(phi_tilde)) = 0;
    phi1 = u2\phi_tilde;
    phi1(~isfinite(phi1)) = 0;


    w6 = phi1-phi;
    err6 = norm(w6,2)^2 <= r*norm(phi,2)^2;
    phi = phi1;
    k = k+1;

end

phi2 = phi.*(abs(phi)>=2.5*1e-4);
opt_phi = phi2;

% %% ============================================================
% k = 1; err7 = 0; tk5 = 4;
% while k<=1000 && err7==0
% L_prime2 = -(sum(status.*(B-cumsum((safe_exp(Z*beta+B*bspline_beta).*B))./cumsum(safe_exp(Z*beta+B*bspline_beta)))))'/n;
%     bspline_beta1 =  bspline_beta- L_prime2/tk5;
% 
%     w7 = bspline_beta1-bspline_beta;
%     err7 = norm(w7,2)^2 <= r*norm(bspline_beta,2)^2;
%     bspline_beta = bspline_beta1;
% 
%     k = k+1;
% 
% end
% 
% %% ============================================================
% k = 1; err8 = 0; tk6 = 4;
% while k<=1000 && err8==0
% 
%     L_prime4 = -(sum((1-status).*(B-cumsum((exp(Z*beta+B*bspline_phi).*B))./cumsum(exp(Z*beta+B*bspline_phi)))))'/n;
%     bspline_phi1 =  bspline_phi- L_prime4/tk6;
% 
%     w8 = bspline_phi1-bspline_phi;
%     err8 = norm(w8,2)^2 <= r*norm(bspline_phi,2)^2;
%     bspline_phi = bspline_phi1;
% 
%     k = k+1;
% 
% end


initial_beta = opt_beta;
initial_phi = opt_phi;
initial_spline_beta = bspline_beta;
initial_spline_phi = bspline_phi;

end


function [beta, phi] = initial_linear_lasso_pilot(Z, offset_beta, offset_phi, status, pilot_scale)
% Margin-specific decomposable diagnostic on the screened working set.
n = size(Z,1);
p = size(Z,2);
lambda_beta = pilot_scale*sqrt(max(sum(status),1)*log(max(p,2)))/n;
lambda_phi = pilot_scale*sqrt(max(sum(1-status),1)*log(max(p,2)))/n;
beta = cox_lasso_pilot(Z, offset_beta, status, lambda_beta);
phi = cox_lasso_pilot(Z, offset_phi, 1-status, lambda_phi);
end


function [bspline_beta,bspline_phi] = initial_spline_group_lasso_pilot( ...
    B,offset_beta,offset_phi,status,q,pilot_scale)
% Convex margin pilots with group norms defined by fitted-effect RMS.
n = size(B,1);
num_basis = size(B,2)/q;
rank_per_group = max(num_basis-1,1);
lambda_beta = pilot_scale*sqrt(max(sum(status),1)* ...
    (rank_per_group+log(max(q,2))))/n;
lambda_phi = pilot_scale*sqrt(max(sum(1-status),1)* ...
    (rank_per_group+log(max(q,2))))/n;
bspline_beta = cox_group_lasso_pilot( ...
    B,offset_beta,status,lambda_beta,num_basis,q);
bspline_phi = cox_group_lasso_pilot( ...
    B,offset_phi,1-status,lambda_phi,num_basis,q);
end


function coef = cox_group_lasso_pilot(B,offset,events,lambda,num_basis,q)
metric = prepare_spline_group_metric(B,num_basis,q);
coef = zeros(size(B,2),1);
curvature = 4;
for iter = 1:500
    [loss,gradient] = cox_partial_loss_gradient(B,offset,events,coef);
    accepted = false;
    candidate = coef;
    for backtrack = 1:20
        candidate = apply_spline_group_lasso_metric( ...
            metric,coef,gradient,num_basis,q,curvature,lambda);
        step = candidate-coef;
        upper_bound = loss+gradient'*step+ ...
            0.5*curvature*spline_metric_squared_norm(metric,step,num_basis,q);
        candidate_loss = cox_partial_loss_gradient(B,offset,events,candidate);
        if isfinite(candidate_loss) && candidate_loss <= upper_bound+1e-10
            accepted = true;
            break;
        end
        curvature = 2*curvature;
    end
    if ~accepted
        break;
    end
    step_norm = sqrt(spline_metric_squared_norm(metric,candidate-coef,num_basis,q));
    coef_norm = sqrt(spline_metric_squared_norm(metric,coef,num_basis,q));
    coef = candidate;
    if step_norm <= 1e-6*(1+coef_norm)
        break;
    end
    curvature = max(curvature/1.25,1e-6);
end
coef(~isfinite(coef)) = 0;
end


function value = spline_metric_squared_norm(metric,coef,num_basis,q)
value = 0;
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    value = value+real(coef(idx)'*metric.gram{j}*coef(idx));
end
value = max(value,0);
end


function coef = cox_lasso_pilot(Z, offset, events, lambda)
coef = zeros(size(Z,2),1);
curvature = 4;
for iter = 1:500
    [loss, gradient] = cox_partial_loss_gradient(Z, offset, events, coef);
    accepted = false;
    candidate = coef;
    for backtrack = 1:20
        trial = coef - gradient/curvature;
        candidate = sign(trial).*max(abs(trial)-lambda/curvature, 0);
        candidate_loss = cox_partial_loss_gradient(Z, offset, events, candidate);
        step = candidate-coef;
        upper_bound = loss + gradient'*step + 0.5*curvature*(step'*step);
        if isfinite(candidate_loss) && candidate_loss <= upper_bound + 1e-10
            accepted = true;
            break;
        end
        curvature = 2*curvature;
    end
    if ~accepted
        break;
    end
    if norm(candidate-coef,2) <= 1e-6*(1+norm(coef,2))
        coef = candidate;
        break;
    end
    coef = candidate;
end
coef(~isfinite(coef)) = 0;
end


function [loss, gradient] = cox_partial_loss_gradient(Z, offset, events, coef)
n = size(Z,1);
eta = offset + Z*coef;
risk = safe_exp(eta);
denominator = max(cumsum(risk), 1e-12);
loss = -sum(events.*(eta-log(denominator)))/n;
if nargout > 1
    risk_mean = cumsum(risk.*Z)./denominator;
    gradient = -(sum(events.*(Z-risk_mean)))'/n;
    gradient(~isfinite(gradient)) = 0;
end
end


function [beta, phi] = initial_linear_mic_refine(beta, phi, bspline_beta, bspline_phi, Z, B, status, a, lambda0_T, lambda0_C, blocks, p, block_size, r, tk_beta, tk_phi, max_iter)
%% One MIC/SELO refinement stage for the initial linear coefficients.
beta1 = beta;
phi1 = phi;
err_beta = 0;
err_phi = 0;
k = 1;

while k <= max_iter && (err_beta == 0 || err_phi == 0)
    risk_beta = safe_exp(Z*beta+B*bspline_beta);
    denom_beta = max(cumsum(risk_beta), 1e-8);
    L_beta = -(sum(status.*(Z-cumsum(risk_beta.*Z)./denom_beta)))'/size(Z,1);
    L_beta(~isfinite(L_beta)) = 0;

    beta1 = mic_linear_cd_step(beta,L_beta,a,lambda0_T, ...
        blocks,p,block_size,tk_beta);
    beta1(~isfinite(beta1)) = 0;
    w_beta = beta1 - beta;
    err_beta = norm(w_beta,2)^2 <= r*(norm(beta,2)^2 + 1e-8);
    beta = beta1;

    risk_phi = safe_exp(Z*phi+B*bspline_phi);
    denom_phi = max(cumsum(risk_phi), 1e-8);
    L_phi = -(sum((1-status).*(Z-cumsum(risk_phi.*Z)./denom_phi)))'/size(Z,1);
    L_phi(~isfinite(L_phi)) = 0;

    phi1 = mic_linear_cd_step(phi,L_phi,a,lambda0_C, ...
        blocks,p,block_size,tk_phi);
    phi1(~isfinite(phi1)) = 0;
    w_phi = phi1 - phi;
    err_phi = norm(w_phi,2)^2 <= r*(norm(phi,2)^2 + 1e-8);
    phi = phi1;
    k = k + 1;
end
end


function coef_new = mic_linear_cd_step(coef, grad, a, lambda0, blocks, p, block_size, tk)
coef_new = coef;
if p <= 30
    W = diag(smooth_mic_curvature_weight(coef, a, lambda0, 1e-8));
    u = eye(p) + 2*W/tk + 1e-6*eye(p);
    coef_tilde = coef - grad/tk;
    coef_new = u\coef_tilde;
    return;
end

block_order = randperm(numel(blocks));
for b = block_order
    current_block = blocks{b};
    coef_block = coef(current_block);
    W = diag(smooth_mic_curvature_weight(coef_block, a, lambda0, 1e-8));
    u = eye(length(current_block)) + 2*W/tk;
    coef_tilde = coef_block - grad(current_block)/tk;
    coef_new(current_block,1) = u\coef_tilde;
end
delta = norm(coef - coef_new,2)^2/(norm(coef,2)^2 + 1e-8);
if delta >= 1e-6
    for iter = 2:50
        randperm(numel(blocks));
    end
end
end


%% ===================================================
%                 discrBinNA()
% ============================================================

% % Discretizing a sample of n survival times into sub-intervals
% % with no. of 'binCount' subjects in each sub-interval

%-- binwv: the bin widths
%-- ID: the bin ID of each subject
%-- binedg
%% Discretizing a sample of n survival times into sub-intervals with no. of 'binCount' subjects in each sub-interval
function [ID,binedg,binwv] = discrBinNA(T,binCount)
%% Discretize ordered survival times into flexible bins.
% This version does NOT require binCount to divide the sample size.
% The first m-1 bins have approximately binCount observations, and
% the last bin contains the remaining observations.
%
% Input:
%   T        : ordered survival times
%   binCount : target number of observations per bin
%
% Output:
%   ID     : bin ID for each observation
%   binedg : bin edges
%   binwv  : bin widths

len = length(T);
m = ceil(len/binCount);

ID = zeros(len,1);
binedg = zeros(m+1,1);
binedg(1) = min(0,T(1));

for j = 1:m
    start_idx = (j-1)*binCount + 1;
    end_idx = min(j*binCount, len);

    ID(start_idx:end_idx) = j;
    if j > 1
        % Put an interior boundary between adjacent ordered observations.
        % This gives every observation in bin j nonnegative partial exposure
        % in its current interval and avoids assigning an event to both bins.
        binedg(j) = 0.5*(T(start_idx-1)+T(start_idx));
    end
end

% Add a slightly enlarged final endpoint.
if len >= 2
    last_width = abs(T(end) - T(max(1,end-1)));
    binedg(m+1) = T(end) + 1e-2*last_width + 1e-8;
else
    binedg(m+1) = T(end) + 1e-8;
end

% Ensure positive bin widths even when ties occur.
binwv = diff(binedg);
binwv = max(binwv, 1e-8);

end



%% First order difference%---%For piecewise constant
function mat = mat1(psix,T,eps)
m = size(psix,2);
[~,index] = sort(T);
spsix= psix(index,:);
R = (diff(spsix))'*diff(spsix);
R = R+eps*eye(m); % v2: ridge on the diagonal only
mat = R;
end



%% ===========================================================
%                 bspline_basis()
% ============================================================
function [y,x] = bspline_basis(j,n,t,x)
% B-spline basis function value B(j,n) at x.

% Input arguments:
% j: interval index, 0 =< j < numel(t)-n
% n: B-spline order (2 for linear, 3 for quadratic, etc.)
% t: knot vector
% x (optional): value where the basis function is to be evaluated

% Output arguments:
% y: B-spline basis function value, nonzero for a knot span of n
% Copyright 2010 Levente Hunyadi

validateattributes(j, {'numeric'}, {'nonnegative','integer','scalar'});
validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
validateattributes(t, {'numeric'}, {'real','vector'});
assert(all( t(2:end)-t(1:end-1) >= 0 ),'Knot vector values should be nondecreasing.');
if nargin < 4
    x = linspace(t(n), t(end-n+1), 100);  % allocate points uniformly
else
    validateattributes(x, {'numeric'}, {'real','vector'});
end
assert(0 <= j && j < numel(t)-n, ...
    'Invalid interval index j = %d, expected 0 =< j < %d (0 =< j < numel(t)-n).', j, numel(t)-n);
y = bspline_basis_recurrence(j,n,t,x);

%%  bspline_basis_recurrence()
    function y = bspline_basis_recurrence(j,n,t,x)
        y = zeros(size(x));
        if n > 1
            b = bspline_basis(j,n-1,t,x);
            dn = x - t(j+1);
            dd = t(j+n) - t(j+1);
            if dd ~= 0  % indeterminate forms 0/0 are deemed to be zero
                y = y + b.*(dn./dd);
            end
            b = bspline_basis(j+1,n-1,t,x);
            dn = t(j+n+1) - x;
            dd = t(j+n+1) - t(j+1+1);
            if dd ~= 0
                y = y + b.*(dn./dd);
            end
        elseif t(j+2) < t(end)  % treat last element of knot vector as a special case
        	y(t(j+1) <= x & x < t(j+2)) = 1;
        else
        	y(t(j+1) <= x) = 1;
        end

    end

end



%% ===================================================
%                 MPL()
% ============================================================
function [opt_beta,opt_phi,opt_bspline_beta,opt_bspline_phi,opt_seta,opt_gama,H_T,H_C,tau_hat,alph_hat,mpl_diag,pilot_state,mpl_inference] =...
    MPL(n,q,n0,ini_beta,ini_phi,ini_bspline_beta,ini_bspline_phi,Z,X,B,T,r,Delta,a,copula_family,tau,block_size,q_size,threshold_linear_beta,threshold_linear_phi,threshold_spline,max_outer_iter,mic_penalty_scale,mic_penalty_normalization,mpl_step_scale,collect_diagnostics,spline_threshold_stage,spline_selection_metric,spline_group_penalty_scale,spline_group_penalty_type,spline_group_penalty_normalization,spline_group_mic_sharpness,spline_group_mic_continuation,spline_group_mic_sharpness_path,estimate_tau,linear_mic_form,linear_kink_coefficient,linear_l1_type,linear_l1_scale,linear_ridge_scale,adaptive_weight_gamma,adaptive_weight_epsilon,adaptive_weight_cap,spline_group_mic_form,spline_group_kink_coefficient,spline_group_penalty_scale_T,spline_group_penalty_scale_C,spline_group_mcp_concavity,spline_group_iht_threshold_T,spline_group_iht_threshold_C,spline_group_iht_warmup_iterations,post_selection_refit_type,spline_group_iht_continuation_iterations,spline_group_iht_start_fraction,mpl_step_options,decomposable_pilot_options,linear_mask_T,linear_mask_C,group_mask_T,group_mask_C,post_selection_refit_max_iter,baseline_roughness_scale,mic_penalty_scale_T,mic_penalty_scale_C,mic_penalty_margin_specific,linear_penalty_factor)
if nargin < 20 || isempty(threshold_linear_beta)
    threshold_linear_beta = 5e-4;
end
if nargin < 21 || isempty(threshold_linear_phi)
    threshold_linear_phi = threshold_linear_beta;
end
if nargin < 23 || isempty(mic_penalty_scale)
    mic_penalty_scale = 1;
end
if nargin < 24 || isempty(mic_penalty_normalization)
    mic_penalty_normalization = "sum";
end
if nargin < 25 || isempty(mpl_step_scale)
    mpl_step_scale = 1;
end
if nargin < 26 || isempty(collect_diagnostics)
    collect_diagnostics = false;
end
if nargin < 27 || isempty(spline_threshold_stage)
    spline_threshold_stage = "iterative";
end
spline_threshold_stage = lower(string(spline_threshold_stage));
if nargin < 28 || isempty(spline_selection_metric)
    spline_selection_metric = "coef";
end
spline_selection_metric = lower(string(spline_selection_metric));
if nargin < 29 || isempty(spline_group_penalty_scale)
    spline_group_penalty_scale = 15;
end
if nargin < 30 || isempty(spline_group_penalty_type)
    spline_group_penalty_type = "group_mic";
end
spline_group_penalty_type = lower(string(spline_group_penalty_type));
if nargin < 31 || isempty(spline_group_penalty_normalization)
    spline_group_penalty_normalization = "sum";
end
spline_group_penalty_normalization = lower(string(spline_group_penalty_normalization));
if nargin < 32 || isempty(spline_group_mic_sharpness)
    spline_group_mic_sharpness = 90;
end
if nargin < 33 || isempty(spline_group_mic_continuation)
    spline_group_mic_continuation = true;
end
if nargin < 34
    spline_group_mic_sharpness_path = [];
end
if nargin < 35 || isempty(estimate_tau)
    estimate_tau = true;
end
if nargin < 36 || isempty(linear_mic_form)
    linear_mic_form = "smooth";
end
if nargin < 37 || isempty(linear_kink_coefficient)
    linear_kink_coefficient = 0;
end
if nargin < 38 || isempty(linear_l1_type)
    linear_l1_type = "none";
end
if nargin < 39 || isempty(linear_l1_scale)
    linear_l1_scale = 0;
end
if nargin < 40 || isempty(linear_ridge_scale)
    linear_ridge_scale = 0;
end
if nargin < 41 || isempty(adaptive_weight_gamma)
    adaptive_weight_gamma = 1;
end
if nargin < 42 || isempty(adaptive_weight_epsilon)
    adaptive_weight_epsilon = 0.05;
end
if nargin < 43 || isempty(adaptive_weight_cap)
    adaptive_weight_cap = 100;
end
if nargin < 44 || isempty(spline_group_mic_form)
    spline_group_mic_form = "smooth";
end
if nargin < 45 || isempty(spline_group_kink_coefficient)
    spline_group_kink_coefficient = 0;
end
if nargin < 46 || isempty(spline_group_penalty_scale_T)
    spline_group_penalty_scale_T = spline_group_penalty_scale;
end
if nargin < 47 || isempty(spline_group_penalty_scale_C)
    spline_group_penalty_scale_C = spline_group_penalty_scale;
end
if nargin < 48 || isempty(spline_group_mcp_concavity)
    spline_group_mcp_concavity = 3;
end
if nargin < 49 || isempty(spline_group_iht_threshold_T)
    spline_group_iht_threshold_T = 0.1;
end
if nargin < 50 || isempty(spline_group_iht_threshold_C)
    spline_group_iht_threshold_C = spline_group_iht_threshold_T;
end
if nargin < 51 || isempty(spline_group_iht_warmup_iterations)
    spline_group_iht_warmup_iterations = 10;
end
if nargin < 52 || isempty(post_selection_refit_type)
    post_selection_refit_type = "spline";
end
post_selection_refit_type = lower(string(post_selection_refit_type));
if nargin < 53 || isempty(spline_group_iht_continuation_iterations)
    spline_group_iht_continuation_iterations = 0;
end
if nargin < 54 || isempty(spline_group_iht_start_fraction)
    spline_group_iht_start_fraction = 1;
end
if nargin < 55 || isempty(mpl_step_options)
    mpl_step_options = default_mpl_step_options(mpl_step_scale);
end
if nargin < 56 || isempty(decomposable_pilot_options)
    decomposable_pilot_options = struct( ...
        'iterations', 0, 'linear_scale', 1, 'group_scale', 1, ...
        'start', "zero", 'weighting', "none", 'weight_gamma', 1, ...
        'group_weight_gamma', 1, 'weight_epsilon', 0.05, 'weight_cap', 100);
end
linear_mic_form = lower(string(linear_mic_form));
linear_l1_type = lower(string(linear_l1_type));
spline_group_mic_form = lower(string(spline_group_mic_form));
if spline_group_penalty_type == "group_iht" && threshold_spline ~= 0
    error('PLAPH:GroupIHTDoubleThreshold', ...
        'ThresholdSpline must equal zero when group-IHT is active.');
end
if spline_threshold_stage == "final"
    threshold_spline_iter = 0;
else
    threshold_spline_iter = threshold_spline;
end
mpl_diag = empty_mpl_diagnostic();
[~,p] = size(Z);
[~,qs] = size(B);
mpl_inference = empty_mpl_inference(p);
if nargin < 57 || isempty(linear_mask_T)
    linear_mask_T = true(p,1);
end
if nargin < 58 || isempty(linear_mask_C)
    linear_mask_C = true(p,1);
end
if nargin < 59 || isempty(group_mask_T)
    group_mask_T = true(q,1);
end
if nargin < 60 || isempty(group_mask_C)
    group_mask_C = true(q,1);
end
if nargin < 61 || isempty(post_selection_refit_max_iter)
    post_selection_refit_max_iter = 30;
end
if nargin < 62 || isempty(baseline_roughness_scale)
    baseline_roughness_scale = 0.1;
end
if nargin < 63 || isempty(mic_penalty_scale_T)
    mic_penalty_scale_T = mic_penalty_scale;
end
if nargin < 64 || isempty(mic_penalty_scale_C)
    mic_penalty_scale_C = mic_penalty_scale;
end
if nargin < 65 || isempty(mic_penalty_margin_specific)
    mic_penalty_margin_specific = false;
end
linear_mask_T = logical(linear_mask_T(:));
linear_mask_C = logical(linear_mask_C(:));
group_mask_T = logical(group_mask_T(:));
group_mask_C = logical(group_mask_C(:));
if numel(linear_mask_T) ~= p || numel(linear_mask_C) ~= p || ...
        numel(group_mask_T) ~= q || numel(group_mask_C) ~= q
    error('PLAPH:MarginMaskSize', ...
        'Margin-specific screening masks do not match the MPL working dimensions.');
end
spline_mask_T = repelem(group_mask_T,qs/q);
spline_mask_C = repelem(group_mask_C,qs/q);
ini_beta(~linear_mask_T) = 0;
ini_phi(~linear_mask_C) = 0;
ini_bspline_beta(~spline_mask_T) = 0;
ini_bspline_phi(~spline_mask_C) = 0;
pilot_enabled = decomposable_pilot_options.iterations > 0;
pilot_state = empty_decomposable_pilot_state(p, qs, q, pilot_enabled);
pilot_linear_weight_beta = ones(p,1);
pilot_linear_weight_phi = ones(p,1);
pilot_group_weight_T = ones(q,1);
pilot_group_weight_C = ones(q,1);
if pilot_enabled && lower(string(decomposable_pilot_options.weighting)) == "adaptive"
    pilot_linear_weight_beta = min((abs(ini_beta) + ...
        decomposable_pilot_options.weight_epsilon).^ ...
        (-decomposable_pilot_options.weight_gamma), ...
        decomposable_pilot_options.weight_cap);
    pilot_linear_weight_phi = min((abs(ini_phi) + ...
        decomposable_pilot_options.weight_epsilon).^ ...
        (-decomposable_pilot_options.weight_gamma), ...
        decomposable_pilot_options.weight_cap);
    pilot_group_weight_T = min((spline_effect_norms( ...
        B,ini_bspline_beta,qs/q,q) + ...
        decomposable_pilot_options.weight_epsilon).^ ...
        (-decomposable_pilot_options.group_weight_gamma), ...
        decomposable_pilot_options.weight_cap);
    pilot_group_weight_C = min((spline_effect_norms( ...
        B,ini_bspline_phi,qs/q,q) + ...
        decomposable_pilot_options.weight_epsilon).^ ...
        (-decomposable_pilot_options.group_weight_gamma), ...
        decomposable_pilot_options.weight_cap);
end
pilot_state.linear_weight_beta = pilot_linear_weight_beta;
pilot_state.linear_weight_phi = pilot_linear_weight_phi;
pilot_state.group_weight_T = pilot_group_weight_T;
pilot_state.group_weight_C = pilot_group_weight_C;
if nargin < 20 || isempty(max_outer_iter)
    max_outer_iter = 300;
end

% Flexible number of time bins:
% This avoids requiring n0 to exactly divide n.
m = ceil(n/n0);

% tau = 0.2;  % % Kendall's tau
alph = copula_param_from_tau(copula_family,tau);  % copula dependence parameter
tau_current = tau;              % starting value for copula estimation
update_tau_every = 5;           % update tau every few outer iterations for stability
mpl_diag.mpl_profile_tau_enabled = estimate_tau;

beta = ini_beta;
phi = ini_phi;
if pilot_enabled && lower(string(decomposable_pilot_options.start)) == "zero"
    beta(:) = 0;
    phi(:) = 0;
    ini_bspline_beta(:) = 0;
    ini_bspline_phi(:) = 0;
end
[linear_l1_weights_beta,linear_l1_weights_phi] = linear_l1_weights( ...
    ini_beta,ini_phi,linear_l1_type,adaptive_weight_gamma,adaptive_weight_epsilon,adaptive_weight_cap);

beta1 = zeros(p,1);
phi1 = zeros(p,1);

bspline_beta = ini_bspline_beta;
bspline_phi = ini_bspline_phi;

lambda0_T = mic_penalty_lambda(sum(Delta),n,mic_penalty_scale_T, ...
    mic_penalty_normalization);  % MIC/SELO, Su et al. (2016)
if mic_penalty_margin_specific
    lambda0_C = mic_penalty_lambda(sum(1-Delta),n,mic_penalty_scale_C, ...
        mic_penalty_normalization);
else
    lambda0_C = lambda0_T;
end
mpl_diag.mpl_linear_mic_lambda_T = lambda0_T;
mpl_diag.mpl_linear_mic_lambda_C = lambda0_C;
% v2: optional per-coordinate penalty factor (0 = unpenalized, e.g. clinical
% adjustment covariates in the application). Scalar lambdas are kept when
% every factor equals one, so simulation fits are unaffected.
if nargin >= 66 && ~isempty(linear_penalty_factor) && any(linear_penalty_factor(:) ~= 1)
    if linear_mic_form ~= "kinked"
        error('PLAPH:PenaltyFactorForm', ...
            'Per-coordinate penalty factors are implemented for LinearMICForm="kinked" only.');
    end
    linear_penalty_factor = linear_penalty_factor(:);
    if numel(linear_penalty_factor) ~= p
        error('PLAPH:PenaltyFactorSize','Linear penalty factor must have length p.');
    end
    lambda0_T = lambda0_T*linear_penalty_factor;
    lambda0_C = lambda0_C*linear_penalty_factor;
    pilot_linear_weight_beta = pilot_linear_weight_beta.*linear_penalty_factor;
    pilot_linear_weight_phi = pilot_linear_weight_phi.*linear_penalty_factor;
end
spline_lambda_T = spline_group_penalty_lambda(sum(Delta), n, spline_group_penalty_scale_T, spline_group_penalty_normalization);
spline_lambda_C = spline_group_penalty_lambda(sum(1-Delta), n, spline_group_penalty_scale_C, spline_group_penalty_normalization);
% lambda1 = 0.5;
% lambda2 = 0.5;
lambda1 = baseline_roughness_scale*sqrt(n);
lambda2 = baseline_roughness_scale*sqrt(n);


%% Flexible discretization of ordered survival times
% The first m-1 bins contain approximately n0 observations,
% and the last bin contains the remaining observations.
[ID,binedg,binwv] = discrBinNA(T,n0);
m = max(ID);

psix = zeros(n,m);
Psi = zeros(n,m);

%% Piecewise constant basis function for the non-parametric baseline hazard
for i = 1:n
    psix(i,ID(i)) = 1;
end

%% Cumulative piecewise constant basis function
for i = 1:n
    current_bin = ID(i);
    if current_bin > 1
        Psi(i,1:(current_bin-1)) = binwv(1:(current_bin-1));
    end
    Psi(i,current_bin) = max(T(i)-binedg(current_bin),0);
end

%% Initial estimates of theta (piecewise constant estimate of h_{0t}) based on independent censoring assumption
ini_seta = sum( repmat(Delta,1,m).*psix )./( sum(repmat(exp(Z*ini_beta+B*ini_bspline_beta),1,m).*Psi)+1e-6 );

%% Initial estimates of gamma (piecewise constant estimate of h_{0c}) based on independent censoring assumption
ini_gama = sum(repmat(1-Delta,1,m).*psix)./( sum(repmat(exp(Z*ini_phi+B*ini_bspline_phi),1,m).*Psi)+1e-6 );
seta = ini_seta+1e-8;  % 1*m
gama = ini_gama+1e-8;  % 1*m

if collect_diagnostics
    start_linear_eta_T = Z*ini_beta;
    start_spline_eta_T = B*ini_bspline_beta;
    start_linear_eta_C = Z*ini_phi;
    start_spline_eta_C = B*ini_bspline_phi;
    start_eta_T = start_linear_eta_T + start_spline_eta_T;
    start_eta_C = start_linear_eta_C + start_spline_eta_C;
    mpl_diag.mpl_start_etaT_min = min(start_eta_T);
    mpl_diag.mpl_start_etaT_max = max(start_eta_T);
    mpl_diag.mpl_start_etaC_min = min(start_eta_C);
    mpl_diag.mpl_start_etaC_max = max(start_eta_C);
    mpl_diag.mpl_start_exp_etaT_max = max(safe_exp(start_eta_T));
    mpl_diag.mpl_start_exp_etaC_max = max(safe_exp(start_eta_C));
    mpl_diag.mpl_start_linear_etaT_max_abs = max(abs(start_linear_eta_T));
    mpl_diag.mpl_start_spline_etaT_max_abs = max(abs(start_spline_eta_T));
    mpl_diag.mpl_start_linear_etaC_max_abs = max(abs(start_linear_eta_C));
    mpl_diag.mpl_start_spline_etaC_max_abs = max(abs(start_spline_eta_C));
end

% % revised on 26/02/2025
% seta = max(seta, 1e-8);
% gama = max(gama, 1e-8);

RT = mat1(psix,T,1e-5); % First order difference -- For piecewise constant


%% Set up block coordinate descent parameters
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% Construct block indices
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


%% Set up block coordinate descent parameters--spline
num_basis = qs/q;
block_size_spline = block_size*num_basis;  % num_basis  % number of spline basis functions
num_blocks_spline = ceil(q / block_size);
blocks_spline = cell(num_blocks_spline, 1);

%% Construct block indices
for b = 1:num_blocks_spline
    start_idx1 = (b-1)*block_size_spline + 1;
    end_idx1 = min(b*block_size_spline, q*num_basis);
    blocks_spline{b} = start_idx1:end_idx1;
end


% %% Weighted update revised on 07/03/2025
% in_interval1 = (X(:,1) > 1.5) & (X(:,1) < -1.5);
% Weight1 = ones(n) * 1;  
% Weight1(~in_interval1) = 10;   
% Weight1 = diag(Weight1);
% Delta = Delta.*Weight1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1; err1=0; err2=0; err3=0; err4=0; err5=0; err6=0;
% tk1= 4; tk2= 200;
% tk1 = 100; tk2 = 200; tk3 = 100; tk4 = 100;
tk1 = 200*mpl_step_scale; tk2 = 200*mpl_step_scale;
tk3 = 100*mpl_step_scale; tk4 = 100*mpl_step_scale;
step_strategy = lower(string(mpl_step_options.strategy));
tk1_state = tk1; tk2_state = tk2; tk3_state = tk3; tk4_state = tk4;
mpl_diag.mpl_step_strategy = step_strategy;
warned_baseline_armijo = false;
group_iht_enabled = spline_group_penalty_type == "group_iht";
spline_group_metric = [];
if group_iht_enabled || pilot_enabled
    spline_group_metric = prepare_spline_group_metric(B, qs/q, q);
end
pilot_rank = max(qs/q-1, 1);
pilot_linear_lambda_T = decomposable_pilot_options.linear_scale * ...
    sqrt(max(sum(Delta),1)*log(max(p,2)))/n;
pilot_linear_lambda_C = decomposable_pilot_options.linear_scale * ...
    sqrt(max(sum(1-Delta),1)*log(max(p,2)))/n;
pilot_group_lambda_T = decomposable_pilot_options.group_scale * ...
    sqrt(max(sum(Delta),1)*(pilot_rank+log(max(q,2))))/n;
pilot_group_lambda_C = decomposable_pilot_options.group_scale * ...
    sqrt(max(sum(1-Delta),1)*(pilot_rank+log(max(q,2))))/n;
pilot_linear_penalty_T = pilot_linear_lambda_T*pilot_linear_weight_beta;
pilot_linear_penalty_C = pilot_linear_lambda_C*pilot_linear_weight_phi;
pilot_group_penalty_T = pilot_group_lambda_T*pilot_group_weight_T;
pilot_group_penalty_C = pilot_group_lambda_C*pilot_group_weight_C;
pilot_objective_min_delta = Inf;
target_warmup_iterations = decomposable_pilot_options.iterations + ...
    spline_group_iht_warmup_iterations;
mpl_diag.mpl_decomposable_pilot_enabled = pilot_enabled;
mpl_diag.mpl_decomposable_pilot_iterations = decomposable_pilot_options.iterations;
mpl_diag.mpl_decomposable_pilot_linear_lambda_T = pilot_linear_lambda_T;
mpl_diag.mpl_decomposable_pilot_linear_lambda_C = pilot_linear_lambda_C;
mpl_diag.mpl_decomposable_pilot_group_lambda_T = pilot_group_lambda_T;
mpl_diag.mpl_decomposable_pilot_group_lambda_C = pilot_group_lambda_C;
iht_lambda_T_final = 0.5 * tk3 * spline_group_iht_threshold_T^2;
iht_lambda_C_final = 0.5 * tk4 * spline_group_iht_threshold_C^2;
iht_lambda_T = iht_lambda_T_final;
iht_lambda_C = iht_lambda_C_final;
iht_objective_min_delta = Inf;

% v2: the iteration budget applies to the target (K-MIC/IHT) stage only;
% pilot iterations are added on top of max_outer_iter.
max_total_iter = max_outer_iter;
if pilot_enabled
    max_total_iter = max_outer_iter + decomposable_pilot_options.iterations;
end
gate_applicable = group_iht_enabled && linear_mic_form == "kinked";
if gate_applicable && (linear_l1_scale ~= 0 || linear_ridge_scale ~= 0 || ...
        step_strategy ~= "fixed" || spline_group_iht_continuation_iterations > 0)
    error('PLAPH:GateScope', ['The v2 selection gate assumes no external L1/ridge ' ...
        'terms, MPLStepStrategy="fixed" and no IHT threshold continuation.']);
end
selection_gate_trace=zeros(0,8);
selection_gate=struct('status','not_evaluated','pass',~gate_applicable,'maximum',NaN);
% v2: an objective increase no longer terminates the sparse stage. The
% iterate is restored and all block curvatures are doubled (step halving).
% The group-l0 penalty levels are held fixed, so the IHT thresholds become
% kappa*sqrt(t_base/t), as stated in the manuscript.
tk1_base = tk1; tk2_base = tk2; tk3_base = tk3; tk4_base = tk4;
step_doublings = 0;
max_step_doublings = 40;
% v2 acceleration in the target stage: the K-MIC fixed points do not depend
% on the curvature, and for group IHT only the retention DECISION depends on
% it. Linear blocks and retained spline groups therefore use adaptive
% curvatures (shrunk by 1.5 after an accepted sweep, never above the decision
% curvature, floor = decision/100); retention is always decided with the
% decision curvature tk3/tk4, so the set of fixed points is unchanged.
tk1_acc = tk1; tk2_acc = tk2; tk3_acc = tk3; tk4_acc = tk4;
accelerate_target = group_iht_enabled;

iht_support_changes = 0;
iht_support_stable_iterations = 0;
mpl_diag.mpl_group_iht_enabled = group_iht_enabled;
if group_iht_enabled
    mpl_diag.mpl_group_iht_warmup_iterations = spline_group_iht_warmup_iterations;
    mpl_diag.mpl_group_iht_continuation_iterations = spline_group_iht_continuation_iterations;
    mpl_diag.mpl_group_iht_start_fraction = spline_group_iht_start_fraction;
    mpl_diag.mpl_group_iht_threshold_T = spline_group_iht_threshold_T;
    mpl_diag.mpl_group_iht_threshold_C = spline_group_iht_threshold_C;
    mpl_diag.mpl_group_iht_lambda_T = iht_lambda_T;
    mpl_diag.mpl_group_iht_lambda_C = iht_lambda_C;
end

% beta = ini_beta; phi = ini_phi;
while k<=max_total_iter && ((err1==0 ||err2==0 || err3==0 || err4==0 || err5==0 || err6==0) || ...
        (pilot_enabled && k <= decomposable_pilot_options.iterations) || ...
        (group_iht_enabled && k <= target_warmup_iterations + ...
        spline_group_iht_continuation_iterations + 1))
    %k
    pilot_active_iteration = pilot_enabled && k <= decomposable_pilot_options.iterations;
    iht_active_iteration = group_iht_enabled && k > target_warmup_iterations;
    first_target_linear_iteration = pilot_enabled && ...
        k == decomposable_pilot_options.iterations + 1;
    first_target_group_iteration = pilot_enabled && group_iht_enabled && ...
        k == target_warmup_iterations + 1;
    spline_penalty_type_iter = spline_group_penalty_type;
    spline_lambda_T_iter = spline_lambda_T;
    spline_lambda_C_iter = spline_lambda_C;
    if pilot_active_iteration
        spline_penalty_type_iter = "group_lasso_metric";
        spline_lambda_T_iter = pilot_group_penalty_T;
        spline_lambda_C_iter = pilot_group_penalty_C;
    elseif group_iht_enabled && ~iht_active_iteration
        spline_penalty_type_iter = "group_mic";
    end
    spline_group_iht_threshold_T_iter = spline_group_iht_threshold_T*sqrt(tk3_base/tk3);
    spline_group_iht_threshold_C_iter = spline_group_iht_threshold_C*sqrt(tk4_base/tk4);
    if iht_active_iteration && spline_group_iht_continuation_iterations > 0
        iht_stage = k - target_warmup_iterations;
        progress = min(max((iht_stage-1)/spline_group_iht_continuation_iterations, 0), 1);
        threshold_fraction = spline_group_iht_start_fraction + ...
            (1-spline_group_iht_start_fraction)*progress;
        spline_group_iht_threshold_T_iter = threshold_fraction*spline_group_iht_threshold_T*sqrt(tk3_base/tk3);
        spline_group_iht_threshold_C_iter = threshold_fraction*spline_group_iht_threshold_C*sqrt(tk4_base/tk4);
    end
    iht_lambda_T = 0.5 * tk3 * spline_group_iht_threshold_T_iter^2;
    iht_lambda_C = 0.5 * tk4 * spline_group_iht_threshold_C_iter^2;
    beta_old_outer = beta;
    phi_old_outer = phi;
    bspline_beta_old_outer = bspline_beta;
    bspline_phi_old_outer = bspline_phi;
    seta_old_outer = seta;
    gama_old_outer = gama;
    L_outer_old = mpl_loglik_alpha(copula_family, alph, Z, B, Delta, ...
        beta, phi, bspline_beta, bspline_phi, seta, gama, Psi, psix);
    if pilot_active_iteration
        pilot_objective_old = mpl_decomposable_pilot_objective( ...
            copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
            seta,gama,Psi,psix,RT,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
            pilot_group_penalty_T,pilot_group_penalty_C,lambda1,lambda2,qs/q,q);
        if isnan(mpl_diag.mpl_decomposable_pilot_objective_initial)
            mpl_diag.mpl_decomposable_pilot_objective_initial = pilot_objective_old;
        end
    end
    if iht_active_iteration
        iht_objective_old = mpl_iht_penalized_objective( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, RT, a, lambda0_T,lambda0_C,linear_mic_form, ...
            linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
            linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
            iht_lambda_T, iht_lambda_C, qs/q, q);
        if isnan(mpl_diag.mpl_group_iht_objective_initial) && ...
                spline_group_iht_threshold_T_iter == spline_group_iht_threshold_T && ...
                spline_group_iht_threshold_C_iter == spline_group_iht_threshold_C
            mpl_diag.mpl_group_iht_objective_initial = iht_objective_old;
        end
    end

    %% Update regression coefficients
    H_T = sum( safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi),2 );  % Cumulative Risk - Failure
    H_C = sum( safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi),2 );   % Cumulative Risk - Censoring
    H_T = min(max(H_T, 0), 700);
    H_C = min(max(H_C, 0), 700);

    S_T = exp(-H_T);  % Survival Probability - Failure
    S_C = exp(-H_C);  % Survival Probability - Censoring

    [K1,K2,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic(copula_family,S_T,S_C,alph,Delta);
    DL_beta = -(sum( (Delta-Delta.*H_T-LAM1.*S_T.*H_T).*Z,1 ) )'/n;
    DL_beta(~isfinite(DL_beta)) = 0;
    DL_beta(~linear_mask_T) = 0;

 %% Select the update strategy according to the dimension
 if step_strategy ~= "fixed"
     [mic_weight_beta,l1_coefficient_beta,ridge_scale_beta] = ...
         linear_penalty_terms_for_iteration( ...
         beta,a,lambda0_T,linear_mic_form,linear_kink_coefficient, ...
         linear_l1_scale,linear_l1_weights_beta,linear_ridge_scale, ...
         pilot_active_iteration,pilot_linear_penalty_T);
     beta_loss = @(candidate) -mpl_loglik_alpha( ...
         copula_family,alph,Z,B,Delta,candidate,phi,bspline_beta, ...
         bspline_phi,seta,gama,Psi,psix);
     [beta1,tk1_state,step_backtracks,step_fallback] = ...
         adaptive_linear_prox_update(beta,DL_beta,mic_weight_beta, ...
         l1_coefficient_beta,ridge_scale_beta,tk1,tk1_state, ...
         mpl_step_options,beta_loss);
     mpl_diag.mpl_step_backtracks_beta = ...
         mpl_diag.mpl_step_backtracks_beta+step_backtracks;
     mpl_diag.mpl_step_fallbacks_beta = ...
         mpl_diag.mpl_step_fallbacks_beta+double(step_fallback);
 elseif p<=30
     [mic_weight_beta,l1_coefficient_beta,ridge_scale_beta] = ...
         linear_penalty_terms_for_iteration( ...
         beta,a,lambda0_T,linear_mic_form,linear_kink_coefficient, ...
         linear_l1_scale,linear_l1_weights_beta,linear_ridge_scale, ...
         pilot_active_iteration,pilot_linear_penalty_T);
     tk1_now = tk1; if accelerate_target && iht_active_iteration, tk1_now = tk1_acc; end
     beta_tilde = beta - DL_beta/tk1_now;
     beta_tilde(~isfinite(beta_tilde)) = 0;
     beta1 = hybrid_linear_prox(beta_tilde,mic_weight_beta,tk1_now,l1_coefficient_beta, ...
         ridge_scale_beta);
     beta1(~isfinite(beta1)) = 0;
 else
     % DL_beta and beta stay fixed inside the former 50-pass block loop,
     % and the blocks are disjoint. One pass therefore gives exactly the
     % same beta1. Consume the unused permutations to preserve the legacy
     % random-number stream for checkpoint-compatible formal runs.
     block_order1 = randperm(num_blocks);
     for b1 = block_order1
         current_block1 = blocks{b1};
         beta_block = beta(current_block1);
         [mic_weight_beta,l1_coefficient_beta,ridge_scale_beta] = ...
             linear_penalty_terms_for_iteration( ...
             beta_block,a,subset_penalty(lambda0_T,current_block1),linear_mic_form,linear_kink_coefficient, ...
             linear_l1_scale,linear_l1_weights_beta(current_block1), ...
             linear_ridge_scale,pilot_active_iteration, ...
             pilot_linear_penalty_T(current_block1));
         tk1_now = tk1; if accelerate_target && iht_active_iteration, tk1_now = tk1_acc; end
         beta_tilde = beta_block - DL_beta(current_block1)/tk1_now;
         beta_new1 = hybrid_linear_prox(beta_tilde,mic_weight_beta,tk1_now, ...
             l1_coefficient_beta,ridge_scale_beta);
         beta1(current_block1,1) = beta_new1;
     end
     delta1 = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
     if delta1 >= 1e-6
         for ite = 2:50
             randperm(num_blocks);
         end
     end
 end
        
    beta1(~linear_mask_T) = 0;
    if first_target_linear_iteration
        target_beta_L = tk1;
        if step_strategy ~= "fixed"
            target_beta_L = tk1_state;
        end
        [~,target_l1_beta,~] = linear_penalty_terms_for_iteration( ...
            beta_old_outer,a,lambda0_T,linear_mic_form,linear_kink_coefficient, ...
            linear_l1_scale,linear_l1_weights_beta,linear_ridge_scale,false,0);
        pilot_state.first_target_beta_candidate = ...
            beta_old_outer-DL_beta/target_beta_L;
        pilot_state.first_target_beta_threshold = target_l1_beta/target_beta_L;
    end
    w1 = beta1-beta;
    err1 = norm(w1,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;

    % Revised on 31/12/2024
    H_T = sum( safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi),2 );  % Cumulative Risk - Failure
    H_C = sum( safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi),2 );   % Cumulative Risk - Censoring
    H_T = min(max(H_T, 0), 700);
    H_C = min(max(H_C, 0), 700);

    S_T = exp(-H_T);  % Survival Probability - Failure
    S_C = exp(-H_C);  % Survival Probability - Censoring

    [K1,K2,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic(copula_family,S_T,S_C,alph,Delta);

    DL_phi = -(sum( ((1-Delta)-(1-Delta).*H_C-LAM2.*S_C.*H_C).*Z,1 ))'/n;
    DL_phi(~isfinite(DL_phi)) = 0;
    DL_phi(~linear_mask_C) = 0;

    if step_strategy ~= "fixed"
        [mic_weight_phi,l1_coefficient_phi,ridge_scale_phi] = ...
            linear_penalty_terms_for_iteration( ...
            phi,a,lambda0_C,linear_mic_form,linear_kink_coefficient, ...
            linear_l1_scale,linear_l1_weights_phi,linear_ridge_scale, ...
            pilot_active_iteration,pilot_linear_penalty_C);
        phi_loss = @(candidate) -mpl_loglik_alpha( ...
            copula_family,alph,Z,B,Delta,beta,candidate,bspline_beta, ...
            bspline_phi,seta,gama,Psi,psix);
        [phi1,tk2_state,step_backtracks,step_fallback] = ...
            adaptive_linear_prox_update(phi,DL_phi,mic_weight_phi, ...
            l1_coefficient_phi,ridge_scale_phi,tk2,tk2_state, ...
            mpl_step_options,phi_loss);
        mpl_diag.mpl_step_backtracks_phi = ...
            mpl_diag.mpl_step_backtracks_phi+step_backtracks;
        mpl_diag.mpl_step_fallbacks_phi = ...
            mpl_diag.mpl_step_fallbacks_phi+double(step_fallback);
    elseif p<=30
        [mic_weight_phi,l1_coefficient_phi,ridge_scale_phi] = ...
            linear_penalty_terms_for_iteration( ...
            phi,a,lambda0_C,linear_mic_form,linear_kink_coefficient, ...
            linear_l1_scale,linear_l1_weights_phi,linear_ridge_scale, ...
            pilot_active_iteration,pilot_linear_penalty_C);
        tk2_now = tk2; if accelerate_target && iht_active_iteration, tk2_now = tk2_acc; end
        phi_tilde = phi - DL_phi/tk2_now;
        phi_tilde(~isfinite(phi_tilde)) = 0;
        phi1 = hybrid_linear_prox(phi_tilde,mic_weight_phi,tk2_now,l1_coefficient_phi, ...
            ridge_scale_phi);
        phi1(~isfinite(phi1)) = 0;
    else
        block_order2 = randperm(num_blocks);
        for b2 = block_order2
            current_block2 = blocks{b2};
            phi_block = phi(current_block2);
            [mic_weight_phi,l1_coefficient_phi,ridge_scale_phi] = ...
                linear_penalty_terms_for_iteration( ...
                phi_block,a,subset_penalty(lambda0_C,current_block2),linear_mic_form,linear_kink_coefficient, ...
                linear_l1_scale,linear_l1_weights_phi(current_block2), ...
                linear_ridge_scale,pilot_active_iteration, ...
                pilot_linear_penalty_C(current_block2));
            tk2_now = tk2; if accelerate_target && iht_active_iteration, tk2_now = tk2_acc; end
            phi_tilde = phi_block - DL_phi(current_block2)/tk2_now;
            phi_new1 = hybrid_linear_prox(phi_tilde,mic_weight_phi,tk2_now, ...
                l1_coefficient_phi,ridge_scale_phi);
            phi1(current_block2,1) = phi_new1;
        end
        delta2 = norm(phi - phi1,2)^2 / (norm(phi,2)^2 + 1e-8);
        if delta2 >= 1e-6
            for it = 2:50
                randperm(num_blocks);
            end
        end
    end

    phi1(~linear_mask_C) = 0;
    if first_target_linear_iteration
        target_phi_L = tk2;
        if step_strategy ~= "fixed"
            target_phi_L = tk2_state;
        end
        [~,target_l1_phi,~] = linear_penalty_terms_for_iteration( ...
            phi_old_outer,a,lambda0_C,linear_mic_form,linear_kink_coefficient, ...
            linear_l1_scale,linear_l1_weights_phi,linear_ridge_scale,false,0);
        pilot_state.first_target_phi_candidate = ...
            phi_old_outer-DL_phi/target_phi_L;
        pilot_state.first_target_phi_threshold = target_l1_phi/target_phi_L;
    end
    w2 = phi1-phi;
    err2 = norm(w2,2)^2 <= r*norm(phi,2)^2;
    phi = phi1;


    %% Update spline coefficients
    H_T = sum( safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi),2 );  % Cumulative Risk - Failure
    H_C = sum( safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi),2 );   % Cumulative Risk - Censoring
    H_T = min(max(H_T, 0), 700);
    H_C = min(max(H_C, 0), 700);

    S_T = exp(-H_T);  % Survival Probability - Failure
    S_C = exp(-H_C);  % Survival Probability - Censoring

    [K1,K2,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic(copula_family,S_T,S_C,alph,Delta);


    DL_beta_spline =  -(sum( (Delta-Delta.*H_T-LAM1.*S_T.*H_T).*B,1) )'/n;  % first partial derivative - beta
    DL_beta_spline(~isfinite(DL_beta_spline)) = 0;
    DL_beta_spline(~spline_mask_T) = 0;
    % bspline_beta1 = bspline_beta - DL_beta_spline/(tk1/2);
    if spline_penalty_type_iter == "group_lasso_metric"
        bspline_beta1 = apply_spline_group_lasso_metric( ...
            spline_group_metric,bspline_beta,DL_beta_spline,qs/q,q,tk3, ...
            spline_lambda_T_iter);
    elseif spline_penalty_type_iter == "group_iht" && step_strategy ~= "fixed"
        stage_lambda_T = 0.5*tk3*spline_group_iht_threshold_T_iter^2;
        spline_T_loss = @(candidate) -mpl_loglik_alpha( ...
            copula_family,alph,Z,B,Delta,beta,phi,candidate, ...
            bspline_phi,seta,gama,Psi,psix);
        [bspline_beta1,tk3_state,step_backtracks,step_fallback] = ...
            adaptive_group_iht_update(spline_group_metric,bspline_beta, ...
            DL_beta_spline,qs/q,q,tk3,tk3_state,stage_lambda_T, ...
            mpl_step_options,spline_T_loss);
        mpl_diag.mpl_step_backtracks_spline_T = ...
            mpl_diag.mpl_step_backtracks_spline_T+step_backtracks;
        mpl_diag.mpl_step_fallbacks_spline_T = ...
            mpl_diag.mpl_step_fallbacks_spline_T+double(step_fallback);
    elseif spline_penalty_type_iter == "group_iht"
        bspline_beta1 = apply_spline_group_iht_metric_acc( ...
            spline_group_metric, bspline_beta, DL_beta_spline, qs/q, q, tk3, ...
            spline_group_iht_threshold_T_iter, tk3_acc);
    elseif q<=q_size
        bspline_beta1 = bspline_beta - DL_beta_spline/tk3;
    else
        block_order_spline = randperm(num_blocks_spline);
        for b = block_order_spline
            current_block_spline = blocks_spline{b};
            bspline_beta_block = bspline_beta(current_block_spline);
            bspline_beta_new1 = bspline_beta_block - ...
                DL_beta_spline(current_block_spline)/tk3;
            bspline_beta1(current_block_spline,1) = bspline_beta_new1;
        end
        delta = norm(bspline_beta - bspline_beta1,2)^2 / ...
            (norm(bspline_beta,2)^2 + 1e-8);
        if delta >= 1e-6
            for iter = 2:50
                randperm(num_blocks_spline);
            end
        end
    end

    bspline_beta1(~spline_mask_T) = 0;
    if first_target_group_iteration
        target_group_T_L = tk3;
        if step_strategy ~= "fixed"
            target_group_T_L = tk3_state;
        end
        pilot_state.first_target_group_T_candidate_norm = ...
            spline_group_candidate_effect_norms(spline_group_metric, ...
            bspline_beta_old_outer,DL_beta_spline,qs/q,q,target_group_T_L);
        pilot_state.first_target_group_T_threshold = ...
            sqrt(max(2*iht_lambda_T/target_group_T_L,0));
    end
    if ~any(spline_penalty_type_iter == ["group_iht" "group_lasso_metric"])
        bspline_beta1 = apply_spline_group_penalty_step( ...
            bspline_beta1, qs/q, q, spline_lambda_T_iter, tk3, spline_penalty_type_iter, ...
            spline_group_mic_sharpness, spline_group_mic_continuation, spline_group_mic_sharpness_path, ...
            spline_group_mic_form, spline_group_kink_coefficient, spline_group_mcp_concavity, ...
            spline_group_iht_threshold_T);
    end
    bspline_beta1 = apply_spline_group_threshold_metric( ...
        B, bspline_beta1, qs/q, q, threshold_spline_iter, spline_selection_metric);

    w3 = bspline_beta1-bspline_beta;
    err3 = norm(w3,2)^2 <= r*norm(bspline_beta,2)^2;
    bspline_beta = bspline_beta1;

    DL_phi_spline =   -(sum(((1-Delta)-(1-Delta).*H_C-LAM2.*S_C.*H_C).*B,1))'/n;  % first partial derivative - phi
    DL_phi_spline(~isfinite(DL_phi_spline)) = 0;
    DL_phi_spline(~spline_mask_C) = 0;
    % bspline_phi1 = bspline_phi - DL_phi_spline/(tk2/2);

    if spline_penalty_type_iter == "group_lasso_metric"
        bspline_phi1 = apply_spline_group_lasso_metric( ...
            spline_group_metric,bspline_phi,DL_phi_spline,qs/q,q,tk4, ...
            spline_lambda_C_iter);
    elseif spline_penalty_type_iter == "group_iht" && step_strategy ~= "fixed"
        stage_lambda_C = 0.5*tk4*spline_group_iht_threshold_C_iter^2;
        spline_C_loss = @(candidate) -mpl_loglik_alpha( ...
            copula_family,alph,Z,B,Delta,beta,phi,bspline_beta, ...
            candidate,seta,gama,Psi,psix);
        [bspline_phi1,tk4_state,step_backtracks,step_fallback] = ...
            adaptive_group_iht_update(spline_group_metric,bspline_phi, ...
            DL_phi_spline,qs/q,q,tk4,tk4_state,stage_lambda_C, ...
            mpl_step_options,spline_C_loss);
        mpl_diag.mpl_step_backtracks_spline_C = ...
            mpl_diag.mpl_step_backtracks_spline_C+step_backtracks;
        mpl_diag.mpl_step_fallbacks_spline_C = ...
            mpl_diag.mpl_step_fallbacks_spline_C+double(step_fallback);
    elseif spline_penalty_type_iter == "group_iht"
        bspline_phi1 = apply_spline_group_iht_metric_acc( ...
            spline_group_metric, bspline_phi, DL_phi_spline, qs/q, q, tk4, ...
            spline_group_iht_threshold_C_iter, tk4_acc);
    elseif q<=q_size
        bspline_phi1 = bspline_phi - DL_phi_spline/tk4;
    else
        block_order_spline = randperm(num_blocks_spline);
        for b = block_order_spline
            current_block_spline = blocks_spline{b};
            bspline_phi_block = bspline_phi(current_block_spline);
            bspline_phi_new1 = bspline_phi_block - ...
                DL_phi_spline(current_block_spline)/tk4;
            bspline_phi1(current_block_spline,1) = bspline_phi_new1;
        end
        delta = norm(bspline_phi - bspline_phi1,2)^2 / ...
            (norm(bspline_phi,2)^2 + 1e-8);
        if delta >= 1e-6
            for iter = 2:50
                randperm(num_blocks_spline);
            end
        end
    end

    bspline_phi1(~spline_mask_C) = 0;
    if first_target_group_iteration
        target_group_C_L = tk4;
        if step_strategy ~= "fixed"
            target_group_C_L = tk4_state;
        end
        pilot_state.first_target_group_C_candidate_norm = ...
            spline_group_candidate_effect_norms(spline_group_metric, ...
            bspline_phi_old_outer,DL_phi_spline,qs/q,q,target_group_C_L);
        pilot_state.first_target_group_C_threshold = ...
            sqrt(max(2*iht_lambda_C/target_group_C_L,0));
    end
    if ~any(spline_penalty_type_iter == ["group_iht" "group_lasso_metric"])
        bspline_phi1 = apply_spline_group_penalty_step( ...
            bspline_phi1, qs/q, q, spline_lambda_C_iter, tk4, spline_penalty_type_iter, ...
            spline_group_mic_sharpness, spline_group_mic_continuation, spline_group_mic_sharpness_path, ...
            spline_group_mic_form, spline_group_kink_coefficient, spline_group_mcp_concavity, ...
            spline_group_iht_threshold_C);
    end
    bspline_phi1 = apply_spline_group_threshold_metric( ...
        B, bspline_phi1, qs/q, q, threshold_spline_iter, spline_selection_metric);

    w4 = bspline_phi1-bspline_phi;
    err4 = norm(w4,2)^2 <= r*norm(bspline_phi,2)^2;
    bspline_phi = bspline_phi1;


    %% Robust Armijo line search for baseline hazards
    % This guarded version caps the number of backtracking steps and falls back
    % to the previous baseline hazards if no acceptable step is found.

    uu = 0.5;
    sigm = 0.5;
    max_armijo_iter = 30;
    uu_min = 1e-10;

    H_T0 = sum( safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi),2 );
    H_C0 = sum( safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi),2 );
    H_T0 = min(max(H_T0, 0), 700);
    H_C0 = min(max(H_C0, 0), 700);
    S_T0 = exp(-H_T0);
    S_C0 = exp(-H_C0);

    [K1_0,K2_0,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic(copula_family,S_T0,S_C0,alph,Delta);

    DL_seta = sum( (Delta.*psix)./(sum(seta.*psix,2)+1e-8) - ...
        (safe_exp(Z*beta+B*bspline_beta).*(Delta+LAM1.*S_T0)).*Psi,1)/n - lambda1*seta*RT;

    DL_gama = sum( ((1-Delta).*psix)./(sum(gama.*psix,2)+1e-8) - ...
        (safe_exp(Z*phi+B*bspline_phi).*(1-Delta+LAM2.*S_C0)).*Psi,1)/n - lambda2*gama*RT;

    Ps1 = sum((Delta.*safe_exp(Z*beta+B*bspline_beta)).*Psi + ...
        ((LAM11.*S_T0).*safe_exp(Z*beta+B*bspline_beta)).*Psi,1)/n + ...
        lambda1*max(seta*RT,0) + 1e-6;

    Ps2 = sum(((1-Delta).*safe_exp(Z*phi+B*bspline_phi)).*Psi + ...
        ((LAM22.*S_C0).*safe_exp(Z*phi+B*bspline_phi)).*Psi,1)/n + ...
        lambda2*max(gama*RT,0) + 1e-6;

    L_T0 = log(sum(seta.*psix,2)+1e-8) + Z*beta + B*bspline_beta - H_T0 + log(K1_0+1e-8);
    L_C0 = log(sum(gama.*psix,2)+1e-8) + Z*phi  + B*bspline_phi  - H_C0 + log(K2_0+1e-8);
    L_0 = mean(Delta.*L_T0 + (1-Delta).*L_C0);
    if pilot_active_iteration
        pilot_objective_before_baseline = mpl_decomposable_pilot_objective( ...
            copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
            seta,gama,Psi,psix,RT,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
            pilot_group_penalty_T,pilot_group_penalty_C,lambda1,lambda2,qs/q,q);
    end
    if iht_active_iteration
        iht_objective_before_baseline = mpl_iht_penalized_objective( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, RT, a, lambda0_T,lambda0_C,linear_mic_form, ...
            linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
            linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
            iht_lambda_T, iht_lambda_C, qs/q, q);
    end

    if collect_diagnostics && k == 1
        precond_seta_step = DL_seta.*(seta./Ps1);
        precond_gama_step = DL_gama.*(gama./Ps2);
        mpl_diag.mpl_diag_collected = true;
        mpl_diag.mpl_first_L_before_coef = L_outer_old;
        mpl_diag.mpl_first_L_after_coef = L_0;
        mpl_diag.mpl_first_coef_delta = L_0 - L_outer_old;
        mpl_diag.mpl_first_DL_seta_norm = norm(DL_seta);
        mpl_diag.mpl_first_DL_gama_norm = norm(DL_gama);
        mpl_diag.mpl_first_Ps1_min = min(Ps1);
        mpl_diag.mpl_first_Ps2_min = min(Ps2);
        mpl_diag.mpl_first_Ps1_nonpos = sum(Ps1 <= 0);
        mpl_diag.mpl_first_Ps2_nonpos = sum(Ps2 <= 0);
        mpl_diag.mpl_first_seta_min = min(seta);
        mpl_diag.mpl_first_gama_min = min(gama);
        mpl_diag.mpl_first_HT_max = max(H_T0);
        mpl_diag.mpl_first_HC_max = max(H_C0);
        mpl_diag.mpl_first_ST_min = min(S_T0);
        mpl_diag.mpl_first_SC_min = min(S_C0);
        eta_T_diag = Z*beta + B*bspline_beta;
        eta_C_diag = Z*phi + B*bspline_phi;
        mpl_diag.mpl_first_etaT_min = min(eta_T_diag);
        mpl_diag.mpl_first_etaT_max = max(eta_T_diag);
        mpl_diag.mpl_first_etaC_min = min(eta_C_diag);
        mpl_diag.mpl_first_etaC_max = max(eta_C_diag);
        mpl_diag.mpl_first_exp_etaT_max = max(safe_exp(eta_T_diag));
        mpl_diag.mpl_first_exp_etaC_max = max(safe_exp(eta_C_diag));
        mpl_diag.mpl_first_K1_min = min(K1_0);
        mpl_diag.mpl_first_K2_min = min(K2_0);
        mpl_diag.mpl_first_LAM1_abs_max = max(abs(LAM1));
        mpl_diag.mpl_first_LAM2_abs_max = max(abs(LAM2));
        mpl_diag.mpl_first_LAM11_abs_max = max(abs(LAM11));
        mpl_diag.mpl_first_LAM22_abs_max = max(abs(LAM22));
        mpl_diag.mpl_first_precond_plus_1e4 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, precond_seta_step, precond_gama_step, 1e-4, L_0);
        mpl_diag.mpl_first_precond_minus_1e4 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, precond_seta_step, precond_gama_step, -1e-4, L_0);
        mpl_diag.mpl_first_raw_plus_1e4 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, DL_seta, DL_gama, 1e-4, L_0);
        mpl_diag.mpl_first_raw_minus_1e4 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, DL_seta, DL_gama, -1e-4, L_0);
        mpl_diag.mpl_first_precond_plus_1e2 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, precond_seta_step, precond_gama_step, 1e-2, L_0);
        mpl_diag.mpl_first_precond_minus_1e2 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, precond_seta_step, precond_gama_step, -1e-2, L_0);
        mpl_diag.mpl_first_raw_plus_1e2 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, DL_seta, DL_gama, 1e-2, L_0);
        mpl_diag.mpl_first_raw_minus_1e2 = baseline_step_delta( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, DL_seta, DL_gama, -1e-2, L_0);
    end

    seta_old = seta;
    gama_old = gama;
    accepted_step = false;

    baseline_directions = {
        DL_seta.*(seta_old./Ps1), DL_gama.*(gama_old./Ps2);
        DL_seta,                 DL_gama
    };
    for direction_mode = 1:size(baseline_directions, 1)
        step_seta = baseline_directions{direction_mode, 1};
        step_gama = baseline_directions{direction_mode, 2};
        for direction_sign = [1 -1]
            uu = 0.5;
            for armijo_iter = 1:max_armijo_iter
                seta_trial = max(seta_old + direction_sign*uu*step_seta, 1e-8);
                gama_trial = max(gama_old + direction_sign*uu*step_gama, 1e-8);

                H_T_trial = sum( safe_exp(Z*beta+B*bspline_beta).*(seta_trial.*Psi),2 );
                H_C_trial = sum( safe_exp(Z*phi+B*bspline_phi).*(gama_trial.*Psi),2 );
                H_T_trial = min(max(H_T_trial, 0), 700);
                H_C_trial = min(max(H_C_trial, 0), 700);
                S_T_trial = exp(-H_T_trial);
                S_C_trial = exp(-H_C_trial);

                [K1_trial,K2_trial,~,~,~,~] = copula_terms_generic(copula_family,S_T_trial,S_C_trial,alph,Delta);

                L_T_trial = log(sum(seta_trial.*psix,2)+1e-8) + Z*beta + B*bspline_beta - H_T_trial + log(K1_trial+1e-8);
                L_C_trial = log(sum(gama_trial.*psix,2)+1e-8) + Z*phi  + B*bspline_phi  - H_C_trial + log(K2_trial+1e-8);
                L_new = mean(Delta.*L_T_trial + (1-Delta).*L_C_trial);

                accept_baseline_trial = isfinite(L_new) && (L_new >= L_0 - 1e-8);
                if pilot_active_iteration
                    pilot_objective_trial = mpl_decomposable_pilot_objective( ...
                        copula_family,alph,Z,B,Delta,beta,phi, ...
                        bspline_beta,bspline_phi,seta_trial,gama_trial, ...
                        Psi,psix,RT,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
                        pilot_group_penalty_T,pilot_group_penalty_C, ...
                        lambda1,lambda2,qs/q,q);
                    accept_baseline_trial = isfinite(pilot_objective_trial) && ...
                        pilot_objective_trial >= pilot_objective_before_baseline - 1e-8;
                end
                if iht_active_iteration
                    iht_objective_trial = mpl_iht_penalized_objective( ...
                        copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
                        seta_trial, gama_trial, Psi, psix, RT, a, ...
                        lambda0_T,lambda0_C,linear_mic_form, ...
                        linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
                        linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
                        iht_lambda_T, iht_lambda_C, qs/q, q);
                    accept_baseline_trial = isfinite(iht_objective_trial) && ...
                        iht_objective_trial >= iht_objective_before_baseline - 1e-8;
                end
                if accept_baseline_trial
                    accepted_step = true;
                    break;
                end

                uu = uu * sigm;
                if uu < uu_min
                    break;
                end
            end
            if accepted_step
                break;
            end
        end
        if accepted_step
            break;
        end
    end

    if accepted_step
        seta = seta_trial;
        gama = gama_trial;
    else
        seta = seta_old;
        gama = gama_old;
        accept_coef_only = isfinite(L_0) && isfinite(L_outer_old) && L_0 >= L_outer_old - 1e-8;
        if pilot_active_iteration
            pilot_objective_coef = mpl_decomposable_pilot_objective( ...
                copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
                seta,gama,Psi,psix,RT,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
                pilot_group_penalty_T,pilot_group_penalty_C,lambda1,lambda2,qs/q,q);
            accept_coef_only = isfinite(pilot_objective_coef) && ...
                isfinite(pilot_objective_old) && ...
                pilot_objective_coef >= pilot_objective_old - 1e-8;
        end
        if iht_active_iteration
            iht_objective_coef = mpl_iht_penalized_objective( ...
                copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
                seta, gama, Psi, psix, RT, a,lambda0_T,lambda0_C, ...
                linear_mic_form, ...
                linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
                linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
                iht_lambda_T, iht_lambda_C, qs/q, q);
            accept_coef_only = isfinite(iht_objective_coef) && isfinite(iht_objective_old) && ...
                iht_objective_coef >= iht_objective_old - 1e-8;
        end
        if accept_coef_only
            if ~warned_baseline_armijo
                warning('Baseline Armijo failed; accepting regression/spline updates with unchanged baseline hazards in this MPL call.');
                warned_baseline_armijo = true;
            end
            mpl_diag.mpl_baseline_armijo_failed = true;
        else
            beta = beta_old_outer;
            phi = phi_old_outer;
            bspline_beta = bspline_beta_old_outer;
            bspline_phi = bspline_phi_old_outer;
            seta = seta_old_outer;
            gama = gama_old_outer;
            if accelerate_target && iht_active_iteration && ...
                    (tk1_acc < tk1 || tk2_acc < tk2 || tk3_acc < tk3 || tk4_acc < tk4)
                tk1_acc = min(2*tk1_acc,tk1); tk2_acc = min(2*tk2_acc,tk2);
                tk3_acc = min(2*tk3_acc,tk3); tk4_acc = min(2*tk4_acc,tk4);
                err1 = 0;
                continue;
            end
            [tk1,tk2,tk3,tk4,step_doublings] = plaph_double_curvatures( ...
                tk1,tk2,tk3,tk4,step_doublings);
            tk1_acc = tk1; tk2_acc = tk2; tk3_acc = tk3; tk4_acc = tk4;
            if step_doublings > max_step_doublings
                warning('Armijo backtracking exhausted at outer iteration %d; rejecting the MPL outer update.', k);
                mpl_diag.mpl_armijo_rejected = true;
                break;
            end
            err1 = 0;
            continue;
        end
    end

    if pilot_active_iteration
        pilot_objective_new = mpl_decomposable_pilot_objective( ...
            copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
            seta,gama,Psi,psix,RT,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
            pilot_group_penalty_T,pilot_group_penalty_C,lambda1,lambda2,qs/q,q);
        pilot_objective_delta = pilot_objective_new-pilot_objective_old;
        pilot_objective_min_delta = min(pilot_objective_min_delta,pilot_objective_delta);
        if ~isfinite(pilot_objective_new) || pilot_objective_delta < -1e-8
            beta = beta_old_outer;
            phi = phi_old_outer;
            bspline_beta = bspline_beta_old_outer;
            bspline_phi = bspline_phi_old_outer;
            seta = seta_old_outer;
            gama = gama_old_outer;
            if accelerate_target && iht_active_iteration && ...
                    (tk1_acc < tk1 || tk2_acc < tk2 || tk3_acc < tk3 || tk4_acc < tk4)
                tk1_acc = min(2*tk1_acc,tk1); tk2_acc = min(2*tk2_acc,tk2);
                tk3_acc = min(2*tk3_acc,tk3); tk4_acc = min(2*tk4_acc,tk4);
                err1 = 0;
                continue;
            end
            [tk1,tk2,tk3,tk4,step_doublings] = plaph_double_curvatures( ...
                tk1,tk2,tk3,tk4,step_doublings);
            tk1_acc = tk1; tk2_acc = tk2; tk3_acc = tk3; tk4_acc = tk4;
            if step_doublings > max_step_doublings
                mpl_diag.mpl_decomposable_pilot_objective_rejected = true;
                break;
            end
            err1 = 0;
            continue;
        end
        pilot_state.iterations = k;
        pilot_state.beta = beta;
        pilot_state.phi = phi;
        pilot_state.bspline_beta = bspline_beta;
        pilot_state.bspline_phi = bspline_phi;
        if k == decomposable_pilot_options.iterations
            [pilot_state.stationarity_beta,pilot_state.stationarity_phi, ...
                pilot_state.stationarity_group_T,pilot_state.stationarity_group_C] = ...
                decomposable_pilot_stationarity( ...
                copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
                seta,gama,Psi,pilot_linear_penalty_T,pilot_linear_penalty_C, ...
                pilot_group_penalty_T,pilot_group_penalty_C,spline_group_metric, ...
                qs/q,q,tk1,tk2,tk3,tk4);
            pilot_state.completed = true;
            mpl_diag.mpl_decomposable_pilot_completed = true;
            mpl_diag.mpl_decomposable_pilot_objective_final = pilot_objective_new;
            mpl_diag.mpl_decomposable_pilot_stationarity_beta = ...
                pilot_state.stationarity_beta;
            mpl_diag.mpl_decomposable_pilot_stationarity_phi = ...
                pilot_state.stationarity_phi;
            mpl_diag.mpl_decomposable_pilot_stationarity_group_T = ...
                pilot_state.stationarity_group_T;
            mpl_diag.mpl_decomposable_pilot_stationarity_group_C = ...
                pilot_state.stationarity_group_C;
        end
    end

    if iht_active_iteration
        iht_objective_new = mpl_iht_penalized_objective( ...
            copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
            seta, gama, Psi, psix, RT, a,lambda0_T,lambda0_C, ...
            linear_mic_form, ...
            linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
            linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
            iht_lambda_T, iht_lambda_C, qs/q, q);
        iht_objective_delta = iht_objective_new - iht_objective_old;
        iht_objective_min_delta = min(iht_objective_min_delta, iht_objective_delta);
        if ~isfinite(iht_objective_new) || iht_objective_delta < -1e-8
            beta = beta_old_outer;
            phi = phi_old_outer;
            bspline_beta = bspline_beta_old_outer;
            bspline_phi = bspline_phi_old_outer;
            seta = seta_old_outer;
            gama = gama_old_outer;
            if accelerate_target && iht_active_iteration && ...
                    (tk1_acc < tk1 || tk2_acc < tk2 || tk3_acc < tk3 || tk4_acc < tk4)
                tk1_acc = min(2*tk1_acc,tk1); tk2_acc = min(2*tk2_acc,tk2);
                tk3_acc = min(2*tk3_acc,tk3); tk4_acc = min(2*tk4_acc,tk4);
                err1 = 0;
                continue;
            end
            [tk1,tk2,tk3,tk4,step_doublings] = plaph_double_curvatures( ...
                tk1,tk2,tk3,tk4,step_doublings);
            tk1_acc = tk1; tk2_acc = tk2; tk3_acc = tk3; tk4_acc = tk4;
            if step_doublings > max_step_doublings
                mpl_diag.mpl_group_iht_objective_rejected = true;
                break;
            end
            err1 = 0;
            continue;
        end

        support_T_old = spline_group_support(bspline_beta_old_outer, qs/q, q);
        support_C_old = spline_group_support(bspline_phi_old_outer, qs/q, q);
        support_T_new = spline_group_support(bspline_beta, qs/q, q);
        support_C_new = spline_group_support(bspline_phi, qs/q, q);
        support_change_count = sum(xor(support_T_old, support_T_new)) + ...
            sum(xor(support_C_old, support_C_new));
        iht_support_changes = iht_support_changes + support_change_count;
        if support_change_count == 0
            iht_support_stable_iterations = iht_support_stable_iterations + 1;
        else
            iht_support_stable_iterations = 0;
        end
    end

    err5 = norm(seta-seta_old,2)^2 <= r*(norm(seta_old,2)^2 + 1e-8);
    err6 = norm(gama-gama_old,2)^2 <= r*(norm(gama_old,2)^2 + 1e-8);

    % Update the copula parameter by profiling the current observed-data log-likelihood.
    if estimate_tau && mod(k, update_tau_every) == 0
        mpl_diag.mpl_profile_tau_attempts = ...
            mpl_diag.mpl_profile_tau_attempts+1;
        try
            profile_loglik_old = mpl_loglik_alpha(copula_family,alph,Z,B, ...
                Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix);
            [tau_new, alph_new] = update_copula_tau( ...
                copula_family, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, seta, gama, Psi, psix);
            profile_loglik_new = mpl_loglik_alpha(copula_family,alph_new,Z,B, ...
                Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix);
            profile_delta = profile_loglik_new-profile_loglik_old;
            if isnan(mpl_diag.mpl_profile_tau_min_loglik_delta)
                mpl_diag.mpl_profile_tau_min_loglik_delta = profile_delta;
            else
                mpl_diag.mpl_profile_tau_min_loglik_delta = min( ...
                    mpl_diag.mpl_profile_tau_min_loglik_delta,profile_delta);
            end
            if isfinite(tau_new) && isfinite(alph_new) && ...
                    isfinite(profile_loglik_new) && ...
                    isfinite(profile_loglik_old) && profile_delta >= -1e-10
                tau_current = tau_new;
                alph = alph_new;
                mpl_diag.mpl_profile_tau_accepted = ...
                    mpl_diag.mpl_profile_tau_accepted+1;
                mpl_diag.mpl_profile_tau_boundary_hit = ...
                    mpl_diag.mpl_profile_tau_boundary_hit || ...
                    profile_tau_boundary_hit(copula_family,tau_new);
            else
                mpl_diag.mpl_profile_tau_rejected = ...
                    mpl_diag.mpl_profile_tau_rejected+1;
            end
        catch ME
            mpl_diag.mpl_profile_tau_failures = ...
                mpl_diag.mpl_profile_tau_failures+1;
            warning('Copula parameter update failed at outer iteration %d: %s. Keeping previous value.', k, ME.message);
        end
    end

    if accelerate_target && iht_active_iteration
        tk1_acc = max(tk1_acc/1.5,tk1/100); tk2_acc = max(tk2_acc/1.5,tk2/100);
        tk3_acc = max(tk3_acc/1.5,tk3/100); tk4_acc = max(tk4_acc/1.5,tk4/100);
    end
    small_selection_change=err1 && err2 && err3 && err4 && err5 && err6;
    if gate_applicable && k>target_warmup_iterations && ...
            (small_selection_change || mod(k,25)==0 || k==max_total_iter)
        selection_gate=selection_gate_v2(copula_family,alph,Z,B,Delta,beta,phi, ...
            bspline_beta,bspline_phi,seta,gama,Psi,psix,qs/q,q,spline_group_metric, ...
            tk3,tk4,spline_group_iht_threshold_T*sqrt(tk3_base/tk3), ...
            spline_group_iht_threshold_C*sqrt(tk4_base/tk4), ...
            lambda0_T,lambda0_C,a,linear_kink_coefficient,linear_mask_T,linear_mask_C, ...
            spline_mask_T,spline_mask_C,RT,lambda1,lambda2);
        selection_gate_trace(end+1,:)=[k,selection_gate.linear_T,selection_gate.linear_C, ...
            selection_gate.baseline_T,selection_gate.baseline_C, ...
            selection_gate.group_T,selection_gate.group_C,selection_gate.pass]; %#ok<AGROW>
        if small_selection_change && ~selection_gate.pass
            err5=0;
        end
    end

    mpl_diag.mpl_outer_iterations = k;
    k = k+1;
end

% Passive diagnostic snapshot: never used by optimization.
pilot_state.diagnostic_pre_refit_state = struct( ...
    'beta',beta,'phi',phi,'bspline_beta',bspline_beta, ...
    'bspline_phi',bspline_phi,'seta',seta,'gama',gama,'alph',alph, ...
    'q',q,'num_basis',qs/q,'tk3',tk3,'tk4',tk4, ...
    'threshold_T',spline_group_iht_threshold_T, ...
    'threshold_C',spline_group_iht_threshold_C, ...
    'fit_copula',copula_family,'outer_iterations',mpl_diag.mpl_outer_iterations);

if gate_applicable
    selection_gate=selection_gate_v2(copula_family,alph,Z,B,Delta,beta,phi, ...
        bspline_beta,bspline_phi,seta,gama,Psi,psix,qs/q,q,spline_group_metric, ...
        tk3,tk4,spline_group_iht_threshold_T*sqrt(tk3_base/tk3), ...
        spline_group_iht_threshold_C*sqrt(tk4_base/tk4), ...
        lambda0_T,lambda0_C,a,linear_kink_coefficient,linear_mask_T,linear_mask_C, ...
        spline_mask_T,spline_mask_C,RT,lambda1,lambda2);
end
selection_status_code = 4;
if mpl_diag.mpl_armijo_rejected || mpl_diag.mpl_decomposable_pilot_objective_rejected || ...
        mpl_diag.mpl_group_iht_objective_rejected
    selection_gate.status='update_rejected';
    selection_status_code = 3;
elseif selection_gate.pass
    selection_gate.status='necessary_checks_pass';
    selection_status_code = 1;
elseif k>max_total_iter
    selection_gate.status='iteration_budget_exhausted';
    selection_status_code = 2;
else
    selection_gate.status='stopped_without_necessary_checks';
end
mpl_diag.mpl_selection_gate_pass = logical(selection_gate.pass) && selection_status_code == 1;
mpl_diag.mpl_selection_gate_maximum = selection_gate.maximum;
mpl_diag.mpl_selection_status_code = selection_status_code;
mpl_diag.mpl_step_doublings = step_doublings;
mpl_diag.mpl_target_iterations = max(k-1-decomposable_pilot_options.iterations*pilot_enabled,0);
pilot_state.diagnostic_selection_gate=selection_gate;
pilot_state.diagnostic_selection_gate_trace=selection_gate_trace;
mpl_diag.mpl_step_final_L_beta = tk1_state;
mpl_diag.mpl_step_final_L_phi = tk2_state;
mpl_diag.mpl_step_final_L_spline_T = tk3_state;
mpl_diag.mpl_step_final_L_spline_C = tk4_state;
if pilot_enabled && isfinite(pilot_objective_min_delta)
    mpl_diag.mpl_decomposable_pilot_objective_min_delta = ...
        pilot_objective_min_delta;
end

%% Restore the final group-L0 weights after the finite continuation stage.
iht_lambda_T = iht_lambda_T_final;
iht_lambda_C = iht_lambda_C_final;

% v2: a single numerical-zero tolerance is applied once, before the refit;
% all reported supports are taken from the frozen masks below.
beta2 = beta.*(abs(beta)>=threshold_linear_beta);
phi2 = phi.*(abs(phi)>=threshold_linear_phi);
if nargin >= 66 && ~isempty(linear_penalty_factor) && any(linear_penalty_factor(:) == 0)
    % Unpenalized adjustment covariates are always in the refit model.
    keep = linear_penalty_factor(:) == 0;
    beta2(keep) = beta(keep) + (beta(keep) == 0)*1e-8;
    phi2(keep) = phi(keep) + (phi(keep) == 0)*1e-8;
end
opt_seta = seta;
opt_gama = gama;

if spline_threshold_stage == "final"
    bspline_beta = apply_spline_group_threshold_metric( ...
        B, bspline_beta, qs/q, q, threshold_spline, spline_selection_metric);
    bspline_phi = apply_spline_group_threshold_metric( ...
        B, bspline_phi, qs/q, q, threshold_spline, spline_selection_metric);
end

% Freeze the selected model before the unpenalized refit. The refit updates
% only these indices, so reporting must not re-threshold refitted magnitudes.
pilot_state.final_selected_linear_T_mask = abs(beta2) > 0;
pilot_state.final_selected_linear_C_mask = abs(phi2) > 0;
pilot_state.final_selected_additive_T_mask = ...
    spline_group_support(bspline_beta,qs/q,q);
pilot_state.final_selected_additive_C_mask = ...
    spline_group_support(bspline_phi,qs/q,q);

if post_selection_refit_type == "joint_coefficients"
    beta = beta2;
    phi = phi2;
    [beta, phi, bspline_beta, bspline_phi, refit_diag] = ...
        post_selection_joint_coefficient_refit_MPL( ...
        copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
        seta, gama, Psi, psix, post_selection_refit_max_iter, ...
        max(tk1,100), max(tk2,100), ...
        max(tk3,100), max(tk4,100), q);
    mpl_diag.mpl_joint_refit_enabled = true;
    mpl_diag.mpl_joint_refit_initial_loglik = refit_diag.initial_loglik;
    mpl_diag.mpl_joint_refit_final_loglik = refit_diag.final_loglik;
    mpl_diag.mpl_joint_refit_iterations = refit_diag.iterations;
    mpl_diag.mpl_joint_refit_rejected = refit_diag.rejected;
elseif any(post_selection_refit_type == ["joint_oracle" "joint_oracle_alpha"])
    beta = beta2;
    phi = phi2;
    [beta, phi, bspline_beta, bspline_phi, seta, gama, alph, refit_diag] = ...
        post_selection_joint_oracle_refit_MPL( ...
        copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
        seta, gama, Psi, psix, post_selection_refit_max_iter, ...
        max(tk1,100), max(tk2,100), ...
        max(tk3,100), max(tk4,100), q, ...
        post_selection_refit_type == "joint_oracle_alpha",mpl_step_options);
    mpl_diag.mpl_joint_refit_enabled = true;
    mpl_diag.mpl_joint_refit_initial_loglik = refit_diag.initial_loglik;
    mpl_diag.mpl_joint_refit_final_loglik = refit_diag.final_loglik;
    mpl_diag.mpl_joint_refit_iterations = refit_diag.iterations;
    mpl_diag.mpl_joint_refit_rejected = refit_diag.rejected;
    mpl_diag.mpl_oracle_refit_enabled = true;
    mpl_diag.mpl_oracle_refit_converged = refit_diag.converged;
    mpl_diag.mpl_oracle_refit_coefficient_rejections = refit_diag.coefficient_rejections;
    mpl_diag.mpl_oracle_refit_baseline_rejections = refit_diag.baseline_rejections;
    mpl_diag.mpl_oracle_refit_alpha_updates = refit_diag.alpha_updates;
    mpl_diag.mpl_oracle_refit_alpha_rejections = refit_diag.alpha_rejections;
    mpl_diag.mpl_oracle_refit_alpha_boundary_hit = ...
        isfinite(refit_diag.final_tau) && abs(refit_diag.final_tau) >= 0.94;
    mpl_diag.mpl_oracle_refit_coefficient_backtracks = ...
        refit_diag.coefficient_backtracks;
    mpl_diag.mpl_oracle_refit_baseline_backtracks = ...
        refit_diag.baseline_backtracks;
    mpl_diag.mpl_oracle_refit_final_coefficient_step = ...
        refit_diag.final_coefficient_step;
    mpl_diag.mpl_oracle_refit_final_baseline_step = ...
        refit_diag.final_baseline_step;
    mpl_diag.mpl_oracle_refit_min_loglik_delta = refit_diag.min_loglik_delta;
    mpl_diag.mpl_oracle_refit_final_relative_change = refit_diag.final_relative_change;
    mpl_diag.mpl_oracle_refit_final_relative_loglik_change = ...
        refit_diag.final_relative_loglik_change;
    mpl_diag.mpl_oracle_refit_final_coefficient_score_inf = ...
        refit_diag.final_coefficient_score_inf;
    mpl_diag.mpl_oracle_refit_coefficient_score_inf_raw = ...
        refit_diag.final_coefficient_score_inf_raw;
    mpl_diag.mpl_oracle_refit_final_baseline_kkt_inf = ...
        refit_diag.final_baseline_kkt_inf;
    mpl_diag.mpl_oracle_refit_final_stationarity_inf = ...
        refit_diag.final_stationarity_inf;
    mpl_diag.mpl_oracle_refit_stationarity_pass = ...
        refit_diag.stationarity_pass;
    mpl_diag.mpl_oracle_refit_stagnation_handoff = ...
        refit_diag.stagnation_handoff;
    mpl_diag.mpl_oracle_refit_fmincon_attempted = ...
        refit_diag.fmincon_polish_attempted;
    mpl_diag.mpl_oracle_refit_fmincon_accepted = ...
        refit_diag.fmincon_polish_accepted;
    mpl_diag.mpl_oracle_refit_fmincon_exitflag = refit_diag.fmincon_exitflag;
    mpl_diag.mpl_oracle_refit_fmincon_iterations = refit_diag.fmincon_iterations;
    mpl_diag.mpl_oracle_refit_fmincon_firstorderopt = ...
        refit_diag.fmincon_firstorderopt;
    mpl_diag.mpl_oracle_refit_gradient_check_max_abs = ...
        refit_diag.gradient_check_max_abs;
    mpl_diag.mpl_oracle_refit_gradient_check_max_rel = ...
        refit_diag.gradient_check_max_rel;
    mpl_diag.mpl_oracle_refit_gradient_check_max_index = ...
        refit_diag.gradient_check_max_index;
    mpl_diag.mpl_oracle_refit_gradient_check_analytic = ...
        refit_diag.gradient_check_analytic;
    mpl_diag.mpl_oracle_refit_gradient_check_numeric = ...
        refit_diag.gradient_check_numeric;
    mpl_diag.mpl_oracle_refit_gradient_check_block = ...
        refit_diag.gradient_check_block;
    if mpl_step_options.collect_model_based_inference
        if post_selection_refit_type == "joint_oracle" && ...
                refit_diag.stationarity_pass
            mpl_inference = fixed_support_model_inference_MPL( ...
                copula_family,alph,Z,B,Delta,beta,phi,bspline_beta, ...
                bspline_phi,seta,gama,Psi,psix, ...
                find(pilot_state.final_selected_linear_T_mask), ...
                find(pilot_state.final_selected_linear_C_mask), ...
                spline_group_columns(find( ...
                    pilot_state.final_selected_additive_T_mask),qs/q), ...
                spline_group_columns(find( ...
                    pilot_state.final_selected_additive_C_mask),qs/q));
        elseif post_selection_refit_type ~= "joint_oracle"
            mpl_inference.attempted = true;
            mpl_inference.reason = "requires_fixed_copula_joint_oracle_refit";
        else
            mpl_inference.attempted = true;
            mpl_inference.reason = "fixed_support_refit_not_stationary";
        end
    end
    if post_selection_refit_type == "joint_oracle_alpha" && ...
            isfinite(refit_diag.final_tau)
        tau_current = refit_diag.final_tau;
    end
    opt_seta = seta;
    opt_gama = gama;
    H_T = sum(safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi), 2);
    H_C = sum(safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi), 2);
    H_T = min(max(H_T, 0), 700);
    H_C = min(max(H_C, 0), 700);
else
    % Archived behavior: condition on the penalized linear coefficients and
    % refit selected spline groups only.
    [bspline_beta, bspline_phi] = post_selection_spline_refit_MPL( ...
        copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
        seta, gama, Psi, 30, max(tk3,100), max(tk4,100), q);
end

mpl_diag.mpl_inference_attempted = mpl_inference.attempted;
mpl_diag.mpl_inference_valid = mpl_inference.valid;
mpl_diag.mpl_inference_reason = mpl_inference.reason;
mpl_diag.mpl_inference_dimension = mpl_inference.dimension;
mpl_diag.mpl_inference_rank = mpl_inference.rank;
mpl_diag.mpl_inference_condition_number = mpl_inference.condition_number;
mpl_diag.mpl_inference_min_eigenvalue = mpl_inference.min_eigenvalue;

final_mean_loglik = mpl_loglik_alpha(copula_family, alph, Z, B, Delta, ...
    beta, phi, bspline_beta, bspline_phi, seta, gama, Psi, psix);
mpl_diag.mpl_final_mean_loglik = final_mean_loglik;
mpl_diag.mpl_final_total_loglik = n * final_mean_loglik;
if group_iht_enabled
    mpl_diag.mpl_group_iht_objective_final = mpl_iht_penalized_objective( ...
        copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
        seta, gama, Psi, psix, RT, a,lambda0_T,lambda0_C, ...
        linear_mic_form, ...
        linear_kink_coefficient, linear_l1_scale, linear_l1_weights_beta, ...
        linear_l1_weights_phi, linear_ridge_scale, lambda1, lambda2, ...
        iht_lambda_T, iht_lambda_C, qs/q, q);
    if isfinite(iht_objective_min_delta)
        mpl_diag.mpl_group_iht_objective_min_delta = iht_objective_min_delta;
    end
    mpl_diag.mpl_group_iht_support_changes = iht_support_changes;
    mpl_diag.mpl_group_iht_support_stable_iterations = iht_support_stable_iterations;
end

opt_bspline_beta = bspline_beta;
opt_bspline_phi = bspline_phi;
if any(post_selection_refit_type == ["joint_coefficients" "joint_oracle" "joint_oracle_alpha"])
    opt_beta = beta;
    opt_phi = phi;
else
    opt_beta = beta2;
    opt_phi = phi2;
end

tau_hat = tau_current;
alph_hat = alph;
if estimate_tau
    mpl_diag.mpl_profile_tau_boundary_hit = ...
        mpl_diag.mpl_profile_tau_boundary_hit || ...
        profile_tau_boundary_hit(copula_family,tau_current);
end

% H_T = sum(opt_seta.*Psi,2);
% H_C = sum(opt_gama.*Psi,2);

end


%% ===================================================
%          Generic copula helper functions
% ====================================================
function theta = copula_param_from_tau(family,tau_val)
%% Convert Kendall's tau to the copula dependence parameter.
% Supported families: Frank, Gaussian, Clayton, and Gumbel.

    family = lower(family);
    switch family
        case 'frank'
            theta = copulaparam('Frank',tau_val);
        case 'gaussian'
            theta = sin(pi*tau_val/2);  % correlation parameter rho
        case 'clayton'
            tau_val = max(min(tau_val,0.95),1e-4);
            theta = 2*tau_val/(1-tau_val);
        case 'gumbel'
            tau_val = max(min(tau_val,0.95),1e-4);
            theta = 1/(1-tau_val);
        otherwise
            error('Unsupported copula family: %s',family);
    end
end


function U = sample_copula(family,theta,n)
%% Generate n observations from the selected bivariate copula.

    family_l = lower(family);
    switch family_l
        case 'frank'
            if abs(theta) < 1e-6
                U = rand(n,2);
            else
                U = copularnd('Frank',theta,n);
            end
        case 'gaussian'
            R = [1 theta; theta 1];
            U = copularnd('Gaussian',R,n);
        case 'clayton'
            U = copularnd('Clayton',theta,n);
        case 'gumbel'
            U = copularnd('Gumbel',theta,n);
        otherwise
            error('Unsupported copula family: %s',family);
    end
end


function C = copula_cdf_generic(family,u,v,theta)
%% Evaluate a bivariate copula C(u,v).

    eps0 = 1e-8;
    u = min(max(u,eps0),1-eps0);
    v = min(max(v,eps0),1-eps0);
    family_l = lower(family);

    switch family_l
        case 'frank'
            A = exp(-theta*u);
            B = exp(-theta*v);
            D = exp(-theta) - 1;
            if abs(theta) < 1e-6
                C = u.*v;
            else
                C = -log(1 + ((A-1).*(B-1))./D)./theta;
            end
        case 'gaussian'
            R = [1 theta; theta 1];
            C = copulacdf('Gaussian',[u(:) v(:)],R);
            C = reshape(C,size(u));
        case 'clayton'
            C = copulacdf('Clayton',[u(:) v(:)],theta);
            C = reshape(C,size(u));
        case 'gumbel'
            C = copulacdf('Gumbel',[u(:) v(:)],theta);
            C = reshape(C,size(u));
        otherwise
            error('Unsupported copula family: %s',family);
    end
end


function [K1,K2,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic(family,S_T,S_C,theta,Delta)
%% Copula derivative terms for the MPL score equations.
% Let u = S_T and v = S_C. For observed failure and censoring likelihood
% contributions, we need
%   K1 = C_u(u,v),  K2 = C_v(u,v),
% and the derivative ratios
%   LAM1 = Delta*C_uu/C_u + (1-Delta)*C_uv/C_v,
%   LAM2 = Delta*C_uv/C_u + (1-Delta)*C_vv/C_v.
% These terms enter the beta/phi/spline/baseline score equations.
% All supported families use analytic first derivatives and derivative
% ratios. This avoids nested finite differences near the survival boundaries,
% where the previous non-Frank implementation could produce inaccurate scores.

    eps0 = 1e-8;
    active_u = S_T > eps0 & S_T < 1-eps0;
    active_v = S_C > eps0 & S_C < 1-eps0;
    u = min(max(S_T,eps0),1-eps0);
    v = min(max(S_C,eps0),1-eps0);
    family_l = lower(family);

    [K1,K2] = copula_K_only(family_l,u,v,theta);

    switch family_l
        case 'frank'
            if abs(theta) < 1e-6
                dlogK1_du = zeros(size(u));
                dlogK2_du = 1./u;
                dlogK1_dv = 1./v;
                dlogK2_dv = zeros(size(v));
            else
                A = exp(-theta*u);
                B = exp(-theta*v);
                Aminus1 = signed_denominator(expm1(-theta*u),-sign(theta));
                Bminus1 = signed_denominator(expm1(-theta*v),-sign(theta));
                dlogK1_du = theta.*(K1-1);
                dlogK2_du = -theta.*A./Aminus1+theta.*K1;
                dlogK1_dv = -theta.*B./Bminus1+theta.*K2;
                dlogK2_dv = theta.*(K2-1);
            end
        case 'gaussian'
            rho = theta;
            scale = sqrt(max(1-rho^2,eps0));
            zu = norminv(u);
            zv = norminv(v);
            a = (zv-rho.*zu)./scale;
            b = (zu-rho.*zv)./scale;
            mills_a = normpdf(a)./max(K1,realmin('double'));
            mills_b = normpdf(b)./max(K2,realmin('double'));
            phi_u = max(normpdf(zu),realmin('double'));
            phi_v = max(normpdf(zv),realmin('double'));
            dlogK1_du = -rho.*mills_a./(scale.*phi_u);
            dlogK1_dv = mills_a./(scale.*phi_v);
            dlogK2_du = mills_b./(scale.*phi_u);
            dlogK2_dv = -rho.*mills_b./(scale.*phi_v);
        case 'clayton'
            A = u.^(-theta)+v.^(-theta)-1;
            common = theta+1;
            du = u.^(-theta-1)./A;
            dv = v.^(-theta-1)./A;
            dlogK1_du = common.*(du-1./u);
            dlogK1_dv = common.*dv;
            dlogK2_du = common.*du;
            dlogK2_dv = common.*(dv-1./v);
        case 'gumbel'
            x = max(-log(u),eps0);
            y = max(-log(v),eps0);
            A = x.^theta+y.^theta;
            common_x = A.^(1/theta-1).*x.^(theta-1)./u;
            common_y = A.^(1/theta-1).*y.^(theta-1)./v;
            curvature_x = (theta-1).*x.^(theta-1)./(u.*A);
            curvature_y = (theta-1).*y.^(theta-1)./(v.*A);
            dlogK1_du = common_x+curvature_x-(theta-1)./(u.*x)-1./u;
            dlogK1_dv = common_y+curvature_y;
            dlogK2_du = common_x+curvature_x;
            dlogK2_dv = common_y+curvature_y-(theta-1)./(v.*y)-1./v;
        otherwise
            error('Unsupported copula family: %s',family);
    end

    if ~strcmp(family_l,'frank')
        % The finite implementation evaluates log(K+eps0). Apply its exact
        % chain-rule factor so the analytic score matches that objective.
        weight1 = K1./(K1+eps0);
        weight2 = K2./(K2+eps0);
        dlogK1_du = weight1.*dlogK1_du;
        dlogK1_dv = weight1.*dlogK1_dv;
        dlogK2_du = weight2.*dlogK2_du;
        dlogK2_dv = weight2.*dlogK2_dv;
    end

    LAM1 = Delta.*dlogK1_du + (1-Delta).*dlogK2_du;
    LAM2 = Delta.*dlogK1_dv + (1-Delta).*dlogK2_dv;
    % The implemented likelihood clips survival probabilities before
    % evaluating K1 and K2. Outside the unclipped interval, the derivative
    % of that stabilized copula term with respect to the raw survival
    % probability is exactly zero.
    LAM1 = active_u.*LAM1;
    LAM2 = active_v.*LAM2;

    % The bounded survival arguments keep these ratios finite. Do not cap
    % them: a cap would change the score without changing the likelihood.

    % Nonnegative versions used only in the multiplicative baseline update.
    LAM11 = max(LAM1,0);
    LAM22 = max(LAM2,0);
end

function [K1,K2] = copula_K_only(family,u,v,theta)
%% Helper for first partial derivatives K1 and K2.

    eps0 = 1e-8;
    u = min(max(u,eps0),1-eps0);
    v = min(max(v,eps0),1-eps0);

    switch lower(family)
        case 'frank'
            if abs(theta) < 1e-6
                K1 = v;
                K2 = u;
                return;
            end
            A = exp(-theta*u);
            B = exp(-theta*v);
            Aminus1 = expm1(-theta*u);
            Bminus1 = expm1(-theta*v);
            D = signed_denominator(expm1(-theta)+Aminus1.*Bminus1,-sign(theta));
            K1 = max(A.*Bminus1./D,realmin('double'));
            K2 = max(B.*Aminus1./D,realmin('double'));
        case 'gaussian'
            rho = theta;
            scale = sqrt(max(1-rho^2,eps0));
            zu = norminv(u);
            zv = norminv(v);
            K1 = max(normcdf((zv-rho.*zu)./scale),realmin('double'));
            K2 = max(normcdf((zu-rho.*zv)./scale),realmin('double'));
        case 'clayton'
            A = u.^(-theta)+v.^(-theta)-1;
            K1 = max(A.^(-1/theta-1).*u.^(-theta-1),realmin('double'));
            K2 = max(A.^(-1/theta-1).*v.^(-theta-1),realmin('double'));
        case 'gumbel'
            x = max(-log(u),eps0);
            y = max(-log(v),eps0);
            A = x.^theta+y.^theta;
            radial = A.^(1/theta);
            common = exp(-radial).*A.^(1/theta-1);
            K1 = max(common.*x.^(theta-1)./u,realmin('double'));
            K2 = max(common.*y.^(theta-1)./v,realmin('double'));
        otherwise
            error('Unsupported copula family: %s',family);
    end
end


function [tau_hat, theta_hat] = update_copula_tau(family,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix)
%% Estimate Kendall's tau by profile likelihood for the selected copula.

    family_l = lower(family);
    if any(strcmp(family_l,{'clayton','gumbel'}))
        tau_lower = 1e-4;
    else
        tau_lower = -0.95;
    end
    tau_upper = 0.95;

    obj = @(tau_val) -copula_loglik_tau(family_l,tau_val,Z,B,Delta,beta,phi, ...
        bspline_beta,bspline_phi,seta,gama,Psi,psix);

    options = optimset('TolX',1e-4,'Display','off');
    tau_hat = fminbnd(obj, tau_lower, tau_upper, options);
    theta_hat = copula_param_from_tau(family_l,tau_hat);
end


function loglik = copula_loglik_tau(family,tau_val,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix)
%% Observed-data log-likelihood contribution as a function of Kendall's tau.

    theta = copula_param_from_tau(family,tau_val);

    eta_T = Z*beta + B*bspline_beta;
    eta_C = Z*phi  + B*bspline_phi;

    exp_eta_T = exp(min(max(eta_T,-700),700));
    exp_eta_C = exp(min(max(eta_C,-700),700));

    H_T = sum(exp_eta_T .* (seta .* Psi), 2);
    H_C = sum(exp_eta_C .* (gama .* Psi), 2);

    S_T = exp(-H_T);
    S_C = exp(-H_C);

    [K1,K2,~,~,~,~] = copula_terms_generic(family,S_T,S_C,theta,Delta);

    L_T = log(sum(seta.*psix,2)+1e-8) + eta_T - H_T + log(K1+1e-8);
    L_C = log(sum(gama.*psix,2)+1e-8) + eta_C - H_C + log(K2+1e-8);

    loglik = mean(Delta.*L_T + (1-Delta).*L_C);
end


function value = signed_denominator(value,default_sign)
% Preserve the analytic sign while protecting only machine-zero denominators.
tiny = realmin('double');
near_zero = abs(value) < tiny;
if any(near_zero,'all')
    replacement_sign = sign(value(near_zero));
    replacement_sign(replacement_sign == 0) = default_sign;
    value(near_zero) = replacement_sign.*tiny;
end
end


function hit = profile_tau_boundary_hit(family,tau_val)
family = lower(string(family));
if any(family == ["clayton" "gumbel"])
    hit = tau_val <= 0.01 || tau_val >= 0.94;
else
    hit = abs(tau_val) >= 0.94;
end
end


function loglik = mpl_loglik_alpha(family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix)
%% Observed-data MPL log-likelihood at a fixed copula parameter.
eta_T = Z*beta + B*bspline_beta;
eta_C = Z*phi  + B*bspline_phi;

loglik = mpl_loglik_from_eta( ...
    family,theta,Delta,eta_T,eta_C,seta,gama,Psi,psix);
end


function loglik = mpl_loglik_from_eta( ...
    family,theta,Delta,eta_T,eta_C,seta,gama,Psi,psix)
% Evaluate the same likelihood while reusing an available linear predictor.

exp_eta_T = safe_exp(eta_T);
exp_eta_C = safe_exp(eta_C);

H_T = sum(exp_eta_T .* (seta .* Psi), 2);
H_C = sum(exp_eta_C .* (gama .* Psi), 2);
H_T = min(max(H_T, 0), 700);
H_C = min(max(H_C, 0), 700);

S_T = exp(-H_T);
S_C = exp(-H_C);

[K1,K2,~,~,~,~] = copula_terms_generic(family,S_T,S_C,theta,Delta);

L_T = log(sum(seta.*psix,2)+1e-8) + eta_T - H_T + log(K1+1e-8);
L_C = log(sum(gama.*psix,2)+1e-8) + eta_C - H_C + log(K2+1e-8);

loglik = mean(Delta.*L_T + (1-Delta).*L_C);
end


function objective = mpl_decomposable_pilot_objective( ...
    family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama, ...
    Psi,psix,RT,linear_lambda_T,linear_lambda_C,group_lambda_T, ...
    group_lambda_C,baseline_lambda_T,baseline_lambda_C,num_basis,q)
loglik = mpl_loglik_alpha(family,theta,Z,B,Delta,beta,phi, ...
    bspline_beta,bspline_phi,seta,gama,Psi,psix);
linear_penalty = sum(linear_lambda_T(:).*abs(beta)) + ...
    sum(linear_lambda_C(:).*abs(phi));
group_penalty = sum(group_lambda_T(:).*spline_effect_norms( ...
    B,bspline_beta,num_basis,q)) + ...
    sum(group_lambda_C(:).*spline_effect_norms(B,bspline_phi,num_basis,q));
baseline_penalty = 0.5*baseline_lambda_T*(seta*RT*seta') + ...
    0.5*baseline_lambda_C*(gama*RT*gama');
objective = loglik-linear_penalty-group_penalty-baseline_penalty;
end


function [res_beta,res_phi,res_group_T,res_group_C] = ...
    decomposable_pilot_stationarity( ...
    family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama, ...
    Psi,linear_lambda_T,linear_lambda_C,group_lambda_T,group_lambda_C, ...
    metric,num_basis,q,L_beta,L_phi,L_group_T,L_group_C)
n = size(Z,1);
H_T = sum(safe_exp(Z*beta+B*bspline_beta).*(seta.*Psi),2);
H_C = sum(safe_exp(Z*phi+B*bspline_phi).*(gama.*Psi),2);
H_T = min(max(H_T,0),700);
H_C = min(max(H_C,0),700);
S_T = exp(-H_T);
S_C = exp(-H_C);
[~,~,LAM1,LAM2,~,~] = copula_terms_generic( ...
    family,S_T,S_C,theta,Delta);

gradient_beta = -(sum((Delta-Delta.*H_T-LAM1.*S_T.*H_T).*Z,1))'/n;
gradient_phi = -(sum(((1-Delta)-(1-Delta).*H_C- ...
    LAM2.*S_C.*H_C).*Z,1))'/n;
gradient_group_T = -(sum((Delta-Delta.*H_T- ...
    LAM1.*S_T.*H_T).*B,1))'/n;
gradient_group_C = -(sum(((1-Delta)-(1-Delta).*H_C- ...
    LAM2.*S_C.*H_C).*B,1))'/n;
gradient_beta(~isfinite(gradient_beta)) = 0;
gradient_phi(~isfinite(gradient_phi)) = 0;
gradient_group_T(~isfinite(gradient_group_T)) = 0;
gradient_group_C(~isfinite(gradient_group_C)) = 0;

if isscalar(linear_lambda_T)
    linear_lambda_T = linear_lambda_T*ones(size(beta));
end
if isscalar(linear_lambda_C)
    linear_lambda_C = linear_lambda_C*ones(size(phi));
end
prox_beta = hybrid_linear_prox(beta-gradient_beta/L_beta, ...
    zeros(size(beta)),L_beta,linear_lambda_T(:),0);
prox_phi = hybrid_linear_prox(phi-gradient_phi/L_phi, ...
    zeros(size(phi)),L_phi,linear_lambda_C(:),0);
prox_group_T = apply_spline_group_lasso_metric(metric,bspline_beta, ...
    gradient_group_T,num_basis,q,L_group_T,group_lambda_T);
prox_group_C = apply_spline_group_lasso_metric(metric,bspline_phi, ...
    gradient_group_C,num_basis,q,L_group_C,group_lambda_C);

res_beta = L_beta*max(abs(beta-prox_beta));
res_phi = L_phi*max(abs(phi-prox_phi));
res_group_T = L_group_T*max_spline_metric_effect_norm( ...
    metric,bspline_beta-prox_group_T,num_basis,q);
res_group_C = L_group_C*max_spline_metric_effect_norm( ...
    metric,bspline_phi-prox_group_C,num_basis,q);
end


function objective = mpl_iht_penalized_objective( ...
    family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix,RT, ...
    sharpness,lambda0_T,lambda0_C,mic_form,kink_coefficient,l1_scale,l1_weights_beta, ...
    l1_weights_phi,ridge_scale,baseline_lambda_T,baseline_lambda_C, ...
    iht_lambda_T,iht_lambda_C,num_basis,q)
% Mean observed log-likelihood minus the penalties used by the group-IHT
% candidate estimator. This check is isolated from the legacy group-MIC path.
loglik = mpl_loglik_alpha(family,theta,Z,B,Delta,beta,phi, ...
    bspline_beta,bspline_phi,seta,gama,Psi,psix);

linear_mic_penalty = sum(linear_penalty_value(beta,sharpness,lambda0_T, ...
    mic_form,kink_coefficient)) + ...
    sum(linear_penalty_value(phi,sharpness,lambda0_C, ...
    mic_form,kink_coefficient));
linear_l1_penalty = l1_scale * (sum(l1_weights_beta.*abs(beta)) + ...
    sum(l1_weights_phi.*abs(phi)));
linear_ridge_penalty = 0.5 * ridge_scale * (sum(beta.^2) + sum(phi.^2));

support_T = spline_group_support(bspline_beta, num_basis, q);
support_C = spline_group_support(bspline_phi, num_basis, q);
group_penalty = iht_lambda_T*sum(support_T) + iht_lambda_C*sum(support_C);
baseline_penalty = 0.5*baseline_lambda_T*(seta*RT*seta') + ...
    0.5*baseline_lambda_C*(gama*RT*gama');

objective = loglik - linear_mic_penalty - linear_l1_penalty - ...
    linear_ridge_penalty - group_penalty - baseline_penalty;
end


function value = linear_penalty_value(coef,sharpness,lambda0,form,slope_multiplier)
form = lower(string(form));
abs_coef = abs(coef);
switch form
    case "smooth"
        value = lambda0*tanh(sharpness*coef.^2);
    case "kinked"
        value = lambda0.*tanh(sharpness*coef.^2 + ...
            slope_multiplier*abs_coef);
    case "scad"
        penalty_lambda = lambda0*slope_multiplier;
        scad_shape = 3.7;
        value = zeros(size(coef));
        first = abs_coef <= penalty_lambda;
        middle = abs_coef > penalty_lambda & ...
            abs_coef <= scad_shape*penalty_lambda;
        last = abs_coef > scad_shape*penalty_lambda;
        value(first) = penalty_lambda*abs_coef(first);
        value(middle) = (-abs_coef(middle).^2 + ...
            2*scad_shape*penalty_lambda*abs_coef(middle) - ...
            penalty_lambda^2)/(2*(scad_shape-1));
        value(last) = 0.5*(scad_shape+1)*penalty_lambda^2;
    case "mcp"
        penalty_lambda = lambda0*slope_multiplier;
        mcp_concavity = 3;
        value = 0.5*mcp_concavity*penalty_lambda^2*ones(size(coef));
        active = abs_coef <= mcp_concavity*penalty_lambda;
        value(active) = penalty_lambda*abs_coef(active) - ...
            abs_coef(active).^2/(2*mcp_concavity);
    otherwise
        error('Unsupported LinearMICForm: %s',form);
end
end


function delta = baseline_step_delta(family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix,step_seta,step_gama,step_size,base_loglik)
seta_trial = max(seta + step_size*step_seta, 1e-8);
gama_trial = max(gama + step_size*step_gama, 1e-8);
trial_loglik = mpl_loglik_alpha(family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi,seta_trial,gama_trial,Psi,psix);
delta = trial_loglik - base_loglik;
end


function [beta, phi, bspline_beta, bspline_phi, refit_diag] = ...
    post_selection_joint_coefficient_refit_MPL( ...
    copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
    seta, gama, Psi, psix, n_refit_iter, tk_beta, tk_phi, tk_spline_beta, ...
    tk_spline_phi, q)
% Unpenalized coefficient refit on the support selected before this call.
% Baseline coefficients and the copula parameter are held at their final MPL
% values. Backtracking guarantees that the observed likelihood cannot decrease.

num_basis = size(B, 2) / q;
linear_beta_selected = find(abs(beta) > 0);
linear_phi_selected = find(abs(phi) > 0);
spline_beta_selected = selected_spline_indices(bspline_beta, num_basis, q);
spline_phi_selected = selected_spline_indices(bspline_phi, num_basis, q);
refit_diag = struct('initial_loglik', NaN, 'final_loglik', NaN, ...
    'iterations', 0, 'rejected', false);
base_loglik = mpl_loglik_alpha(copula_family, alph, Z, B, Delta, ...
    beta, phi, bspline_beta, bspline_phi, seta, gama, Psi, psix);
refit_diag.initial_loglik = base_loglik;
if isempty(linear_beta_selected) && isempty(linear_phi_selected) && ...
        isempty(spline_beta_selected) && isempty(spline_phi_selected)
    refit_diag.final_loglik = base_loglik;
    return;
end

n = size(Z, 1);
step_state = 1;
for rr = 1:n_refit_iter
    eta_T = Z*beta + B*bspline_beta;
    eta_C = Z*phi + B*bspline_phi;
    exp_eta_T = exp(min(max(eta_T, -700), 700));
    exp_eta_C = exp(min(max(eta_C, -700), 700));
    H_T = sum(exp_eta_T.*(seta.*Psi), 2);
    H_C = sum(exp_eta_C.*(gama.*Psi), 2);
    S_T = exp(-H_T);
    S_C = exp(-H_C);
    [~,~,LAM1,LAM2,~,~] = copula_terms_generic( ...
        copula_family, S_T, S_C, alph, Delta);
    residual_T = Delta - Delta.*H_T - LAM1.*S_T.*H_T;
    residual_C = (1-Delta) - (1-Delta).*H_C - LAM2.*S_C.*H_C;
    grad_beta = -(sum(residual_T.*Z, 1))'/n;
    grad_phi = -(sum(residual_C.*Z, 1))'/n;
    grad_spline_beta = -(sum(residual_T.*B, 1))'/n;
    grad_spline_phi = -(sum(residual_C.*B, 1))'/n;
    selected_gradient = [grad_beta(linear_beta_selected); ...
        grad_phi(linear_phi_selected); ...
        grad_spline_beta(spline_beta_selected); ...
        grad_spline_phi(spline_phi_selected)];
    gradient_inf_norm = max(abs(selected_gradient), [], 'omitnan');
    if isempty(gradient_inf_norm) || gradient_inf_norm < 1e-6
        break;
    end

    accepted = false;
    start_step = min(1, 1.25*step_state);
    fallback_step = 1/max([tk_beta tk_phi tk_spline_beta tk_spline_phi]);
    direction_scales = unique([start_step fallback_step], 'stable');
    for direction_scale = direction_scales
        for backtrack = 0:20
            step = direction_scale*0.5^backtrack;
            beta_trial = beta;
            phi_trial = phi;
            spline_beta_trial = bspline_beta;
            spline_phi_trial = bspline_phi;
            beta_trial(linear_beta_selected) = beta(linear_beta_selected) - ...
                step*grad_beta(linear_beta_selected);
            phi_trial(linear_phi_selected) = phi(linear_phi_selected) - ...
                step*grad_phi(linear_phi_selected);
            spline_beta_trial(spline_beta_selected) = bspline_beta(spline_beta_selected) - ...
                step*grad_spline_beta(spline_beta_selected);
            spline_phi_trial(spline_phi_selected) = bspline_phi(spline_phi_selected) - ...
                step*grad_spline_phi(spline_phi_selected);
            trial_loglik = mpl_loglik_alpha(copula_family, alph, Z, B, Delta, ...
                beta_trial, phi_trial, spline_beta_trial, spline_phi_trial, ...
                seta, gama, Psi, psix);
            if isfinite(trial_loglik) && trial_loglik >= base_loglik - 1e-12
                accepted = true;
                break;
            end
        end
        if accepted
            break;
        end
    end
    if ~accepted
        refit_diag.rejected = true;
        break;
    end

    old_vector = [beta(linear_beta_selected); phi(linear_phi_selected); ...
        bspline_beta(spline_beta_selected); bspline_phi(spline_phi_selected)];
    new_vector = [beta_trial(linear_beta_selected); phi_trial(linear_phi_selected); ...
        spline_beta_trial(spline_beta_selected); spline_phi_trial(spline_phi_selected)];
    beta = beta_trial;
    phi = phi_trial;
    bspline_beta = spline_beta_trial;
    bspline_phi = spline_phi_trial;
    old_loglik = base_loglik;
    base_loglik = trial_loglik;
    step_state = step;
    refit_diag.iterations = rr;
    relative_change = norm(new_vector-old_vector)/(norm(old_vector)+1e-8);
    relative_loglik_change = abs(base_loglik-old_loglik)/(abs(old_loglik)+1);
    if relative_change < 1e-7 && relative_loglik_change < 1e-10
        break;
    end
end
refit_diag.final_loglik = base_loglik;
end


function [beta, phi, bspline_beta, bspline_phi, seta, gama, alph, refit_diag] = ...
    post_selection_joint_oracle_refit_MPL( ...
    copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
    seta, gama, Psi, psix, n_refit_iter, tk_beta, tk_phi, tk_spline_beta, ...
    tk_spline_phi, q, refit_estimate_alpha,step_options)
% Unpenalized observed-likelihood refit on a fixed selected support.
% Regression, spline, and both nonnegative baseline-hazard blocks are
% alternated with likelihood backtracking. The default keeps the working
% copula fixed; joint_oracle_alpha adds an accepted profile update.

num_basis = size(B, 2) / q;
if nargin < 21 || isempty(step_options)
    step_options = default_mpl_step_options(1);
end
step_strategy = lower(string(step_options.strategy));
linear_beta_selected = find(abs(beta) > 0);
linear_phi_selected = find(abs(phi) > 0);
spline_beta_selected = selected_spline_indices(bspline_beta, num_basis, q);
spline_phi_selected = selected_spline_indices(bspline_phi, num_basis, q);
n = size(Z, 1);
base_loglik = mpl_loglik_alpha(copula_family, alph, Z, B, Delta, ...
    beta, phi, bspline_beta, bspline_phi, seta, gama, Psi, psix);
refit_diag = struct('initial_loglik', base_loglik, 'final_loglik', base_loglik, ...
    'iterations', 0, 'rejected', false, 'converged', false, ...
    'coefficient_rejections', 0, 'baseline_rejections', 0, ...
    'alpha_updates', 0, 'alpha_rejections', 0, ...
    'coefficient_backtracks', 0, 'baseline_backtracks', 0, ...
    'final_coefficient_step', NaN, 'final_baseline_step', NaN, ...
    'min_loglik_delta', Inf, 'final_relative_change', NaN, ...
    'final_relative_loglik_change', NaN, 'final_tau', NaN, ...
    'final_alpha', alph, 'final_coefficient_score_inf', NaN, ...
    'final_baseline_kkt_inf', NaN, 'final_stationarity_inf', NaN, ...
    'stationarity_pass', false, 'stagnation_handoff', false, ...
    'fmincon_polish_attempted', false, ...
    'fmincon_polish_accepted', false, 'fmincon_exitflag', NaN, ...
    'fmincon_iterations', NaN, 'fmincon_firstorderopt', NaN, ...
    'gradient_check_max_abs', NaN, 'gradient_check_max_rel', NaN);
refit_diag.gradient_check_max_index = NaN;
refit_diag.gradient_check_analytic = NaN;
refit_diag.gradient_check_numeric = NaN;
refit_diag.gradient_check_block = "not_audited";

coefficient_step_state = 1;
baseline_step_state = 0.5;

for rr = 1:n_refit_iter
    old_vector = [beta(linear_beta_selected); phi(linear_phi_selected); ...
        bspline_beta(spline_beta_selected); bspline_phi(spline_phi_selected); ...
        seta(:); gama(:)];
    if refit_estimate_alpha
        old_vector = [old_vector; alph];
    end
    iteration_loglik = base_loglik;

    eta_T = Z*beta + B*bspline_beta;
    eta_C = Z*phi + B*bspline_phi;
    exp_eta_T = safe_exp(eta_T);
    exp_eta_C = safe_exp(eta_C);
    H_T_raw = sum(exp_eta_T.*(seta.*Psi), 2);
    H_C_raw = sum(exp_eta_C.*(gama.*Psi), 2);
    active_H_T = H_T_raw > 0 & H_T_raw < 700;
    active_H_C = H_C_raw > 0 & H_C_raw < 700;
    H_T = min(max(H_T_raw, 0), 700);
    H_C = min(max(H_C_raw, 0), 700);
    S_T = exp(-H_T);
    S_C = exp(-H_C);
    [~,~,LAM1,LAM2,~,~] = copula_terms_generic( ...
        copula_family, S_T, S_C, alph, Delta);
    residual_T = Delta-active_H_T.*( ...
        Delta.*H_T+LAM1.*S_T.*H_T);
    residual_C = (1-Delta)-active_H_C.*( ...
        (1-Delta).*H_C+LAM2.*S_C.*H_C);
    grad_beta = -(sum(residual_T.*Z, 1))'/n;
    grad_phi = -(sum(residual_C.*Z, 1))'/n;
    grad_spline_beta = -(sum(residual_T.*B, 1))'/n;
    grad_spline_phi = -(sum(residual_C.*B, 1))'/n;

    coefficient_accepted = isempty(linear_beta_selected) && ...
        isempty(linear_phi_selected) && isempty(spline_beta_selected) && ...
        isempty(spline_phi_selected);
    selected_gradient = [grad_beta(linear_beta_selected); ...
        grad_phi(linear_phi_selected); ...
        grad_spline_beta(spline_beta_selected); ...
        grad_spline_phi(spline_phi_selected)];
    coefficient_directional_slope = sum(selected_gradient.^2);
    if step_strategy == "auto"
        coefficient_start = min(1, ...
            coefficient_step_state*step_options.auto_expansion);
        coefficient_scales = [coefficient_start, ...
            1/max([tk_beta tk_phi tk_spline_beta tk_spline_phi])];
    else
        coefficient_scales = [1, ...
            1/max([tk_beta tk_phi tk_spline_beta tk_spline_phi])];
    end
    for direction_mode = 1:numel(coefficient_scales)
        if coefficient_accepted
            break;
        end
        direction_scale = coefficient_scales(direction_mode);
        for backtrack = 0:20
            step = direction_scale*0.5^backtrack;
            beta_trial = beta;
            phi_trial = phi;
            spline_beta_trial = bspline_beta;
            spline_phi_trial = bspline_phi;
            beta_trial(linear_beta_selected) = beta(linear_beta_selected) - ...
                step*grad_beta(linear_beta_selected);
            phi_trial(linear_phi_selected) = phi(linear_phi_selected) - ...
                step*grad_phi(linear_phi_selected);
            spline_beta_trial(spline_beta_selected) = bspline_beta(spline_beta_selected) - ...
                step*grad_spline_beta(spline_beta_selected);
            spline_phi_trial(spline_phi_selected) = bspline_phi(spline_phi_selected) - ...
                step*grad_spline_phi(spline_phi_selected);
            eta_T_trial = Z*beta_trial + B*spline_beta_trial;
            eta_C_trial = Z*phi_trial + B*spline_phi_trial;
            trial_loglik = mpl_loglik_from_eta(copula_family, alph, Delta, ...
                eta_T_trial, eta_C_trial, seta, gama, Psi, psix);
            sufficient_increase = 1e-4*step*coefficient_directional_slope;
            if isfinite(trial_loglik) && ...
                    trial_loglik >= base_loglik+sufficient_increase-1e-12
                coefficient_accepted = true;
                beta = beta_trial;
                phi = phi_trial;
                bspline_beta = spline_beta_trial;
                bspline_phi = spline_phi_trial;
                eta_T = eta_T_trial;
                eta_C = eta_C_trial;
                base_loglik = trial_loglik;
                coefficient_step_state = step;
                break;
            end
            refit_diag.coefficient_backtracks = ...
                refit_diag.coefficient_backtracks+1;
        end
    end
    if ~coefficient_accepted
        refit_diag.coefficient_rejections = refit_diag.coefficient_rejections + 1;
    end

    exp_eta_T = safe_exp(eta_T);
    exp_eta_C = safe_exp(eta_C);
    H_T_raw = sum(exp_eta_T.*(seta.*Psi), 2);
    H_C_raw = sum(exp_eta_C.*(gama.*Psi), 2);
    active_H_T = H_T_raw > 0 & H_T_raw < 700;
    active_H_C = H_C_raw > 0 & H_C_raw < 700;
    H_T = min(max(H_T_raw, 0), 700);
    H_C = min(max(H_C_raw, 0), 700);
    S_T = exp(-H_T);
    S_C = exp(-H_C);
    [~,~,LAM1,LAM2,LAM11,LAM22] = copula_terms_generic( ...
        copula_family, S_T, S_C, alph, Delta);
    score_seta = sum((Delta.*psix)./(sum(seta.*psix,2)+1e-8) - ...
        (active_H_T.*exp_eta_T.*(Delta+LAM1.*S_T)).*Psi, 1)/n;
    score_gama = sum(((1-Delta).*psix)./(sum(gama.*psix,2)+1e-8) - ...
        (active_H_C.*exp_eta_C.*(1-Delta+LAM2.*S_C)).*Psi, 1)/n;
    baseline_accepted = false;
    if abs(alph) < 1e-6
        % Under independence, each piecewise baseline block has a closed-form
        % profile MLE for fixed regression effects.
        seta_trial = max(sum(Delta.*psix,1) ./ ...
            (sum(exp_eta_T.*Psi,1)+1e-8), 1e-8);
        gama_trial = max(sum((1-Delta).*psix,1) ./ ...
            (sum(exp_eta_C.*Psi,1)+1e-8), 1e-8);
        trial_loglik = mpl_loglik_from_eta(copula_family, alph, Delta, ...
            eta_T, eta_C, seta_trial, gama_trial, Psi, psix);
        if isfinite(trial_loglik) && trial_loglik >= base_loglik - 1e-12
            baseline_accepted = true;
            seta = seta_trial;
            gama = gama_trial;
            base_loglik = trial_loglik;
            baseline_step_state = 1;
        end
    end
    if ~baseline_accepted
        % The score is for the mean log-likelihood, so its preconditioner
        % must use the same 1/n scale.
        Ps1 = sum((active_H_T.*Delta.*exp_eta_T).*Psi + ...
            ((active_H_T.*LAM11.*S_T).*exp_eta_T).*Psi, 1)/n + 1e-6;
        Ps2 = sum((active_H_C.*(1-Delta).*exp_eta_C).*Psi + ...
            ((active_H_C.*LAM22.*S_C).*exp_eta_C).*Psi, 1)/n + 1e-6;
        baseline_directions = {
            score_seta,              score_gama;
            score_seta.*(seta./Ps1), score_gama.*(gama./Ps2)
        };
        baseline_accepted = norm([score_seta score_gama]) <= 1e-8;
        if step_strategy == "auto"
            baseline_start = min(0.5,baseline_step_state*step_options.auto_expansion);
        else
            baseline_start = 0.5;
        end
        for direction_mode = 1:size(baseline_directions, 1)
            if baseline_accepted
                break;
            end
            step_seta = baseline_directions{direction_mode, 1};
            step_gama = baseline_directions{direction_mode, 2};
            for direction_sign = [1 -1]
                directional_slope = direction_sign*( ...
                    score_seta*step_seta'+score_gama*step_gama');
                if ~isfinite(directional_slope) || directional_slope <= 0
                    continue;
                end
                for backtrack = 0:29
                    step = baseline_start*0.5^backtrack;
                    seta_trial = max(seta + direction_sign*step*step_seta, 1e-8);
                    gama_trial = max(gama + direction_sign*step*step_gama, 1e-8);
                    trial_loglik = mpl_loglik_from_eta(copula_family, alph, Delta, ...
                        eta_T, eta_C, seta_trial, gama_trial, Psi, psix);
                    sufficient_increase = 1e-4*step*directional_slope;
                    if isfinite(trial_loglik) && ...
                            trial_loglik >= base_loglik+sufficient_increase-1e-12
                        baseline_accepted = true;
                        seta = seta_trial;
                        gama = gama_trial;
                        base_loglik = trial_loglik;
                        baseline_step_state = step;
                        break;
                    end
                    refit_diag.baseline_backtracks = ...
                        refit_diag.baseline_backtracks+1;
                end
                if baseline_accepted
                    break;
                end
            end
        end
    end
    if ~baseline_accepted
        refit_diag.baseline_rejections = refit_diag.baseline_rejections + 1;
    end

    if refit_estimate_alpha
        try
            [tau_trial, alph_trial] = update_copula_tau( ...
                copula_family, Z, B, Delta, beta, phi, bspline_beta, ...
                bspline_phi, seta, gama, Psi, psix);
            trial_loglik = mpl_loglik_from_eta(copula_family, alph_trial, Delta, ...
                eta_T, eta_C, seta, gama, Psi, psix);
            if isfinite(trial_loglik) && trial_loglik >= base_loglik - 1e-12
                alph = alph_trial;
                base_loglik = trial_loglik;
                refit_diag.alpha_updates = refit_diag.alpha_updates + 1;
                refit_diag.final_tau = tau_trial;
            else
                refit_diag.alpha_rejections = refit_diag.alpha_rejections + 1;
            end
        catch
            refit_diag.alpha_rejections = refit_diag.alpha_rejections + 1;
        end
    end

    new_vector = [beta(linear_beta_selected); phi(linear_phi_selected); ...
        bspline_beta(spline_beta_selected); bspline_phi(spline_phi_selected); ...
        seta(:); gama(:)];
    if refit_estimate_alpha
        new_vector = [new_vector; alph];
    end
    loglik_delta = base_loglik - iteration_loglik;
    refit_diag.min_loglik_delta = min(refit_diag.min_loglik_delta, loglik_delta);
    refit_diag.iterations = rr;
    if ~isfinite(base_loglik) || loglik_delta < -1e-10
        refit_diag.rejected = true;
        break;
    end
    relative_change = norm(new_vector-old_vector)^2/(norm(old_vector)^2+1e-8);
    relative_loglik_change = abs(loglik_delta)/(1+abs(iteration_loglik));
    refit_diag.final_relative_change = relative_change;
    refit_diag.final_relative_loglik_change = relative_loglik_change;
    coefficient_score_iteration = max(abs(selected_gradient),[],'omitnan');
    if isempty(coefficient_score_iteration)
        coefficient_score_iteration = 0;
    end
    baseline_value_iteration = [seta(:); gama(:)];
    baseline_score_iteration = [score_seta(:); score_gama(:)];
    positive_iteration = baseline_value_iteration > 1e-6;
    baseline_kkt_iteration = zeros(size(baseline_score_iteration));
    baseline_kkt_iteration(positive_iteration) = ...
        abs(baseline_score_iteration(positive_iteration));
    baseline_kkt_iteration(~positive_iteration) = ...
        max(baseline_score_iteration(~positive_iteration),0);
    stationarity_iteration = max(coefficient_score_iteration, ...
        max(baseline_kkt_iteration,[],'omitnan'));
    small_iteration_change = relative_change < 1e-6 && ...
        relative_loglik_change < 1e-5;
    if small_iteration_change
        if stationarity_iteration <= 1e-4
            refit_diag.converged = true;
            break;
        elseif step_options.fmincon_polish && rr >= 20
            % The alternating baseline update has stagnated. Hand the same
            % fixed-support likelihood to the constrained polish instead of
            % exhausting the outer iteration cap.
            refit_diag.stagnation_handoff = true;
            break;
        end
    end
end

if step_options.fmincon_polish && exist('fmincon','file') == 2
    refit_diag.fmincon_polish_attempted = true;
    [beta_trial,phi_trial,spline_beta_trial,spline_phi_trial, ...
        seta_trial,gama_trial,polish] = oracle_refit_fmincon_polish( ...
        copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
        seta,gama,Psi,psix,linear_beta_selected,linear_phi_selected, ...
        spline_beta_selected,spline_phi_selected, ...
        step_options.fmincon_algorithm, ...
        step_options.fmincon_scale_problem, ...
        step_options.fmincon_gradient_audit);
    refit_diag.fmincon_exitflag = polish.exitflag;
    refit_diag.fmincon_iterations = polish.iterations;
    refit_diag.fmincon_firstorderopt = polish.firstorderopt;
    refit_diag.gradient_check_max_abs = polish.gradient_check_max_abs;
    refit_diag.gradient_check_max_rel = polish.gradient_check_max_rel;
    refit_diag.gradient_check_max_index = polish.gradient_check_max_index;
    refit_diag.gradient_check_analytic = polish.gradient_check_analytic;
    refit_diag.gradient_check_numeric = polish.gradient_check_numeric;
    refit_diag.gradient_check_block = polish.gradient_check_block;
    if polish.accepted && polish.loglik >= base_loglik-1e-10
        beta = beta_trial;
        phi = phi_trial;
        bspline_beta = spline_beta_trial;
        bspline_phi = spline_phi_trial;
        seta = seta_trial;
        gama = gama_trial;
        base_loglik = polish.loglik;
        refit_diag.fmincon_polish_accepted = true;
    end
end
refit_diag.final_loglik = base_loglik;
refit_diag.final_alpha = alph;
refit_diag.final_coefficient_step = coefficient_step_state;
refit_diag.final_baseline_step = baseline_step_state;

% A small iterate or likelihood change can also occur after backtracking has
% reduced a step almost to zero. Record the fixed-support score/KKT residual
% so numerical convergence is distinguishable from stagnation.
eta_T = Z*beta + B*bspline_beta;
eta_C = Z*phi + B*bspline_phi;
exp_eta_T = safe_exp(eta_T);
exp_eta_C = safe_exp(eta_C);
H_T_raw = sum(exp_eta_T.*(seta.*Psi), 2);
H_C_raw = sum(exp_eta_C.*(gama.*Psi), 2);
active_H_T = H_T_raw > 0 & H_T_raw < 700;
active_H_C = H_C_raw > 0 & H_C_raw < 700;
H_T = min(max(H_T_raw, 0), 700);
H_C = min(max(H_C_raw, 0), 700);
S_T = exp(-H_T);
S_C = exp(-H_C);
[~,~,LAM1,LAM2,~,~] = copula_terms_generic( ...
    copula_family, S_T, S_C, alph, Delta);
residual_T = Delta-active_H_T.*(Delta.*H_T+LAM1.*S_T.*H_T);
residual_C = (1-Delta)-active_H_C.*( ...
    (1-Delta).*H_C+LAM2.*S_C.*H_C);
score_beta = sum(residual_T.*Z,1)'/n;
score_phi = sum(residual_C.*Z,1)'/n;
score_spline_beta = sum(residual_T.*B,1)'/n;
score_spline_phi = sum(residual_C.*B,1)'/n;
selected_score_raw = [score_beta(linear_beta_selected); ...
    score_phi(linear_phi_selected); ...
    score_spline_beta(spline_beta_selected); ...
    score_spline_phi(spline_phi_selected)];
% v2: spline scores are evaluated in G_j-orthonormal (estimable) coordinates,
% c_j = G_j^{1/2} b_j, so the 1e-4 rule does not depend on the B-spline
% parameterization. Linear scores are unchanged.
selected_score = [score_beta(linear_beta_selected); ...
    score_phi(linear_phi_selected); ...
    whitened_group_scores(B,score_spline_beta,spline_beta_selected,num_basis); ...
    whitened_group_scores(B,score_spline_phi,spline_phi_selected,num_basis)];
if isempty(selected_score)
    coefficient_score_inf = 0;
    coefficient_score_inf_raw = 0;
else
    coefficient_score_inf = max(abs(selected_score),[],'omitnan');
    coefficient_score_inf_raw = max(abs(selected_score_raw),[],'omitnan');
end
refit_diag.final_coefficient_score_inf_raw = coefficient_score_inf_raw;
score_seta = sum((Delta.*psix)./(sum(seta.*psix,2)+1e-8) - ...
    (active_H_T.*exp_eta_T.*(Delta+LAM1.*S_T)).*Psi,1)/n;
score_gama = sum(((1-Delta).*psix)./(sum(gama.*psix,2)+1e-8) - ...
    (active_H_C.*exp_eta_C.*(1-Delta+LAM2.*S_C)).*Psi,1)/n;
baseline_score = [score_seta(:); score_gama(:)];
baseline_value = [seta(:); gama(:)];
positive = baseline_value > 1e-6;
baseline_kkt = zeros(size(baseline_score));
baseline_kkt(positive) = abs(baseline_score(positive));
baseline_kkt(~positive) = max(baseline_score(~positive),0);
baseline_kkt_inf = max(baseline_kkt,[],'omitnan');
refit_diag.final_coefficient_score_inf = coefficient_score_inf;
refit_diag.final_baseline_kkt_inf = baseline_kkt_inf;
refit_diag.final_stationarity_inf = max(coefficient_score_inf,baseline_kkt_inf);
refit_diag.stationarity_pass = isfinite(refit_diag.final_stationarity_inf) && ...
    refit_diag.final_stationarity_inf <= 1e-4;
refit_diag.converged = ~refit_diag.rejected && refit_diag.stationarity_pass;
if ~isfinite(refit_diag.min_loglik_delta)
    refit_diag.min_loglik_delta = NaN;
end
end


function inference = fixed_support_model_inference_MPL( ...
        copula_family,alph,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
        seta,gama,Psi,psix,linear_beta_selected,linear_phi_selected, ...
        spline_beta_selected,spline_phi_selected)
% Model-based observed-information inference conditional on fixed support.
% The copula parameter is held fixed. Centered spline null directions and
% baseline hazards on the nonnegativity boundary are held fixed as well.

p = size(Z,2);
inference = empty_mpl_inference(p);
inference.attempted = true;
linear_beta_selected = linear_beta_selected(:);
linear_phi_selected = linear_phi_selected(:);
spline_beta_selected = spline_beta_selected(:);
spline_phi_selected = spline_phi_selected(:);

coefficient_start = [beta(linear_beta_selected);phi(linear_phi_selected); ...
    bspline_beta(spline_beta_selected);bspline_phi(spline_phi_selected)];
x0 = [coefficient_start;seta(:);gama(:)];
n_coefficient = numel(coefficient_start);

R_beta = selected_spline_reparameterization(B,spline_beta_selected);
R_phi = selected_spline_reparameterization(B,spline_phi_selected);
P_coefficient = blkdiag(eye(numel(linear_beta_selected)), ...
    eye(numel(linear_phi_selected)),R_beta,R_phi);
baseline_value = [seta(:);gama(:)];
baseline_active = find(baseline_value > 1e-6);
P = zeros(numel(x0),size(P_coefficient,2)+numel(baseline_active));
P(1:n_coefficient,1:size(P_coefficient,2)) = P_coefficient;
for j = 1:numel(baseline_active)
    P(n_coefficient+baseline_active(j),size(P_coefficient,2)+j) = 1;
end

dimension = size(P,2);
inference.dimension = dimension;
if dimension == 0
    inference.reason = "empty_interior_parameterization";
    return;
end

objective = @(x) oracle_refit_vector_objective(x,copula_family,alph,Z,B, ...
    Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix, ...
    linear_beta_selected,linear_phi_selected,spline_beta_selected, ...
    spline_phi_selected);
H = NaN(dimension,dimension);
for j = 1:dimension
    direction = P(:,j);
    coordinate = direction'*x0;
    step = eps(class(x0))^(1/3)*(1+abs(coordinate));
    baseline_row = find(abs(direction(n_coefficient+1:end)) > 0,1);
    if ~isempty(baseline_row)
        full_row = n_coefficient+baseline_row;
        step = min(step,0.25*(x0(full_row)-1e-8));
    end
    if ~isfinite(step) || step <= 0
        inference.reason = "invalid_finite_difference_step";
        return;
    end
    [~,gradient_plus] = objective(x0+step*direction);
    [~,gradient_minus] = objective(x0-step*direction);
    if any(~isfinite(gradient_plus)) || any(~isfinite(gradient_minus))
        inference.reason = "nonfinite_gradient";
        return;
    end
    H(:,j) = P'*(gradient_plus-gradient_minus)/(2*step);
end
H = (H+H')/2;
if any(~isfinite(H),'all')
    inference.reason = "nonfinite_observed_information";
    return;
end

eigenvalues = eig(H);
max_eigenvalue = max(eigenvalues);
min_eigenvalue = min(eigenvalues);
tolerance = max(1e-10,eps(class(H))*dimension*max(abs(eigenvalues)));
inference.rank = rank(H,tolerance);
inference.min_eigenvalue = min_eigenvalue;
if min_eigenvalue <= tolerance || inference.rank < dimension
    inference.reason = "observed_information_not_positive_definite";
    return;
end
inference.condition_number = max_eigenvalue/min_eigenvalue;
if ~isfinite(inference.condition_number) || inference.condition_number > 1e12
    inference.reason = "observed_information_ill_conditioned";
    return;
end

covariance_reduced = H\eye(dimension);
covariance_full = (P*covariance_reduced*P')/size(Z,1);
variance_full = diag(covariance_full);
if any(~isfinite(variance_full)) || any(variance_full < -1e-10)
    inference.reason = "invalid_covariance_diagonal";
    return;
end
standard_error = sqrt(max(variance_full,0));
cursor = 0;
n_beta = numel(linear_beta_selected);
inference.linear_beta_se(linear_beta_selected) = ...
    standard_error(cursor+(1:n_beta));
cursor = cursor+n_beta;
n_phi = numel(linear_phi_selected);
inference.linear_phi_se(linear_phi_selected) = ...
    standard_error(cursor+(1:n_phi));
inference.valid = true;
inference.reason = "ok";
end


function R = selected_spline_reparameterization(B,selected)
% Remove exact null directions induced by centered spline groups.
selected = selected(:);
if isempty(selected)
    R = zeros(0,0);
    return;
end
[~,S,V] = svd(B(:,selected),'econ');
singular_values = diag(S);
if isempty(singular_values)
    numerical_rank = 0;
else
    numerical_rank = sum(singular_values > ...
        max(size(B(:,selected)))*eps(max(singular_values)));
end
R = V(:,1:numerical_rank);
end


function [beta,phi,bspline_beta,bspline_phi,seta,gama,diagnostics] = ...
    oracle_refit_fmincon_polish(copula_family,alph,Z,B,Delta,beta,phi, ...
    bspline_beta,bspline_phi,seta,gama,Psi,psix,linear_beta_selected, ...
    linear_phi_selected,spline_beta_selected,spline_phi_selected, ...
    algorithm,scale_problem,gradient_audit)
% Constrained fixed-support polish for the same observed likelihood.

coefficient_start = [beta(linear_beta_selected);phi(linear_phi_selected); ...
    bspline_beta(spline_beta_selected);bspline_phi(spline_phi_selected)];
x0 = [coefficient_start;seta(:);gama(:)];
n_coefficient = numel(coefficient_start);
lower_bound = [-Inf(n_coefficient,1);1e-8*ones(numel(seta)+numel(gama),1)];
objective = @(x) oracle_refit_vector_objective(x,copula_family,alph,Z,B, ...
    Delta,beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix, ...
    linear_beta_selected,linear_phi_selected,spline_beta_selected, ...
    spline_phi_selected);
gradient_check_max_abs = NaN;
gradient_check_max_rel = NaN;
gradient_check_max_index = NaN;
gradient_check_analytic = NaN;
gradient_check_numeric = NaN;
gradient_check_block = "not_audited";
if gradient_audit
    [~,analytic_gradient] = objective(x0);
    numeric_gradient = finite_difference_gradient( ...
        objective,x0,lower_bound);
    difference = abs(analytic_gradient-numeric_gradient);
    [gradient_check_max_abs,gradient_check_max_index] = ...
        max(difference,[],'omitnan');
    gradient_check_max_rel = max(difference ./ ...
        max(1,abs(numeric_gradient)),[],'omitnan');
    gradient_check_analytic = analytic_gradient(gradient_check_max_index);
    gradient_check_numeric = numeric_gradient(gradient_check_max_index);
    block_ends = cumsum([numel(linear_beta_selected), ...
        numel(linear_phi_selected),numel(spline_beta_selected), ...
        numel(spline_phi_selected),numel(seta),numel(gama)]);
    block_names = ["linear_T","linear_C","additive_T","additive_C", ...
        "baseline_T","baseline_C"];
    gradient_check_block = block_names(find( ...
        gradient_check_max_index <= block_ends,1));
    audit_steps = 10.^(-3:-1:-7);
    fprintf("Gradient audit block=%s index=%d analytic=%.12g x=%.12g\n", ...
        gradient_check_block,gradient_check_max_index, ...
        gradient_check_analytic,x0(gradient_check_max_index));
    for audit_step = audit_steps
        plus = x0;
        minus = x0;
        plus(gradient_check_max_index) = ...
            plus(gradient_check_max_index)+audit_step;
        minus(gradient_check_max_index) = ...
            minus(gradient_check_max_index)-audit_step;
        numeric_at_step = (objective(plus)-objective(minus))/(2*audit_step);
        fprintf("  central step %.1e numeric %.12g\n", ...
            audit_step,numeric_at_step);
    end
    if gradient_check_block == "linear_T"
        coefficient_position = gradient_check_max_index;
        beta_index = linear_beta_selected(coefficient_position);
        component_step = 1e-5;
        plus = x0;
        minus = x0;
        plus(gradient_check_max_index) = ...
            plus(gradient_check_max_index)+component_step;
        minus(gradient_check_max_index) = ...
            minus(gradient_check_max_index)-component_step;
        [~,~,beta_plus,phi_plus,spline_beta_plus,spline_phi_plus, ...
            seta_plus,gama_plus] = objective(plus);
        [~,~,beta_minus,phi_minus,spline_beta_minus,spline_phi_minus, ...
            seta_minus,gama_minus] = objective(minus);
        [hazard_plus,copula_plus] = failure_linear_objective_components( ...
            copula_family,alph,Z,B,Delta,beta_plus,phi_plus, ...
            spline_beta_plus,spline_phi_plus,seta_plus,gama_plus,Psi);
        [hazard_minus,copula_minus] = failure_linear_objective_components( ...
            copula_family,alph,Z,B,Delta,beta_minus,phi_minus, ...
            spline_beta_minus,spline_phi_minus,seta_minus,gama_minus,Psi);
        numeric_hazard = -(hazard_plus-hazard_minus)/(2*component_step);
        numeric_copula = -(copula_plus-copula_minus)/(2*component_step);

        [~,~,beta_zero,phi_zero,spline_beta_zero,spline_phi_zero, ...
            seta_zero,gama_zero] = objective(x0);
        eta_T_zero = Z*beta_zero+B*spline_beta_zero;
        eta_C_zero = Z*phi_zero+B*spline_phi_zero;
        exp_eta_T_zero = safe_exp(eta_T_zero);
        exp_eta_C_zero = safe_exp(eta_C_zero);
        H_T_raw_zero = sum(exp_eta_T_zero.*(seta_zero.*Psi),2);
        H_C_raw_zero = sum(exp_eta_C_zero.*(gama_zero.*Psi),2);
        active_T_zero = H_T_raw_zero > 0 & H_T_raw_zero < 700;
        H_T_zero = min(max(H_T_raw_zero,0),700);
        H_C_zero = min(max(H_C_raw_zero,0),700);
        S_T_zero = exp(-H_T_zero);
        S_C_zero = exp(-H_C_zero);
        [~,~,LAM1_zero,~,~,~] = copula_terms_generic( ...
            copula_family,S_T_zero,S_C_zero,alph,Delta);
        z_column = Z(:,beta_index);
        analytic_hazard = -mean((Delta-active_T_zero.*Delta.*H_T_zero).* ...
            z_column);
        analytic_copula = mean(active_T_zero.*LAM1_zero.*S_T_zero.* ...
            H_T_zero.*z_column);
        fprintf(['  component linear_T variable=%d hazard analytic %.12g ' ...
            'numeric %.12g; copula analytic %.12g numeric %.12g\n'], ...
            beta_index,analytic_hazard,numeric_hazard, ...
            analytic_copula,numeric_copula);
        if strcmpi(copula_family,'clayton')
            eps0 = 1e-8;
            active_u_zero = S_T_zero > eps0 & S_T_zero < 1-eps0;
            u_zero = min(max(S_T_zero,eps0),1-eps0);
            v_zero = min(max(S_C_zero,eps0),1-eps0);
            A_zero = u_zero.^(-alph)+v_zero.^(-alph)-1;
            du_zero = u_zero.^(-alph-1)./A_zero;
            K1_zero = max(A_zero.^(-1/alph-1).* ...
                u_zero.^(-alph-1),realmin('double'));
            K2_zero = max(A_zero.^(-1/alph-1).* ...
                v_zero.^(-alph-1),realmin('double'));
            raw_dlogK1_du = (K1_zero./(K1_zero+eps0)).* ...
                (alph+1).*(du_zero-1./u_zero);
            raw_dlogK2_du = (K2_zero./(K2_zero+eps0)).* ...
                (alph+1).*du_zero;
            raw_LAM1 = active_u_zero.*(Delta.*raw_dlogK1_du+ ...
                (1-Delta).*raw_dlogK2_du);
            analytic_copula_raw = mean(active_T_zero.*raw_LAM1.* ...
                S_T_zero.*H_T_zero.*z_column);
            fprintf(['  Clayton raw copula analytic %.12g; clipped=%d/%d; ' ...
                'raw |LAM1| max %.12g\n'],analytic_copula_raw, ...
                nnz(abs(raw_LAM1)>1e6),numel(raw_LAM1), ...
                max(abs(raw_LAM1)));
        end
    end
end
options = optimoptions('fmincon','Algorithm',char(algorithm), ...
    'SpecifyObjectiveGradient',true,'Display','off','MaxIterations',1000, ...
    'MaxFunctionEvaluations',50000,'OptimalityTolerance',1e-7, ...
    'StepTolerance',1e-10,'ConstraintTolerance',1e-10, ...
    'ScaleProblem',logical(scale_problem));
[candidate,fval,exitflag,output] = fmincon(objective,x0,[],[],[],[], ...
    lower_bound,[],[],options);
[~,~,beta,phi,bspline_beta,bspline_phi,seta,gama] = ...
    oracle_refit_vector_objective(candidate,copula_family,alph,Z,B,Delta, ...
    beta,phi,bspline_beta,bspline_phi,seta,gama,Psi,psix, ...
    linear_beta_selected,linear_phi_selected,spline_beta_selected, ...
    spline_phi_selected);
diagnostics = struct('accepted',isfinite(fval), 'loglik',-fval, ...
    'exitflag',exitflag,'iterations',output.iterations, ...
    'firstorderopt',output.firstorderopt, ...
    'gradient_check_max_abs',gradient_check_max_abs, ...
    'gradient_check_max_rel',gradient_check_max_rel, ...
    'gradient_check_max_index',gradient_check_max_index, ...
    'gradient_check_analytic',gradient_check_analytic, ...
    'gradient_check_numeric',gradient_check_numeric, ...
    'gradient_check_block',gradient_check_block);
end


function [hazard_component,copula_component] = ...
        failure_linear_objective_components( ...
        family,theta,Z,B,Delta,beta,phi,bspline_beta,bspline_phi, ...
        seta,gama,Psi)
% Components that vary when a failure-margin linear coefficient changes.

eta_T = Z*beta+B*bspline_beta;
eta_C = Z*phi+B*bspline_phi;
H_T = min(max(sum(safe_exp(eta_T).*(seta.*Psi),2),0),700);
H_C = min(max(sum(safe_exp(eta_C).*(gama.*Psi),2),0),700);
S_T = exp(-H_T);
S_C = exp(-H_C);
[K1,K2,~,~,~,~] = copula_terms_generic( ...
    family,S_T,S_C,theta,Delta);
hazard_component = mean(Delta.*(eta_T-H_T));
copula_component = mean(Delta.*log(K1+1e-8)+ ...
    (1-Delta).*log(K2+1e-8));
end


function gradient = finite_difference_gradient(objective,x,lower_bound)
% Central differences except at an active lower bound.

gradient = NaN(size(x));
base_objective = objective(x);
for j = 1:numel(x)
    step = eps(class(x))^(1/3)*(1+abs(x(j)));
    plus = x;
    plus(j) = x(j)+step;
    objective_plus = objective(plus);
    if x(j)-step >= lower_bound(j)
        minus = x;
        minus(j) = x(j)-step;
        objective_minus = objective(minus);
        gradient(j) = (objective_plus-objective_minus)/(2*step);
    else
        gradient(j) = (objective_plus-base_objective)/step;
    end
end
end


function [objective,gradient,beta,phi,bspline_beta,bspline_phi,seta,gama] = ...
    oracle_refit_vector_objective(x,copula_family,alph,Z,B,Delta,beta,phi, ...
    bspline_beta,bspline_phi,seta,gama,Psi,psix,linear_beta_selected, ...
    linear_phi_selected,spline_beta_selected,spline_phi_selected)

cursor = 0;
n_block = numel(linear_beta_selected);
beta(linear_beta_selected) = x(cursor+(1:n_block));
cursor = cursor+n_block;
n_block = numel(linear_phi_selected);
phi(linear_phi_selected) = x(cursor+(1:n_block));
cursor = cursor+n_block;
n_block = numel(spline_beta_selected);
bspline_beta(spline_beta_selected) = x(cursor+(1:n_block));
cursor = cursor+n_block;
n_block = numel(spline_phi_selected);
bspline_phi(spline_phi_selected) = x(cursor+(1:n_block));
cursor = cursor+n_block;
seta = reshape(x(cursor+(1:numel(seta))),size(seta));
cursor = cursor+numel(seta);
gama = reshape(x(cursor+(1:numel(gama))),size(gama));

eta_T = Z*beta+B*bspline_beta;
eta_C = Z*phi+B*bspline_phi;
exp_eta_T = safe_exp(eta_T);
exp_eta_C = safe_exp(eta_C);
H_T_raw = sum(exp_eta_T.*(seta.*Psi),2);
H_C_raw = sum(exp_eta_C.*(gama.*Psi),2);
active_H_T = H_T_raw > 0 & H_T_raw < 700;
active_H_C = H_C_raw > 0 & H_C_raw < 700;
H_T = min(max(H_T_raw,0),700);
H_C = min(max(H_C_raw,0),700);
S_T = exp(-H_T);
S_C = exp(-H_C);
[~,~,LAM1,LAM2,~,~] = copula_terms_generic( ...
    copula_family,S_T,S_C,alph,Delta);
residual_T = Delta-active_H_T.*(Delta.*H_T+LAM1.*S_T.*H_T);
residual_C = (1-Delta)-active_H_C.*( ...
    (1-Delta).*H_C+LAM2.*S_C.*H_C);
n = size(Z,1);
gradient_beta = -(sum(residual_T.*Z,1))'/n;
gradient_phi = -(sum(residual_C.*Z,1))'/n;
gradient_spline_beta = -(sum(residual_T.*B,1))'/n;
gradient_spline_phi = -(sum(residual_C.*B,1))'/n;
score_seta = sum((Delta.*psix)./(sum(seta.*psix,2)+1e-8) - ...
    (active_H_T.*exp_eta_T.*(Delta+LAM1.*S_T)).*Psi,1)/n;
score_gama = sum(((1-Delta).*psix)./(sum(gama.*psix,2)+1e-8) - ...
    (active_H_C.*exp_eta_C.*(1-Delta+LAM2.*S_C)).*Psi,1)/n;
gradient = [gradient_beta(linear_beta_selected); ...
    gradient_phi(linear_phi_selected); ...
    gradient_spline_beta(spline_beta_selected); ...
    gradient_spline_phi(spline_phi_selected);-score_seta(:);-score_gama(:)];
objective = -mpl_loglik_from_eta(copula_family,alph,Delta,eta_T,eta_C, ...
    seta,gama,Psi,psix);
end


function [bspline_beta, bspline_phi] = post_selection_spline_refit_MPL( ...
    copula_family, alph, Z, B, Delta, beta, phi, bspline_beta, bspline_phi, ...
    seta, gama, Psi, n_refit_iter, tk_beta, tk_phi, q)
%% Post-selection refit for spline coefficients in the MPL estimator.
% The penalized algorithm first selects important nonparametric components.
% This step then refits the selected spline groups without additional group
% thresholding, which reduces shrinkage bias and improves fitted curves.

    [n, qs] = size(B);
    num_basis = qs / q;
    beta_selected = selected_spline_indices(bspline_beta, num_basis, q);
    phi_selected  = selected_spline_indices(bspline_phi,  num_basis, q);

    if isempty(beta_selected) && isempty(phi_selected)
        return;
    end

    for rr = 1:n_refit_iter
        eta_T = Z*beta + B*bspline_beta;
        eta_C = Z*phi  + B*bspline_phi;
        exp_eta_T = exp(min(max(eta_T,-700),700));
        exp_eta_C = exp(min(max(eta_C,-700),700));

        H_T = sum(exp_eta_T.*(seta.*Psi),2);
        H_C = sum(exp_eta_C.*(gama.*Psi),2);
        S_T = exp(-H_T);
        S_C = exp(-H_C);

        [~,~,LAM1,LAM2,~,~] = copula_terms_generic(copula_family,S_T,S_C,alph,Delta);

        if ~isempty(beta_selected)
            DL_beta_spline = -(sum((Delta - Delta.*H_T - LAM1.*S_T.*H_T).*B,1))'/n;
            old_beta = bspline_beta(beta_selected);
            bspline_beta(beta_selected) = bspline_beta(beta_selected) - DL_beta_spline(beta_selected)/tk_beta;
            beta_change = norm(bspline_beta(beta_selected)-old_beta)^2/(norm(old_beta)^2+1e-8);
        else
            beta_change = 0;
        end

        if ~isempty(phi_selected)
            DL_phi_spline = -(sum(((1-Delta) - (1-Delta).*H_C - LAM2.*S_C.*H_C).*B,1))'/n;
            old_phi = bspline_phi(phi_selected);
            bspline_phi(phi_selected) = bspline_phi(phi_selected) - DL_phi_spline(phi_selected)/tk_phi;
            phi_change = norm(bspline_phi(phi_selected)-old_phi)^2/(norm(old_phi)^2+1e-8);
        else
            phi_change = 0;
        end

        if max(beta_change, phi_change) < 1e-7
            break;
        end
    end
end


function [beta,bspline_beta,refit_diag] = post_selection_joint_refit_PL( ...
    Z,B,Delta,beta,bspline_beta,n_refit_iter,q)
% Unpenalized Cox partial-likelihood refit on the fixed PL support.

n = size(Z,1);
num_basis = size(B,2)/q;
linear_selected = find(abs(beta)>0);
spline_selected = selected_spline_indices(bspline_beta,num_basis,q);
linear_support = abs(beta)>0;
additive_support = false(q,1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    additive_support(j) = norm(bspline_beta(idx),2)>1e-10; % v2: same tolerance as spline_group_support
end
base_loglik = cox_partial_mean_loglik(Z,B,Delta,beta,bspline_beta);
refit_diag = struct('iterations',0,'converged',false,'rejections',0, ...
    'initial_loglik',base_loglik,'final_loglik',base_loglik, ...
    'final_relative_change',NaN,'final_relative_loglik_change',NaN, ...
    'linear_support',linear_support,'additive_support',additive_support);
if isempty(linear_selected) && isempty(spline_selected)
    refit_diag.converged = true;
    return;
end

for rr = 1:n_refit_iter
    old = [beta(linear_selected);bspline_beta(spline_selected)];
    iteration_loglik = base_loglik;
    eta = Z*beta+B*bspline_beta;
    risk = safe_exp(eta);
    denom = max(cumsum(risk),1e-8);
    grad_beta = -(sum(Delta.*(Z-cumsum(risk.*Z)./denom),1))'/n;
    grad_spline = -(sum(Delta.*(B-cumsum(risk.*B)./denom),1))'/n;
    accepted = false;
    for backtrack = 0:30
        step = 0.5^backtrack;
        beta_trial = beta;
        spline_trial = bspline_beta;
        beta_trial(linear_selected) = beta(linear_selected)- ...
            step*grad_beta(linear_selected);
        spline_trial(spline_selected) = bspline_beta(spline_selected)- ...
            step*grad_spline(spline_selected);
        trial_loglik = cox_partial_mean_loglik( ...
            Z,B,Delta,beta_trial,spline_trial);
        if isfinite(trial_loglik) && trial_loglik >= base_loglik-1e-12
            beta = beta_trial;
            bspline_beta = spline_trial;
            base_loglik = trial_loglik;
            accepted = true;
            break;
        end
    end
    refit_diag.iterations = rr;
    if ~accepted
        refit_diag.rejections = refit_diag.rejections+1;
        break;
    end
    current = [beta(linear_selected);bspline_beta(spline_selected)];
    relative_change = norm(current-old)^2/(norm(old)^2+1e-8);
    relative_loglik_change = abs(base_loglik-iteration_loglik)/ ...
        (1+abs(iteration_loglik));
    refit_diag.final_relative_change = relative_change;
    refit_diag.final_relative_loglik_change = relative_loglik_change;
    if relative_change < 1e-6 && relative_loglik_change < 1e-5
        refit_diag.converged = true;
        break;
    end
end
refit_diag.final_loglik = base_loglik;
end


function value = cox_partial_mean_loglik(Z,B,Delta,beta,bspline_beta)
eta = Z*beta+B*bspline_beta;
risk_sum = max(cumsum(safe_exp(eta)),realmin);
value = sum(Delta.*(eta-log(risk_sum)))/size(Z,1);
end


function bspline_beta = post_selection_spline_refit_PL(Z, B, Delta, beta, bspline_beta, n_refit_iter, tk_spline, q)
%% Post-selection refit for spline coefficients in the PL estimator.
% Only the selected spline groups are updated, and no additional thresholding
% is applied during this refit step.

    [n, qs] = size(B);
    num_basis = qs / q;
    selected = selected_spline_indices(bspline_beta, num_basis, q);

    if isempty(selected)
        return;
    end

    for rr = 1:n_refit_iter
        eta = Z*beta + B*bspline_beta;
        exp_eta = exp(min(max(eta,-700),700));
        denom = cumsum(exp_eta) + 1e-8;
        weighted_B = cumsum(exp_eta.*B)./denom;
        DL_spline = -(sum(Delta.*(B-weighted_B),1))'/n;

        old = bspline_beta(selected);
        bspline_beta(selected) = bspline_beta(selected) - DL_spline(selected)/tk_spline;
        rel_change = norm(bspline_beta(selected)-old)^2/(norm(old)^2+1e-8);
        if rel_change < 1e-7
            break;
        end
    end
end


function selected = selected_spline_indices(bspline_coef, num_basis, q)
%% Return coefficient indices for nonzero selected spline groups.

    selected = [];
    for jj = 1:q
        idx = ((jj-1)*num_basis+1):(jj*num_basis);
        if norm(bspline_coef(idx),2) > 1e-10
            selected = [selected, idx]; %#ok<AGROW>
        end
    end
    selected = selected(:);
end


%% ===================================================
%                 PL()
% ============================================================
function [opt_beta,opt_bspline_beta,refit_diag] = PL(n,q,ini_beta,ini_bspline_beta,Z,X,B,T,r,Delta,a,block_size,q_size,threshold_linear,threshold_spline,mic_penalty_scale,mic_penalty_normalization,spline_threshold_stage,spline_selection_metric,spline_group_penalty_scale,spline_group_penalty_type,spline_group_penalty_normalization,spline_group_mic_sharpness,spline_group_mic_continuation,spline_group_mic_sharpness_path,spline_group_mic_form,spline_group_kink_coefficient,spline_group_mcp_concavity,spline_group_iht_threshold,spline_group_iht_warmup_iterations,linear_mic_form,linear_kink_coefficient,post_selection_refit_max_iter)
if nargin < 14 || isempty(threshold_linear)
    threshold_linear = 5e-4;
end
if nargin < 16 || isempty(mic_penalty_scale)
    mic_penalty_scale = 1;
end
if nargin < 17 || isempty(mic_penalty_normalization)
    mic_penalty_normalization = "sum";
end
if nargin < 18 || isempty(spline_threshold_stage)
    spline_threshold_stage = "iterative";
end
spline_threshold_stage = lower(string(spline_threshold_stage));
if nargin < 19 || isempty(spline_selection_metric)
    spline_selection_metric = "coef";
end
spline_selection_metric = lower(string(spline_selection_metric));
if nargin < 20 || isempty(spline_group_penalty_scale)
    spline_group_penalty_scale = 15;
end
if nargin < 21 || isempty(spline_group_penalty_type)
    spline_group_penalty_type = "group_mic";
end
spline_group_penalty_type = lower(string(spline_group_penalty_type));
if nargin < 22 || isempty(spline_group_penalty_normalization)
    spline_group_penalty_normalization = "sum";
end
spline_group_penalty_normalization = lower(string(spline_group_penalty_normalization));
if nargin < 23 || isempty(spline_group_mic_sharpness)
    spline_group_mic_sharpness = 90;
end
if nargin < 24 || isempty(spline_group_mic_continuation)
    spline_group_mic_continuation = true;
end
if nargin < 25
    spline_group_mic_sharpness_path = [];
end
if nargin < 26 || isempty(spline_group_mic_form)
    spline_group_mic_form = "smooth";
end
if nargin < 27 || isempty(spline_group_kink_coefficient)
    spline_group_kink_coefficient = 0;
end
if nargin < 28 || isempty(spline_group_mcp_concavity)
    spline_group_mcp_concavity = 3;
end
if nargin < 29 || isempty(spline_group_iht_threshold)
    spline_group_iht_threshold = 0.1;
end
if nargin < 30 || isempty(spline_group_iht_warmup_iterations)
    spline_group_iht_warmup_iterations = 10;
end
if nargin < 31 || isempty(linear_mic_form)
    linear_mic_form = "smooth";
end
if nargin < 32 || isempty(linear_kink_coefficient)
    linear_kink_coefficient = 0;
end
if nargin < 33 || isempty(post_selection_refit_max_iter)
    post_selection_refit_max_iter = 30;
end
spline_group_mic_form = lower(string(spline_group_mic_form));
linear_mic_form = lower(string(linear_mic_form));
if spline_group_penalty_type == "group_iht" && threshold_spline ~= 0
    error('PLAPH:GroupIHTDoubleThreshold', ...
        'ThresholdSpline must equal zero when group-IHT is active.');
end
if spline_threshold_stage == "final"
    threshold_spline_iter = 0;
else
    threshold_spline_iter = threshold_spline;
end
[~,p] = size(Z);
[~,qs] = size(B);

beta = ini_beta;
beta1 = zeros(p,1);
bspline_beta = ini_bspline_beta;

lambda0 = mic_penalty_lambda(sum(Delta), n, mic_penalty_scale, mic_penalty_normalization);  % MIC/SELO, Su et al. (2016)
spline_lambda_T = spline_group_penalty_lambda(sum(Delta), n, spline_group_penalty_scale, spline_group_penalty_normalization);


%% Set up block coordinate descent parameters
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% Construct block indices
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


%% Set up block coordinate descent parameters--spline
num_basis = qs/q;
spline_group_metric = [];
if spline_group_penalty_type == "group_iht"
    spline_group_metric = prepare_spline_group_metric(B, num_basis, q);
end
block_size_spline = block_size*num_basis;  % num_basis  % number of spline basis functions
num_blocks_spline = ceil(q / block_size);
blocks_spline = cell(num_blocks_spline, 1);

%% Construct block indices
for b = 1:num_blocks_spline
    start_idx1 = (b-1)*block_size_spline + 1;
    end_idx1 = min(b*block_size_spline, q*num_basis);
    blocks_spline{b} = start_idx1:end_idx1;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1;

tk1 = 200; tk2 = 200; tk3 = 100; tk4 = 100;

% beta = ini_beta; phi = ini_phi;
while k<=1000
    %k
    spline_penalty_type_iter = spline_group_penalty_type;
    spline_lambda_T_iter = spline_lambda_T;
    if spline_group_penalty_type == "group_iht" && k <= spline_group_iht_warmup_iterations
        spline_penalty_type_iter = "group_mic";
    end

    %% 
    risk_beta = safe_exp(Z*beta+B*bspline_beta);
    denom_beta = max(cumsum(risk_beta), 1e-8);
    DL_beta = -(sum(Delta.*(Z-cumsum(risk_beta.*Z)./denom_beta)))'/n;
    DL_beta(~isfinite(DL_beta)) = 0;

 %% Select the update strategy according to the dimension
 if p<=30
     [mic_weight,l1_coefficient] = linear_penalty_local_terms( ...
         beta,a,lambda0,linear_mic_form,linear_kink_coefficient,0,1);
     beta_tilde = beta - DL_beta/tk1;
     beta1 = hybrid_linear_prox(beta_tilde,mic_weight,tk1, ...
         l1_coefficient,tk1*1e-6);
 else
     block_order1 = randperm(num_blocks);
     for b1 = block_order1
         current_block1 = blocks{b1};
         beta_block = beta(current_block1);
         [mic_weight,l1_coefficient] = linear_penalty_local_terms( ...
             beta_block,a,lambda0,linear_mic_form, ...
             linear_kink_coefficient,0,1);
         beta_tilde = beta_block - DL_beta(current_block1)/tk1;
         beta_tilde(~isfinite(beta_tilde)) = 0;
         beta_new1 = hybrid_linear_prox(beta_tilde,mic_weight,tk1, ...
             l1_coefficient,0);
         beta_new1(~isfinite(beta_new1)) = 0;
         beta1(current_block1,1) = beta_new1;
     end
 end
        
    w1 = beta1-beta;
    beta = beta1;

    
    %%
    risk_beta = safe_exp(Z*beta+B*bspline_beta);
    denom_beta = max(cumsum(risk_beta), 1e-8);
    DL_beta_spline = -(sum(Delta.*(B-cumsum(risk_beta.*B)./denom_beta)))'/n;
    DL_beta_spline(~isfinite(DL_beta_spline)) = 0;
    
    % bspline_beta1 = bspline_beta - DL_beta_spline/(tk1/2);
    if spline_penalty_type_iter == "group_iht"
        bspline_beta1 = apply_spline_group_iht_metric( ...
            spline_group_metric, bspline_beta, DL_beta_spline, qs/q, q, tk2, spline_group_iht_threshold);
    elseif q<=q_size
        bspline_beta1 = bspline_beta - DL_beta_spline/tk2;
    else
        for iter = 1:50
            % Randomly permute the block order
            block_order_spline = randperm(num_blocks_spline);
            %%
            for b = block_order_spline
                current_block_spline = blocks_spline{b};
                bspline_beta_block = bspline_beta(current_block_spline);

                bspline_beta_new1 =  bspline_beta_block - DL_beta_spline(current_block_spline)/tk2;
                bspline_beta_new1(~isfinite(bspline_beta_new1)) = 0;

                % Update parameters
                bspline_beta1(current_block_spline,1) = bspline_beta_new1;
            end

            % Check convergence
            delta = norm(bspline_beta - bspline_beta1,2)^2 / (norm(bspline_beta,2)^2 + 1e-8);
            if delta < 1e-6
                break;
            end
        end
    end

    if spline_penalty_type_iter ~= "group_iht"
        bspline_beta1 = apply_spline_group_penalty_step( ...
            bspline_beta1, qs/q, q, spline_lambda_T_iter, tk2, spline_penalty_type_iter, ...
            spline_group_mic_sharpness, spline_group_mic_continuation, spline_group_mic_sharpness_path, ...
            spline_group_mic_form, spline_group_kink_coefficient, spline_group_mcp_concavity, ...
            spline_group_iht_threshold);
    end
    bspline_beta1 = apply_spline_group_threshold_metric( ...
        B, bspline_beta1, qs/q, q, threshold_spline_iter, spline_selection_metric);

    w3 = bspline_beta1-bspline_beta;
    bspline_beta = bspline_beta1;

    k = k+1;

end

%%
beta2 = beta.*(abs(beta)>=threshold_linear);
% Post-selection spline refit for the PL estimator.
if spline_threshold_stage == "final"
    bspline_beta = apply_spline_group_threshold_metric( ...
        B, bspline_beta, qs/q, q, threshold_spline, spline_selection_metric);
end
[beta2,bspline_beta,refit_diag] = post_selection_joint_refit_PL( ...
    Z,B,Delta,beta2,bspline_beta,post_selection_refit_max_iter,q);

opt_beta = beta2;
opt_bspline_beta = bspline_beta;

% H_T = sum(opt_seta.*Psi,2);
% H_C = sum(opt_gama.*Psi,2);

end




%% ===================================================
%                 oracle_value()
% ============================================================
function [opt_beta,opt_spline_beta]= oracle_value(~,Z,X,B,~,status)
% Failure-margin Cox benchmark fitted on the known active support.
% Small nonzero markers declare the fixed support to the common refit routine;
% they are starting values only and are replaced by unpenalized updates.
p_oracle = size(Z,2);
q_oracle = size(X,2);
beta_start = 1e-6*ones(p_oracle,1);
spline_start = 1e-6*ones(size(B,2),1);
[opt_beta,opt_spline_beta] = post_selection_joint_refit_PL( ...
    Z,B,status,beta_start,spline_start,1000,q_oracle);


end

%% ===========================================================
%                 safe_exp()
% ============================================================
function y = safe_exp(x)
% Numerically stable exponential function.
% This prevents overflow and underflow in exp(x).

    y = exp(min(max(x, -700), 700));

end


%% ===========================================================
%                 v2 helper functions
% ============================================================
function order = plaph_spline_order(value)
% Persistent B-spline order used by every basis constructor (4 = cubic,
% 3 = quadratic). Set once at the start of each scenario.
persistent stored_order
if nargin > 0
    if ~(isscalar(value) && any(value == [2 3 4 5]))
        error('PLAPH:SplineOrder','SplineOrder must be 2, 3, 4, or 5.');
    end
    stored_order = value;
end
if isempty(stored_order)
    stored_order = 4;
end
order = stored_order;
end


function opts = plaph_v2_default_options(opts)
% Defaults for options introduced in v2, so older drivers remain usable.
defaults = struct( ...
    'TuningMode', "fixed", ...          % "fixed" | "grid_ebic" | "grid_bic"
    'TuningGridCL', [0.5 1.0 1.5 2.0 2.5], ...
    'TuningGridKappa', [0.03 0.05 0.075 0.10 0.15], ...
    'PLVersion', "v2", ...               % "v2" aligned comparator | "legacy"
    'ExportDataDir', "", ...             % per-replication data for grpsurv
    'LinearPenaltyTuningCriterion', "none", ...
    'LinearPenaltyMultiplierGrid', 25.6685, ...
    'LinearPenaltyEBICGamma', 1, ...
    'LinearL1Type', "none", ...
    'LinearL1Scale', 0);
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(opts,names{j}) || isempty(opts.(names{j}))
        opts.(names{j}) = defaults.(names{j});
    end
end
end


function [tk1,tk2,tk3,tk4,count] = plaph_double_curvatures(tk1,tk2,tk3,tk4,count)
tk1 = 2*tk1; tk2 = 2*tk2; tk3 = 2*tk3; tk4 = 2*tk4;
count = count + 1;
end


function w = trapezoid_weights(m)
w = ones(m,1);
if m > 1
    w([1 m]) = 0.5;
end
w = w/sum(w);
end


function g = selection_gate_v2(copula_family,alph,Z,B,Delta,beta,phi, ...
    bT,bC,seta,gama,Psi,psix,k,q,metric,LT,LC,thT,thC,lambdaT,lambdaC,a,xi, ...
    mask_lin_T,mask_lin_C,mask_spl_T,mask_spl_C,RT,lambda1,lambda2)
if nargin < 29
    RT = []; lambda1 = 0; lambda2 = 0;
end
% Necessary first-order conditions of the penalized selection objective,
% evaluated with the same derivative expressions as the update maps:
% kinked-MIC scalar KKT, nonnegative-baseline KKT, and the group-IHT
% fixed-point residual. Valid for every supported copula family.
n = size(Z,1);
Delta = double(Delta(:));
eta_T = Z*beta + B*bT;
eta_C = Z*phi + B*bC;
exp_T = safe_exp(eta_T);
exp_C = safe_exp(eta_C);
H_T = min(max(sum(exp_T.*(seta.*Psi),2),0),700);
H_C = min(max(sum(exp_C.*(gama.*Psi),2),0),700);
S_T = exp(-H_T);
S_C = exp(-H_C);
[~,~,LAM1,LAM2,~,~] = copula_terms_generic(copula_family,S_T,S_C,alph,Delta);
res_T = Delta - Delta.*H_T - LAM1.*S_T.*H_T;
res_C = (1-Delta) - (1-Delta).*H_C - LAM2.*S_C.*H_C;
g.linear_gradient_T = -(Z'*res_T)/n;
g.linear_gradient_C = -(Z'*res_C)/n;
grad_bT = -(B'*res_T)/n;
grad_bC = -(B'*res_C)/n;
g.linear_gradient_T(~mask_lin_T) = 0;
g.linear_gradient_C(~mask_lin_C) = 0;
grad_bT(~mask_spl_T) = 0;
grad_bC(~mask_spl_C) = 0;
DL_seta = sum((Delta.*psix)./(sum(seta.*psix,2)+1e-8) - ...
    (exp_T.*(Delta+LAM1.*S_T)).*Psi,1)/n;
DL_gama = sum(((1-Delta).*psix)./(sum(gama.*psix,2)+1e-8) - ...
    (exp_C.*((1-Delta)+LAM2.*S_C)).*Psi,1)/n;
if ~isempty(RT)
    % Same baseline roughness term as the update map and the objective.
    DL_seta = DL_seta - lambda1*seta*RT;
    DL_gama = DL_gama - lambda2*gama*RT;
end
g.linear_T = gate_linear_kkt(beta(mask_lin_T),g.linear_gradient_T(mask_lin_T), ...
    subset_penalty(lambdaT,find(mask_lin_T)),a,xi);
g.linear_C = gate_linear_kkt(phi(mask_lin_C),g.linear_gradient_C(mask_lin_C), ...
    subset_penalty(lambdaC,find(mask_lin_C)),a,xi);
g.baseline_T = gate_baseline_kkt(seta(:),DL_seta(:));
g.baseline_C = gate_baseline_kkt(gama(:),DL_gama(:));
g.group_T = gate_group_residual(metric,B,bT,grad_bT,k,q,LT,thT,n);
g.group_C = gate_group_residual(metric,B,bC,grad_bC,k,q,LC,thC,n);
g.maximum = max([g.linear_T,g.linear_C,g.baseline_T,g.baseline_C,g.group_T,g.group_C]);
g.pass = isfinite(g.maximum) && g.maximum <= 1e-4;
g.status = 'evaluated';
end


function x = gate_linear_kkt(b,grad,lambda,a,xi)
if isempty(b)
    x = 0;
    return;
end
inner = a*b.^2 + xi*abs(b);
if isscalar(lambda)
    lambda = lambda*ones(size(b));
end
pgrad = lambda.*(1-tanh(inner).^2).*(2*a*b + xi*sign(b));
r = abs(grad + pgrad);
zero = b == 0;
r(zero) = max(abs(grad(zero)) - lambda(zero)*xi, 0);
x = max(r);
end


function x = gate_baseline_kkt(b,grad)
% grad is the gradient of the mean log-likelihood (to be maximized).
r = abs(grad);
boundary = b <= 1e-6;
r(boundary) = max(grad(boundary),0);
x = max(r);
end


function x = gate_group_residual(metric,B,b,grad,k,q,L,threshold,n)
if q == 0
    x = 0;
    return;
end
mapped = apply_spline_group_iht_metric(metric,b,grad,k,q,L,threshold);
r = zeros(q,1);
current = zeros(q,1);
for j = 1:q
    idx = (j-1)*k + (1:k);
    r(j) = norm(B(:,idx)*(mapped(idx)-b(idx)))/sqrt(n);
    current(j) = norm(B(:,idx)*b(idx))/sqrt(n);
end
x = L*norm(r)/max(1,norm(current));
end


function w = whitened_group_scores(B,score,selected,num_basis)
% Scores in G_j-orthonormal coordinates for the selected spline columns:
% for c_j = G_j^{1/2} b_j the score is G_j^{+1/2} times the b-score.
w = zeros(0,1);
if isempty(selected)
    return;
end
selected = selected(:);
groups = unique(ceil(selected/num_basis));
n = size(B,1);
for j = groups(:)'
    idx = ((j-1)*num_basis+1):(j*num_basis);
    gram = (B(:,idx)'*B(:,idx))/n;
    gram = (gram+gram')/2;
    [V,D] = eig(gram);
    d = real(diag(D));
    tol = max(1e-10, max([d;0])*num_basis*eps);
    inv_sqrt = zeros(size(d));
    keep = d > tol;
    inv_sqrt(keep) = 1./sqrt(d(keep));
    w = [w; V*(inv_sqrt.*(V'*score(idx)))]; %#ok<AGROW>
end
end


function r = joint_selected_rank(Z,B,num_basis,linear_mask,group_mask)
cols_B = spline_group_columns(find(group_mask(:)),num_basis);
X = [Z(:,logical(linear_mask(:))), B(:,cols_B)];
if isempty(X)
    r = 0;
    return;
end
s = svd(X,'econ');
r = sum(s > max(size(X))*eps(max(s)));
end


%% ===========================================================
%                 PL_v2(): aligned partial-likelihood comparator
% ============================================================
function [opt_beta,opt_bspline_beta,diag] = PL_v2(n,q,ini_beta,ini_bspline_beta, ...
    Z,B,Delta,a,threshold_linear,lambda_scale,kink_coefficient,kappa, ...
    step_scale,pilot_options,max_outer_iter,refit_max_iter)
% Cox partial-likelihood counterpart of the PLAPH selection stage. It uses
% the same adaptive weighted lasso/group-lasso pilot, the same K-MIC and
% group hard-thresholding maps, the same base curvatures (200 and 100 times
% the step scale), the same step-doubling safeguard, the same first-order
% stopping gate and the same refit acceptance rule. Rows of Z and B must be
% sorted by decreasing observed time so that cumulative sums give risk sets.
Delta = double(Delta(:));
p = size(Z,2);
qs = size(B,2);
num_basis = 0;
if q > 0
    num_basis = qs/q;
end
n_events = sum(Delta);
lambda0 = mic_penalty_lambda(n_events,n,lambda_scale,"mean");
t_lin_base = 200*step_scale;
t_spl_base = 100*step_scale;
t_lin = t_lin_base;
t_spl = t_spl_base;
t_lin_acc = t_lin;   % same target-stage acceleration as MPL (see MPL comments)
t_spl_acc = t_spl;
lambda_A = 0.5*t_spl_base*kappa^2;
metric = [];
if q > 0
    metric = prepare_spline_group_metric(B,num_basis,q);
end
pilot_enabled = pilot_options.iterations > 0;
w_lin = ones(p,1);
w_grp = ones(q,1);
if pilot_enabled && lower(string(pilot_options.weighting)) == "adaptive"
    w_lin = min((abs(ini_beta(:))+pilot_options.weight_epsilon).^ ...
        (-pilot_options.weight_gamma),pilot_options.weight_cap);
    if q > 0
        w_grp = min((spline_effect_norms(B,ini_bspline_beta,num_basis,q)+ ...
            pilot_options.weight_epsilon).^(-pilot_options.group_weight_gamma), ...
            pilot_options.weight_cap);
    end
end
pilot_rank = max(num_basis-1,1);
pilot_lambda_lin = pilot_options.linear_scale*sqrt(max(n_events,1)*log(max(p,2)))/n;
pilot_lambda_grp = pilot_options.group_scale* ...
    sqrt(max(n_events,1)*(pilot_rank+log(max(q,2))))/n;
beta = ini_beta(:);
bspl = ini_bspline_beta(:);
if pilot_enabled && lower(string(pilot_options.start)) == "zero"
    beta(:) = 0;
    bspl(:) = 0;
end
max_total = max_outer_iter + pilot_enabled*pilot_options.iterations;
doublings = 0;
max_doublings = 40;
gate = struct('pass',false,'maximum',NaN);
k = 1;
status = 4;
while k <= max_total
    in_pilot = pilot_enabled && k <= pilot_options.iterations;
    beta_old = beta;
    bspl_old = bspl;
    if in_pilot
        obj_old = pl_pilot_objective(Z,B,Delta,beta,bspl,pilot_lambda_lin*w_lin, ...
            pilot_lambda_grp*w_grp,metric,num_basis,q);
    else
        obj_old = pl_target_objective(Z,B,Delta,beta,bspl,a,lambda0, ...
            kink_coefficient,lambda_A,num_basis,q);
    end
    [g_beta,~] = pl_gradients(Z,B,Delta,beta,bspl);
    if in_pilot
        z = beta - g_beta/t_lin;
        thr = pilot_lambda_lin*w_lin/t_lin;
        beta = sign(z).*max(abs(z)-thr,0);
    else
        [mic_weight,l1_coefficient] = linear_penalty_local_terms( ...
            beta,a,lambda0,"kinked",kink_coefficient,0,1);
        beta = hybrid_linear_prox(beta-g_beta/t_lin_acc,mic_weight,t_lin_acc, ...
            l1_coefficient,0);
    end
    beta(~isfinite(beta)) = 0;
    if q > 0
        [~,g_spl] = pl_gradients(Z,B,Delta,beta,bspl);
        if in_pilot
            bspl = apply_spline_group_lasso_metric(metric,bspl,g_spl, ...
                num_basis,q,t_spl,pilot_lambda_grp*w_grp);
        else
            bspl = apply_spline_group_iht_metric_acc(metric,bspl,g_spl, ...
                num_basis,q,t_spl,kappa*sqrt(t_spl_base/t_spl),t_spl_acc);
        end
    end
    if in_pilot
        obj_new = pl_pilot_objective(Z,B,Delta,beta,bspl,pilot_lambda_lin*w_lin, ...
            pilot_lambda_grp*w_grp,metric,num_basis,q);
    else
        obj_new = pl_target_objective(Z,B,Delta,beta,bspl,a,lambda0, ...
            kink_coefficient,lambda_A,num_basis,q);
    end
    if ~isfinite(obj_new) || obj_new > obj_old + 1e-8
        beta = beta_old;
        bspl = bspl_old;
        if ~in_pilot && (t_lin_acc < t_lin || t_spl_acc < t_spl)
            t_lin_acc = min(2*t_lin_acc,t_lin);
            t_spl_acc = min(2*t_spl_acc,t_spl);
            continue;
        end
        t_lin = 2*t_lin;
        t_spl = 2*t_spl;
        t_lin_acc = t_lin;
        t_spl_acc = t_spl;
        doublings = doublings + 1;
        if doublings > max_doublings
            status = 3;
            break;
        end
        continue;
    end
    if ~in_pilot
        t_lin_acc = max(t_lin_acc/1.5,t_lin/100);
        t_spl_acc = max(t_spl_acc/1.5,t_spl/100);
        small = norm(beta-beta_old)^2 <= 1e-5*(norm(beta_old)^2+1e-8) && ...
            norm(bspl-bspl_old)^2 <= 1e-5*(norm(bspl_old)^2+1e-8);
        if small || mod(k,25) == 0 || k == max_total
            gate = pl_selection_gate(Z,B,Delta,beta,bspl,a,lambda0, ...
                kink_coefficient,metric,num_basis,q,t_spl, ...
                kappa*sqrt(t_spl_base/t_spl));
            if small && gate.pass
                status = 1;
                break;
            end
        end
    end
    k = k + 1;
end
if status ~= 3
    gate = pl_selection_gate(Z,B,Delta,beta,bspl,a,lambda0,kink_coefficient, ...
        metric,num_basis,q,t_spl,kappa*sqrt(t_spl_base/t_spl));
    if gate.pass
        status = 1;
    elseif k > max_total
        status = 2;
    end
end
beta2 = beta.*(abs(beta) >= threshold_linear);
linear_support = abs(beta2) > 0;
additive_support = false(q,1);
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    additive_support(j) = norm(bspl(idx),2) > 1e-10;
end
bspl(~repelem(additive_support,num_basis)) = 0;
[beta2,bspl,refit] = pl_newton_refit(Z,B,Delta,beta2,bspl,linear_support, ...
    additive_support,num_basis,refit_max_iter);
opt_beta = beta2;
opt_bspline_beta = bspl;
diag = refit;
diag.linear_support = linear_support;
diag.additive_support = additive_support;
diag.selection_gate_pass = logical(gate.pass) && status == 1;
diag.selection_gate_maximum = gate.maximum;
diag.selection_status_code = status;
diag.step_doublings = doublings;
diag.outer_iterations = k;
diag.lambda0 = lambda0;
diag.lambda_A = lambda_A;
diag.accepted = diag.selection_gate_pass && refit.stationarity_pass;
end


function [g_beta,g_spl] = pl_gradients(Z,B,Delta,beta,bspl)
% Gradient of the NEGATIVE mean Cox partial log-likelihood (Breslow, no ties;
% rows sorted by decreasing time so the risk set of row i is rows 1..i).
n = size(Z,1);
eta = Z*beta + B*bspl;
risk = safe_exp(eta);
denom = max(cumsum(risk),1e-300);
g_beta = -(sum(Delta.*(Z - cumsum(risk.*Z)./denom),1))'/n;
if nargout > 1
    if isempty(B)
        g_spl = zeros(0,1);
    else
        g_spl = -(sum(Delta.*(B - cumsum(risk.*B)./denom),1))'/n;
    end
end
g_beta(~isfinite(g_beta)) = 0;
end


function value = pl_pilot_objective(Z,B,Delta,beta,bspl,lin_pen,grp_pen,metric,num_basis,q)
value = -cox_partial_mean_loglik(Z,B,Delta,beta,bspl) + sum(lin_pen(:).*abs(beta(:)));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    value = value + grp_pen(j)*sqrt(max(real(bspl(idx)'*metric.gram{j}*bspl(idx)),0));
end
end


function value = pl_target_objective(Z,B,Delta,beta,bspl,a,lambda0,xi,lambda_A,num_basis,q)
value = -cox_partial_mean_loglik(Z,B,Delta,beta,bspl) + ...
    sum(linear_penalty_value(beta,a,lambda0,"kinked",xi));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    value = value + lambda_A*double(norm(bspl(idx)) > 1e-10);
end
end


function g = pl_selection_gate(Z,B,Delta,beta,bspl,a,lambda0,xi,metric,num_basis,q,L,threshold)
n = size(Z,1);
[g_beta,g_spl] = pl_gradients(Z,B,Delta,beta,bspl);
g.linear = gate_linear_kkt(beta,g_beta,lambda0,a,xi);
g.group = gate_group_residual(metric,B,bspl,g_spl,num_basis,q,L,threshold,n);
g.maximum = max([g.linear g.group]);
g.pass = isfinite(g.maximum) && g.maximum <= 1e-4;
end


function [beta,bspl,diag] = pl_newton_refit(Z,B,Delta,beta,bspl,linear_support, ...
    additive_support,num_basis,max_iter)
% Unpenalized Cox partial-likelihood refit on a fixed support by damped
% Newton-Raphson in an estimable reparameterization of the spline blocks.
% Acceptance: maximum absolute score (linear, and spline in G_j-orthonormal
% coordinates) at most 1e-4, the same rule as for PLAPH.
n = size(Z,1);
lin = find(linear_support(:));
spl = spline_group_columns(find(additive_support(:)),num_basis);
R = selected_spline_reparameterization(B,spl);
X = [Z(:,lin), B(:,spl)*R];
theta = [beta(lin); R'*bspl(spl)];
diag = struct('iterations',0,'converged',false,'rejections',0, ...
    'final_loglik',NaN,'stationarity_inf',NaN,'stationarity_pass',false);
if isempty(X)
    diag.final_loglik = cox_partial_mean_loglik(Z,B,Delta,beta,bspl);
    diag.stationarity_inf = 0;
    diag.stationarity_pass = true;
    diag.converged = true;
    return;
end
d = size(X,2);
loglik = @(th) cox_partial_mean_loglik(X,zeros(n,0),Delta,th,zeros(0,1));
current = loglik(theta);
for it = 1:max_iter
    eta = X*theta;
    risk = safe_exp(eta);
    S0 = cumsum(risk);
    S1 = cumsum(risk.*X);
    xbar = S1./S0;
    score = (sum(Delta.*(X - xbar),1))'/n;
    info = zeros(d,d);
    S2 = zeros(d,d);
    for i = 1:n
        S2 = S2 + risk(i)*(X(i,:)'*X(i,:));
        if Delta(i) > 0
            info = info + (S2/S0(i) - xbar(i,:)'*xbar(i,:));
        end
    end
    info = (info+info')/(2*n);
    step = (info + 1e-10*eye(d))\score;
    accepted = false;
    for bt = 0:30
        trial = theta + 0.5^bt*step;
        value = loglik(trial);
        if isfinite(value) && value >= current - 1e-12
            accepted = true;
            break;
        end
    end
    diag.iterations = it;
    if ~accepted
        diag.rejections = diag.rejections + 1;
        break;
    end
    theta = trial;
    improvement = value - current;
    current = value;
    if max(abs(score)) < 1e-9 || abs(improvement) < 1e-14
        break;
    end
end
beta(:) = 0;
beta(lin) = theta(1:numel(lin));
bspl(:) = 0;
if ~isempty(spl)
    bspl(spl) = R*theta(numel(lin)+1:end);
end
[g_beta,g_spl] = pl_gradients(Z,B,Delta,beta,bspl);
scores = [-g_beta(lin); -whitened_group_scores(B,g_spl,spl,num_basis)];
diag.final_loglik = cox_partial_mean_loglik(Z,B,Delta,beta,bspl);
diag.stationarity_inf = max([abs(scores);0]);
diag.stationarity_pass = isfinite(diag.stationarity_inf) && diag.stationarity_inf <= 1e-4;
diag.converged = diag.stationarity_pass && diag.rejections == 0;
end


function [beta,bspl,diag] = run_pl_v2_tuned(n,q,ini_beta,ini_bspl,Z,B,Delta,a, ...
    threshold_linear,xi,cL_fixed,kappa_fixed,step_scale,pilot_options, ...
    max_outer_iter,refit_max_iter,opts,num_basis,ebic_gamma)
% Fit the aligned PL comparator at the fixed tuning values or, when
% opts.TuningMode is a grid mode, over the same (c_L, kappa) grid as PLAPH
% with the analogous partial-likelihood EBIC (cold start at every point).
tuning_mode = lower(string(opts.TuningMode));
if tuning_mode == "fixed"
    grid = [cL_fixed kappa_fixed];
else
    [gc,gk] = ndgrid(opts.TuningGridCL(:),opts.TuningGridKappa(:));
    grid = [gc(:) gk(:)];
end
best = Inf; best_any = Inf;
fits = cell(size(grid,1),1);
crit = Inf(size(grid,1),1);
valid = false(size(grid,1),1);
for g = 1:size(grid,1)
    [b,s,d] = PL_v2(n,q,ini_beta,ini_bspl,Z,B,Delta,a,threshold_linear, ...
        grid(g,1),xi,grid(g,2),step_scale,pilot_options,max_outer_iter,refit_max_iter);
    df = joint_selected_rank(Z,B,num_basis,d.linear_support,d.additive_support);
    model_space = log_binomial(size(Z,2),sum(d.linear_support)) + ...
        log_binomial(q,sum(d.additive_support));
    value = -2*n*d.final_loglik + log(n)*df;
    if tuning_mode == "grid_ebic"
        value = value + 2*ebic_gamma*model_space;
    end
    crit(g) = value;
    valid(g) = d.accepted && isfinite(value);
    d.tuning_cL = grid(g,1);
    d.tuning_kappa = grid(g,2);
    d.criterion = value;
    fits{g} = struct('beta',b,'bspl',s,'diag',d);
end
candidates = find(valid);
if isempty(candidates)
    candidates = find(isfinite(crit));
end
if isempty(candidates)
    error('PLAPH:PLNoFiniteCandidate','No finite PL candidate.');
end
if tuning_mode == "fixed"
    [~,pos] = min(crit(candidates));
    chosen = candidates(pos);
else
    eligible = false(size(crit));
    eligible(candidates) = true;
    chosen = plaph_tuning_tie_break(crit,eligible,grid(:,1),grid(:,2), ...
        opts.TuningGridCL,opts.TuningGridKappa);
end
beta = fits{chosen}.beta;
bspl = fits{chosen}.bspl;
diag = fits{chosen}.diag;
diag.iterations = diag.iterations;
diag.rejections = diag.rejections;
end


function v = subset_penalty(lambda,idx)
if isscalar(lambda)
    v = lambda;
else
    v = lambda(idx);
end
end


function [C,names] = clinical_adjustment_design(clinical)
% Age (years), male indicator, and AJCC stage II/III/IV indicators (stage I
% reference). Missing age, sex, or stage gives NaN rows (dropped upstream).
n = height(clinical);
needed = ["age" "sex" "stage"];
missing = needed(~ismember(needed,string(clinical.Properties.VariableNames)));
if ~isempty(missing)
    error('PLAPH:ClinicalColumns','Clinical file lacks %s; available: %s', ...
        strjoin(missing,", "),strjoin(string(clinical.Properties.VariableNames),", "));
end
age = double(clinical.age);
sex = lower(strtrim(string(clinical.sex)));
male = NaN(n,1);
male(sex == "male") = 1;
male(sex == "female") = 0;
stage = upper(strtrim(string(clinical.stage)));
s2 = NaN(n,1); s3 = NaN(n,1); s4 = NaN(n,1);
is4 = startsWith(stage,"STAGE IV");
is3 = startsWith(stage,"STAGE III");
is2 = startsWith(stage,"STAGE II") & ~is3;
is1 = startsWith(stage,"STAGE I") & ~is2 & ~is3 & ~is4;
known = is1 | is2 | is3 | is4;
s2(known) = double(is2(known));
s3(known) = double(is3(known));
s4(known) = double(is4(known));
C = [age male s2 s3 s4];
names = ["age";"male";"stage_II";"stage_III";"stage_IV"];
end


function coef_new = apply_spline_group_iht_metric_acc(metric, coef, gradient, ...
    num_basis, q, step_decide, threshold, step_update)
% Group IHT with the retention decision made at curvature step_decide (the
% definition of the estimator) and the update of retained groups made with
% curvature step_update <= step_decide (acceleration). With step_update equal
% to step_decide this is apply_spline_group_iht_metric.
coef_new = zeros(size(coef));
for j = 1:q
    idx = ((j-1)*num_basis+1):(j*num_basis);
    gram = metric.gram{j};
    direction = metric.gram_pinv{j}*gradient(idx);
    candidate = coef(idx) - direction/step_decide;
    candidate(~isfinite(candidate)) = 0;
    effect_norm = sqrt(max(real(candidate' * gram * candidate), 0));
    if effect_norm > threshold
        updated = coef(idx) - direction/step_update;
        updated(~isfinite(updated)) = 0;
        coef_new(idx) = updated;
    end
end
end


function [idx,tie_count,on_boundary] = plaph_tuning_tie_break(value,eligible, ...
    cL,kappa,grid_cL,grid_kappa)
% Grid choice with a deterministic tie rule. Among eligible candidates whose
% criterion is within max(1e-4, 1e-6|min|) of the minimum, choose the largest
% c_L and then the largest kappa. Distinct selected models differ in EBIC by
% at least log(n) - O(1), so the tolerance only absorbs refit rounding.
% on_boundary is true when no tied candidate is interior to the grid.
value = value(:);
eligible = logical(eligible(:)) & isfinite(value);
if ~any(eligible)
    eligible = isfinite(value);
end
m = min(value(eligible));
tol = max(1e-4,1e-6*abs(m));
tie = find(eligible & value <= m + tol);
[~,order] = sortrows([cL(tie) kappa(tie)],[-1 -2]);
idx = tie(order(1));
tie_count = numel(tie);
interior = true(numel(tie),1);
if numel(grid_cL) > 1
    interior = interior & cL(tie) > min(grid_cL) & cL(tie) < max(grid_cL);
end
if numel(grid_kappa) > 1
    interior = interior & kappa(tie) > min(grid_kappa) & kappa(tie) < max(grid_kappa);
end
on_boundary = ~any(interior);
end
