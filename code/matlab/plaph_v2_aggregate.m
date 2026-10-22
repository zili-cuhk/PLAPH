function S = plaph_v2_aggregate(out_root, varargin)
%PLAPH_V2_AGGREGATE Summarize the v2 campaign into CSV and LaTeX fragments.
%
%   S = plaph_v2_aggregate("/path/to/out")
%   S = plaph_v2_aggregate("/path/to/out", "Experiments", ["main" "ablation_tau0"])
%
% Reads every chunk's scenario_001_replicates.csv (plain text, no HDF5),
% the grpsurv CSVs under <out_root>/grpsurv_results/<scenario_id>/*.csv,
% and writes <out_root>/summary/*.csv and *.tex. Metrics follow the
% manuscript definitions; acceptance rates are always reported, and the
% main metrics are computed over all replications (with accepted-only
% versions in separate columns).

parser = inputParser;
parser.addParameter("Experiments", ["main" "ablation_tau0" "highcorr" ...
    "copula_family" "xi_sensitivity" "spline_pilot" "smoke"]);
parser.parse(varargin{:});
o = parser.Results;
out_root = char(out_root);
summary_dir = fullfile(out_root, "summary");
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end
S = struct();
for experiment = string(o.Experiments)
    exp_dir = fullfile(out_root, char(experiment));
    if ~isfolder(exp_dir)
        continue;
    end
    M = plaph_v2_manifest(experiment);
    rows = cell(height(M),1);
    for i = 1:height(M)
        R = read_scenario(fullfile(exp_dir, char(M.scenario_id(i))));
        if isempty(R)
            continue;
        end
        rows{i} = summarize_scenario(R, M(i,:), experiment, out_root);
    end
    rows = rows(~cellfun(@isempty,rows));
    if isempty(rows)
        continue;
    end
    T = vertcat(rows{:});
    writetable(T, fullfile(summary_dir, "summary_" + experiment + ".csv"));
    G = tuning_optimal_share(exp_dir, M);
    if ~isempty(G)
        writetable(G, fullfile(summary_dir, "tuning_optimal_share_" + experiment + ".csv"));
    end
    S.(char(experiment)) = T;
    if experiment == "main"
        write_main_tables(T, summary_dir);
        write_inference_table(exp_dir, M, summary_dir);
    end
end
if isfield(S,"main") && isfield(S,"ablation_tau0")
    A = paired_ablation(out_root);
    writetable(A, fullfile(summary_dir, "ablation_paired.csv"));
    S.ablation_paired = A;
    write_ablation_table(A, summary_dir);
end
end


%% ----------------------------------------------------------------------
function G = tuning_optimal_share(exp_dir, M)
% For each grid point (c_L, kappa): share of replications in which it is
% EBIC-optimal among accepted candidates, ties within max(1e-4,1e-6|min|)
% included (same rule as the engine). Pooled over scenarios and per scenario.
G = table();
parts = {};
for i = 1:height(M)
    files = dir(fullfile(exp_dir, char(M.scenario_id(i)), "chunk_*", ...
        "scenario_*_tuning_path.csv"));
    for f = 1:numel(files)
        P = readtable(fullfile(files(f).folder, files(f).name));
        if ~all(ismember(["c_L" "kappa" "criterion_value" "numerically_valid" "seed"], ...
                string(P.Properties.VariableNames)))
            continue;
        end
        P.scenario_index = repmat(M.scenario_index(i), height(P), 1);
        parts{end+1} = P(:, ["scenario_index" "seed" "c_L" "kappa" ...
            "criterion_value" "numerically_valid"]); %#ok<AGROW>
    end
end
if isempty(parts)
    return;
end
P = vertcat(parts{:});
P.optimal = zeros(height(P),1);
[grp,~] = findgroups(P.scenario_index, P.seed);
for g = 1:max(grp)
    k = find(grp == g);
    v = P.criterion_value(k);
    e = logical(P.numerically_valid(k)) & isfinite(v);
    if ~any(e)
        continue;
    end
    m = min(v(e));
    P.optimal(k) = double(e & v <= m + max(1e-4, 1e-6*abs(m)));
end
pooled = groupsummary(P, ["c_L" "kappa"], "mean", "optimal");
pooled.scenario_index = zeros(height(pooled),1);
per = groupsummary(P, ["scenario_index" "c_L" "kappa"], "mean", "optimal");
G = [pooled(:, ["scenario_index" "c_L" "kappa" "GroupCount" "mean_optimal"]); ...
     per(:, ["scenario_index" "c_L" "kappa" "GroupCount" "mean_optimal"])];
G.Properties.VariableNames = ["scenario_index" "c_L" "kappa" "n_reps" "optimal_share"];
end


function R = read_scenario(scenario_dir)
R = [];
if ~isfolder(scenario_dir)
    return;
end
files = dir(fullfile(scenario_dir, "chunk_*", "scenario_*_replicates.csv"));
parts = cell(numel(files),1);
for f = 1:numel(files)
    parts{f} = numeric_table(readtable(fullfile(files(f).folder, files(f).name), ...
        "TextType","string"));
end
parts = parts(~cellfun(@isempty,parts));
if isempty(parts)
    return;
end
common = parts{1}.Properties.VariableNames;
for f = 2:numel(parts)
    common = intersect(common, parts{f}.Properties.VariableNames, "stable");
end
for f = 1:numel(parts)
    parts{f} = parts{f}(:, common);
end
R = vertcat(parts{:});
R = sortrows(R, "seed");
end


function T = numeric_table(T)
% Keep numeric/logical columns only, converting text columns that hold numbers.
keep = false(1,width(T));
for j = 1:width(T)
    v = T.(j);
    if isnumeric(v) || islogical(v)
        T.(j) = double(v);
        keep(j) = true;
    elseif isstring(v) || iscellstr(v)
        d = str2double(string(v));
        lower_v = lower(strtrim(string(v)));
        d(lower_v == "true") = 1;
        d(lower_v == "false") = 0;
        if all(isfinite(d) | lower_v == "" | lower_v == "nan")
            T.(j) = d;
            keep(j) = true;
        end
    end
end
T = T(:, keep);
end


function x = col(R, name)
if ismember(name, R.Properties.VariableNames)
    x = R.(name);
else
    x = NaN(height(R),1);
end
end


function m = mean_finite(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end


function T = summarize_scenario(R, row, experiment, out_root)
beta_true = [1.5 -1.2 0.9];
gamma_true = [-1.2 1.5 0.9];
exact = @(fp,fn) double(fp == 0 & fn == 0);
acc = col(R,"mpl_accepted") == 1;
T = table();
T.experiment = string(experiment);
T.scenario_index = row.scenario_index;
T.scenario_id = row.scenario_id;
T.n = row.n; T.p = row.p; T.censor_target = row.target_censor_rate;
T.tau = row.tau; T.rho = row.rho; T.fit_tau = row.fit_tau; T.xi = row.xi;
T.generate_copula = row.generate_copula; T.fit_copula = row.fit_copula;
T.R_done = height(R);
T.ok_rate = mean_finite(col(R,"ok"));
T.censor_rate = mean_finite(col(R,"censor_rate"));
T.screen_all_sure = mean_finite(col(R,"screen_all_sure"));
T.accept_rate = mean(acc);
T.selection_gate_rate = mean_finite(col(R,"mpl_selection_gate_pass"));
T.refit_pass_rate = mean_finite(col(R,"mpl_refit_stationarity_pass"));
T.max_refit_residual = max(col(R,"mpl_refit_stationarity_inf"),[],"omitnan");
T.median_target_iterations = median(col(R,"mpl_target_iterations"),"omitnan");
T.step_doubling_rate = mean_finite(double(col(R,"mpl_step_doublings") > 0));
T.pilot_final_agree_add_T = mean_finite(col(R,"mpl_pilot_final_agree_add_T"));
T.pilot_final_agree_add_C = mean_finite(col(R,"mpl_pilot_final_agree_add_C"));
T.pilot_final_agree_lin_T = mean_finite(col(R,"mpl_pilot_final_agree_linear_T"));
T.pilot_final_agree_lin_C = mean_finite(col(R,"mpl_pilot_final_agree_linear_C"));
T.tuning_cL_mode = mode(col(R,"tuning_cL"));
T.tuning_kappa_mode = mode(col(R,"tuning_kappa"));
T.tuning_boundary_rate = mean_finite(col(R,"tuning_at_boundary"));
% Share of replications whose selected model is EBIC-optimal at every
% accepted grid point (the tuning is then immaterial for that replication).
T.tuning_full_tie_rate = mean_finite(double(col(R,"tuning_tie_count") == ...
    col(R,"tuning_valid_count")));
% Support recovery (PLAPH)
T.exact_beta = mean_finite(exact(col(R,"beta_fp"),col(R,"beta_fn")));
T.exact_gamma = mean_finite(exact(col(R,"phi_fp"),col(R,"phi_fn")));
T.exact_phiT = mean_finite(exact(col(R,"mpl_add_T_fp"),col(R,"mpl_add_T_fn")));
T.exact_phiC = mean_finite(exact(col(R,"mpl_add_C_fp"),col(R,"mpl_add_C_fn")));
T.fp_beta = mean_finite(col(R,"beta_fp")); T.fn_beta = mean_finite(col(R,"beta_fn"));
T.fp_gamma = mean_finite(col(R,"phi_fp")); T.fn_gamma = mean_finite(col(R,"phi_fn"));
T.fp_phiT = mean_finite(col(R,"mpl_add_T_fp")); T.fn_phiT = mean_finite(col(R,"mpl_add_T_fn"));
T.fp_phiC = mean_finite(col(R,"mpl_add_C_fp")); T.fn_phiC = mean_finite(col(R,"mpl_add_C_fn"));
T.exact_all = mean_finite(exact(col(R,"beta_fp"),col(R,"beta_fn")) .* ...
    exact(col(R,"phi_fp"),col(R,"phi_fn")) .* ...
    exact(col(R,"mpl_add_T_fp"),col(R,"mpl_add_T_fn")) .* ...
    exact(col(R,"mpl_add_C_fp"),col(R,"mpl_add_C_fn")));
% Estimation (PLAPH)
T.mse_beta = coef_mse(R,"mpl_beta",beta_true);
T.mse_gamma = coef_mse(R,"mpl_phi",gamma_true);
T.imse_phiT = mean_finite(col(R,"imse_mpl_add_T"));
T.imse_phiC = mean_finite(col(R,"imse_mpl_add_C"));
T.imse_phiT_accepted = mean_finite(col(R(acc,:),"imse_mpl_add_T"));
T.imse_phiC_accepted = mean_finite(col(R(acc,:),"imse_mpl_add_C"));
% PL comparator
T.pl_accept_rate = mean_finite(col(R,"pl_accepted"));
T.pl_exact_beta = mean_finite(exact(col(R,"pl_fp"),col(R,"pl_fn")));
T.pl_exact_phiT = mean_finite(exact(col(R,"pl_add_fp"),col(R,"pl_add_fn")));
T.pl_fp_phiT = mean_finite(col(R,"pl_add_fp")); T.pl_fn_phiT = mean_finite(col(R,"pl_add_fn"));
T.pl_mse_beta = coef_mse(R,"pl_beta",beta_true);
T.pl_imse_phiT = mean_finite(col(R,"imse_pl_add"));
d = col(R,"imse_mpl_add_T") - col(R,"imse_pl_add");
T.imse_diff_PL_minus_PLAPH = -mean_finite(d);
T.imse_diff_mcse = std(d(isfinite(d)))/sqrt(sum(isfinite(d)));
T.pr_PLAPH_lt_PL = mean_finite(double(d(isfinite(d)) < 0));
% grpsurv comparator (failure margin)
G = read_grpsurv(out_root, row.scenario_id);
T.grp_R = height(G);
T.grp_exact_beta = NaN; T.grp_exact_phiT = NaN; T.grp_mse_beta = NaN; T.grp_imse_phiT = NaN;
T.grp_fp_phiT = NaN; T.grp_fn_phiT = NaN;
if ~isempty(G)
    T.grp_exact_beta = mean(G.lin_fp == 0 & G.lin_fn == 0);
    T.grp_exact_phiT = mean(G.add_fp == 0 & G.add_fn == 0);
    T.grp_fp_phiT = mean(G.add_fp); T.grp_fn_phiT = mean(G.add_fn);
    T.grp_mse_beta = mean(G.mse_beta_active);
    T.grp_imse_phiT = mean(G.imse_add);
end
end


function m = coef_mse(R, prefix, truth)
v = NaN(1,numel(truth));
for j = 1:numel(truth)
    x = col(R, sprintf("%s_active%d",prefix,j));
    v(j) = mean_finite((x - truth(j)).^2);
end
m = mean(v,"omitnan");
end


function G = read_grpsurv(out_root, scenario_id)
G = [];
grp_scenario = regexprep(string(scenario_id),"^V2_[a-z0-9_]+?_S","V2_grpsurv_export_S");
d = fullfile(out_root, "grpsurv_results", char(grp_scenario));
files = dir(fullfile(d, "rows_*.csv"));
files = files(~endsWith({files.name}, "_curves.csv"));   % fitted curves, not result rows
if isempty(files)
    return;
end
parts = cell(numel(files),1);
for f = 1:numel(files)
    parts{f} = readtable(fullfile(files(f).folder, files(f).name));
end
G = vertcat(parts{:});
end


function write_inference_table(exp_dir, M, summary_dir)
% Coefficient-level bias, ESE, ASE, CP: unconditional and conditional on
% exact recovery of the corresponding linear support.
beta_true = [1.5 -1.2 0.9]; gamma_true = [-1.2 1.5 0.9];
rows = {};
for i = 1:height(M)
    R = read_scenario(fullfile(exp_dir, char(M.scenario_id(i))));
    if isempty(R)
        continue;
    end
    for margin = ["beta" "phi"]
        if margin == "beta"
            truth = beta_true; ok = col(R,"beta_fp") == 0 & col(R,"beta_fn") == 0;
        else
            truth = gamma_true; ok = col(R,"phi_fp") == 0 & col(R,"phi_fn") == 0;
        end
        for j = 1:3
            est = col(R, sprintf("mpl_%s_active%d",margin,j));
            ase = col(R, sprintf("mpl_%s_ase_active%d",margin,j));
            cp = col(R, sprintf("mpl_%s_cp_active%d",margin,j));
            for conditional = [false true]
                sel = true(height(R),1);
                if conditional
                    sel = ok;
                end
                e = est(sel); a = ase(sel); c = cp(sel);
                rows(end+1,:) = {M.scenario_index(i), margin, j, truth(j), ...
                    conditional, sum(sel), mean(e,"omitnan")-truth(j), ...
                    std(e,"omitnan"), mean(a(isfinite(a))), mean(isfinite(a)), ...
                    mean(c(isfinite(c))), mean((e-truth(j)).^2,"omitnan")}; %#ok<AGROW>
            end
        end
    end
end
if isempty(rows)
    return;
end
T = cell2table(rows,"VariableNames",{'scenario_index','margin','coefficient', ...
    'truth','conditional_on_support','n_used','bias','ese','ase', ...
    'ase_available','cp','mse'});
T.ase_over_ese = T.ase./T.ese;
writetable(T, fullfile(summary_dir, "inference_main.csv"));
end


function A = paired_ablation(out_root)
% Same seeds, working association at the truth (main) versus 0 (ablation).
Mm = plaph_v2_manifest("main");
Ma = plaph_v2_manifest("ablation_tau0");
rows = {};
for i = 1:height(Ma)
    j = find(Mm.scenario_index == Ma.scenario_index(i));
    Rm = read_scenario(fullfile(out_root,"main",char(Mm.scenario_id(j))));
    Ra = read_scenario(fullfile(out_root,"ablation_tau0",char(Ma.scenario_id(i))));
    if isempty(Rm) || isempty(Ra)
        continue;
    end
    [~,im,ia] = intersect(Rm.seed, Ra.seed);
    Rm = Rm(im,:); Ra = Ra(ia,:);
    metrics = ["imse_mpl_add_T" "imse_mpl_add_C"];
    for metric = metrics
        d = col(Ra,metric) - col(Rm,metric);  % >0: modelling dependence helps
        rows(end+1,:) = {Ma.scenario_index(i), Ma.n(i), Ma.p(i), ...
            Ma.target_censor_rate(i), Ma.tau(i), metric, numel(im), ...
            mean(col(Rm,metric),"omitnan"), mean(col(Ra,metric),"omitnan"), ...
            mean(d,"omitnan"), std(d,"omitnan")/sqrt(sum(isfinite(d))), ...
            mean(d > 0,"omitnan")}; %#ok<AGROW>
    end
    for comp = ["beta" "phi" "mpl_add_T" "mpl_add_C"]
        fp = comp + "_fp"; fn = comp + "_fn";
        em = double(col(Rm,fp) == 0 & col(Rm,fn) == 0);
        ea = double(col(Ra,fp) == 0 & col(Ra,fn) == 0);
        d = em - ea;
        rows(end+1,:) = {Ma.scenario_index(i), Ma.n(i), Ma.p(i), ...
            Ma.target_censor_rate(i), Ma.tau(i), "exact_" + comp, numel(im), ...
            mean(em), mean(ea), mean(d), std(d)/sqrt(numel(d)), mean(d > 0)}; %#ok<AGROW>
    end
end
A = cell2table(rows,"VariableNames",{'scenario_index','n','p','censor_target', ...
    'tau','metric','R_paired','value_true_tau','value_tau0','difference', ...
    'difference_mcse','share_favoring_true_tau'});
end


%% ----------------------------- LaTeX ------------------------------------
function write_main_tables(T, summary_dir)
principal = T(ismember(T.scenario_index,1:9),:);
principal = sortrows(principal,"scenario_index");
fid = fopen(fullfile(summary_dir,"table1_support_rows.tex"),"w");
for i = 1:height(principal)
    r = principal(i,:);
    fprintf(fid, "%.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %s & %s \\\\\n", ...
        r.tau, r.exact_beta, r.exact_gamma, r.exact_phiT, r.exact_phiC, ...
        r.pl_exact_beta, r.pl_exact_phiT, fmt(r.grp_exact_beta,"%.3f"), ...
        fmt(r.grp_exact_phiT,"%.3f"));
end
fclose(fid);
fid = fopen(fullfile(summary_dir,"table2_estimation_rows.tex"),"w");
for i = 1:height(principal)
    r = principal(i,:);
    fprintf(fid, "%.2f & %.4f & %.4f & %.4f & %.4f & %.4f & %s & %.4f & %.3f \\\\\n", ...
        r.tau, r.mse_beta, r.pl_mse_beta, r.mse_gamma, r.imse_phiT, ...
        r.pl_imse_phiT, fmt(r.grp_imse_phiT,"%.4f"), r.imse_phiC, r.pr_PLAPH_lt_PL);
end
fclose(fid);
dim = sortrows(T(ismember(T.scenario_index,[10 5 11]),:),"p");
fid = fopen(fullfile(summary_dir,"table3_dimension_rows.tex"),"w");
for i = 1:height(dim)
    r = dim(i,:);
    fprintf(fid, "%d & %.3f & %.3f & %.3f & %.3f & %.4f & %.4f & %s \\\\\n", ...
        r.p, r.exact_beta, r.exact_gamma, r.exact_phiT, r.exact_phiC, ...
        r.imse_phiT, r.pl_imse_phiT, fmt(r.grp_imse_phiT,"%.4f"));
end
fclose(fid);
fid = fopen(fullfile(summary_dir,"numerical_acceptance.tex"),"w");
fprintf(fid, "%% scenario & accept & selection gate & refit pass & max refit residual & sure screening\n");
for i = 1:height(T)
    r = T(i,:);
    fprintf(fid, "%d & %.3f & %.3f & %.3f & %.2e & %.3f \\\\\n", r.scenario_index, ...
        r.accept_rate, r.selection_gate_rate, r.refit_pass_rate, ...
        r.max_refit_residual, r.screen_all_sure);
end
fclose(fid);
end


function write_ablation_table(A, summary_dir)
fid = fopen(fullfile(summary_dir,"table_ablation_rows.tex"),"w");
for i = 1:height(A)
    r = A(i,:);
    fprintf(fid, "%d & %.2f & %s & %.4f & %.4f & %.4f (%.4f) & %.3f \\\\\n", ...
        r.n, r.tau, r.metric, r.value_true_tau, r.value_tau0, r.difference, ...
        r.difference_mcse, r.share_favoring_true_tau);
end
fclose(fid);
end


function s = fmt(x, f)
if isfinite(x)
    s = sprintf(f, x);
else
    s = "--";
end
end
