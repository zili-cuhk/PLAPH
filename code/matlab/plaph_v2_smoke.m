function plaph_v2_smoke(out_root, working_set_dir)
%PLAPH_V2_SMOKE End-to-end check of the v2 code before any large submission.
%   1. B-spline basis: order 4 gives 8 columns, centered rank 7; order 3 gives
%      7 columns, centered rank 6; partition of unity holds.
%   2. Two replications of scenario 5 (development seeds) with the full EBIC
%      grid for PLAPH and the aligned PL, with data export for grpsurv.
%   3. One TCGA fit (Frank, tau = 0, c_L = 1.5, kappa = 0.20).
%   4. Aggregation of the smoke output.
% Every step prints PASS/FAIL lines; the job log is the record.
out_root = char(out_root);
fprintf("==== PLAPH v2 smoke test, %s ====\n", string(datetime("now")));

% 1. basis check (independent of the engine's internal helper)
x = sort(rand(400,1));
for order = [4 3]
    knots_in = quantile(x, linspace(0,1,6)); knots_in = knots_in(2:5);
    t = [min(x)*ones(1,order), knots_in, max(x)*ones(1,order)];
    nb = numel(knots_in) + order;
    Bx = zeros(numel(x), nb);
    for j = 1:nb
        Bx(:,j) = cox_de_boor(j-1, order, t, x);
    end
    pu = max(abs(sum(Bx,2) - 1));
    r = rank(Bx - mean(Bx,1), 1e-10);
    fprintf("[%s] order %d: %d columns, centered rank %d (expected %d), partition-of-unity error %.1e\n", ...
        passfail(r == nb-1 && pu < 1e-12), order, nb, r, nb-1, pu);
end

% 2. simulation smoke (2 replications, EBIC grid, PL, export)
t0 = tic;
plaph_v2_run_task("smoke", 1, out_root, "ChunkSize", 2, ...
    "ExportRoot", fullfile(out_root, "grpsurv_data"));
M = plaph_v2_manifest("smoke");
f = dir(fullfile(out_root, "smoke", char(M.scenario_id(1)), "chunk_*", "scenario_001_replicates.csv"));
ok = ~isempty(f);
fprintf("[%s] simulation chunk finished in %.1f min\n", passfail(ok), toc(t0)/60);
if ok
    R = readtable(fullfile(f(1).folder, f(1).name));
    show = intersect({'seed','ok','mpl_accepted','mpl_selection_gate_pass', ...
        'mpl_selection_status_code','mpl_selection_gate_maximum', ...
        'mpl_refit_stationarity_inf','mpl_target_iterations','mpl_step_doublings', ...
        'tuning_cL','tuning_kappa','tuning_tie_count','tuning_at_boundary', ...
        'pl_tuning_cL','pl_tuning_kappa','beta_fp','beta_fn','phi_fp','phi_fn', ...
        'mpl_add_T_fp','mpl_add_T_fn','mpl_add_C_fp','mpl_add_C_fn', ...
        'imse_mpl_add_T','imse_mpl_add_C','pl_accepted','pl_add_fp','pl_add_fn', ...
        'imse_pl_add','spline_order'}, R.Properties.VariableNames, 'stable');
    disp(R(:, show));
    fprintf("[%s] all replications ok\n", passfail(all(R.ok == 1)));
end
g = dir(fullfile(out_root, "grpsurv_data", "*", "data_seed*.csv.gz"));
fprintf("[%s] exported %d data sets for grpsurv\n", passfail(numel(g) == 2), numel(g));

% 3. TCGA single fit
if nargin > 1 && strlength(string(working_set_dir)) > 0
    t0 = tic;
    try
        row = plaph_v2_tcga_task("grid_Frank_0.00", working_set_dir, ...
            fullfile(out_root, "tcga_smoke"), "GridCL", 1.5, "GridKappa", 0.20);
        disp(row);
        fprintf("[PASS] TCGA fit finished in %.1f min (accepted = %d)\n", ...
            toc(t0)/60, row.accepted);
    catch ME
        fprintf("[FAIL] TCGA fit: %s\n", ME.message);
        disp(getReport(ME));
    end
end

% 4. aggregation
try
    S = plaph_v2_aggregate(out_root, "Experiments", "smoke");
    if isfield(S, "smoke")
        disp(S.smoke(:, {'R_done','accept_rate','exact_beta','exact_gamma', ...
            'exact_phiT','exact_phiC','imse_phiT','pl_imse_phiT'}));
    end
    fprintf("[PASS] aggregation\n");
    fprintf("[%s] tuning optimal-share table written\n", passfail(isfile( ...
        fullfile(out_root, "summary", "tuning_optimal_share_smoke.csv"))));
catch ME
    fprintf("[FAIL] aggregation: %s\n", ME.message);
end
fprintf("==== smoke test finished ====\n");
end


function y = cox_de_boor(j, k, t, x)
% Independent Cox-de Boor recurrence (order k), right end included.
if k == 1
    y = double(t(j+1) <= x & x < t(j+2));
    if t(j+2) == t(end) && t(j+1) < t(j+2)
        y(x == t(end)) = 1;
    end
    return;
end
y = zeros(size(x));
d1 = t(j+k) - t(j+1);
if d1 > 0
    y = y + (x - t(j+1))/d1 .* cox_de_boor(j, k-1, t, x);
end
d2 = t(j+k+1) - t(j+2);
if d2 > 0
    y = y + (t(j+k+1) - x)/d2 .* cox_de_boor(j+1, k-1, t, x);
end
end


function s = passfail(c)
if c
    s = "PASS";
else
    s = "FAIL";
end
end
