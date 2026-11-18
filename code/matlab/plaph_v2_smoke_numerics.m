function plaph_v2_smoke_numerics(out_root, tmp_root)
%PLAPH_V2_SMOKE_NUMERICS Check the LegacyNumerics switch (2026-11-17).
%   Refits replication 1 of main D6 (n = 500, p = q = 400, tau = 0.8) twice:
%   with LegacyNumerics = true, which must reproduce the stored result, and
%   with the corrected numerics. Prints the differences from the stored run.
%   plaph_v2_smoke_numerics(OUT_ROOT, fullfile(OUT_ROOT,"smoke_numerics"))
M = plaph_v2_manifest("main");
row = find(M.scenario_index == 6);
task = (row-1)*M.R(row) + 1;            % ChunkSize 1: replication 1 of D6
seed = M.seed_root(row) + 1000*M.seed_scenario(row) + 1;
stored = read_reps(fullfile(out_root, "main", M.scenario_id(row)));
stored = stored(stored.seed == seed, :);
if isempty(stored)
    error("Stored replication with seed %d not found under %s.", seed, out_root);
end
cols = ["beta_fp" "beta_fn" "phi_fp" "phi_fn" "mpl_add_T_fp" "mpl_add_C_fp" ...
    "mpl_accepted" "tuning_cL" "tuning_kappa" ...
    "mpl_beta_active1" "mpl_beta_active2" "mpl_beta_active3" ...
    "mpl_phi_active1" "mpl_phi_active2" "mpl_phi_active3" ...
    "mpl_beta_ase_active1" "mpl_phi_ase_active1" "imse_mpl_add_T" "imse_mpl_add_C"];
cols = cols(ismember(cols, string(stored.Properties.VariableNames)));
for legacy = [true false]
    root = fullfile(char(tmp_root), sprintf("legacy%d", legacy));
    t0 = tic;
    plaph_v2_run_task("main", task, root, "ChunkSize", 1, "SplineOrder", 3, ...
        "LegacyNumerics", legacy);
    new = read_reps(fullfile(root, "main", M.scenario_id(row)));
    new = new(new.seed == seed, :);
    fprintf("\nLegacyNumerics = %d (%.1f min)\n", legacy, toc(t0)/60);
    for c = cols
        a = stored.(c); b = new.(c);
        fprintf("  %-24s stored %12.6g   new %12.6g   diff %10.3g\n", c, a, b, b - a);
    end
end
end

function T = read_reps(scenario_dir)
% Results of the first chunk (it contains replication 1).
f = dir(fullfile(char(scenario_dir), "chunk_001_*", "scenario_001_replicates.csv"));
if isempty(f)
    error("No chunk_001 results under %s.", scenario_dir);
end
T = readtable(fullfile(f(1).folder, f(1).name));
end
