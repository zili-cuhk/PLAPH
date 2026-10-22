function plaph_v2_curve_summary(out_root, varargin)
%PLAPH_V2_CURVE_SUMMARY Function-estimation summaries from the stored curves.
%
%   plaph_v2_curve_summary("/path/to/plaph_v2_runs")
%
% For every scenario of the listed experiments, reads the fitted additive
% curves stored in each chunk's PLAPH_simulation_results.mat
% (CollectFunctionFits = true) and writes to <out_root>/summary:
%   function_error_<exp>.csv  scenario, method, component: ISB^2, IVar, IMSE
%   curves_<exp>.csv          grid, truth, Monte Carlo mean, 2.5%/97.5% bands
% Curves are re-centered with normalized trapezoidal weights, and ISB^2,
% IVar and IMSE are trapezoidal integrals over the 101-point grid, matching
% the per-replication IMSE of the engine (IMSE = ISB^2 + IVar). All
% replications are used; a component that is not selected is the zero curve.

parser = inputParser;
parser.addParameter("Experiments", ["main" "ablation_tau0" "highcorr" "copula_family"]);
parser.parse(varargin{:});
experiments = string(parser.Results.Experiments);
out_root = char(out_root);
summary_dir = fullfile(out_root, "summary");
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end
methods = ["mpl_T" "mpl_C" "pl"];
for experiment = experiments
    exp_dir = fullfile(out_root, char(experiment));
    if ~isfolder(exp_dir)
        warning("PLAPH:CurveSummary", "No folder %s; skipped.", exp_dir);
        continue;
    end
    M = plaph_v2_manifest(experiment);
    err_rows = {};
    curve_rows = {};
    for i = 1:height(M)
        files = dir(fullfile(exp_dir, char(M.scenario_id(i)), "chunk_*", ...
            "PLAPH_simulation_results.mat"));
        if isempty(files)
            continue;
        end
        store = struct();
        for m = methods
            store.(m) = [];
        end
        grid = []; truth_T = []; truth_C = [];
        for f = 1:numel(files)
            S = load(fullfile(files(f).folder, files(f).name), "sim");
            ff = S.sim.results{1}.function_fits;
            if isempty(grid)
                grid = ff.grid; truth_T = ff.true_T; truth_C = ff.true_C;
            end
            for m = methods
                if isfield(ff, m)
                    store.(m) = cat(3, store.(m), ff.(m));
                end
            end
        end
        w = trapezoid_weights(size(grid,1));
        for m = methods
            F = store.(m);
            if isempty(F) || all(isnan(F(:)))
                continue;
            end
            if m == "mpl_C"
                truth = truth_C;
            else
                truth = truth_T;
            end
            for j = 1:size(grid,2)
                fj = squeeze(F(:,j,:));                       % grid x R
                keep = all(isfinite(fj),1);
                fj = fj(:,keep);
                if isempty(fj)
                    continue;
                end
                fj = fj - sum(w.*fj,1);                       % trapezoid centering
                tj = truth(:,j) - sum(w.*truth(:,j));
                mean_curve = mean(fj,2);
                isb2 = sum(w.*(mean_curve - tj).^2);
                ivar = sum(w.*mean((fj - mean_curve).^2,2));
                imse = mean(sum(w.*(fj - tj).^2,1));
                err_rows{end+1} = table(M.scenario_index(i), m, j, size(fj,2), ...
                    isb2, ivar, imse, 'VariableNames', {'scenario_index', ...
                    'method','component','n_curves','isb2','ivar','imse'}); %#ok<AGROW>
                lo = quantile(fj, 0.025, 2);
                hi = quantile(fj, 0.975, 2);
                curve_rows{end+1} = table(repmat(M.scenario_index(i),size(grid,1),1), ...
                    repmat(m,size(grid,1),1), repmat(j,size(grid,1),1), grid(:,j), ...
                    tj, mean_curve, lo, hi, 'VariableNames', {'scenario_index', ...
                    'method','component','x','truth','mean','lo','hi'}); %#ok<AGROW>
            end
        end
        fprintf("%s scenario %d: %d chunks\n", experiment, M.scenario_index(i), numel(files));
    end
    if ~isempty(err_rows)
        writetable(vertcat(err_rows{:}), fullfile(summary_dir, "function_error_" + experiment + ".csv"));
        writetable(vertcat(curve_rows{:}), fullfile(summary_dir, "curves_" + experiment + ".csv"));
    end
end
end


function w = trapezoid_weights(m)
w = ones(m,1);
w([1 m]) = 0.5;
w = w/sum(w);
end
