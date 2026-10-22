function decision = plaph_v2_spline_decision(out_root)
%PLAPH_V2_SPLINE_DECISION Prespecified choice between cubic and quadratic.
%
% Development seeds only (scenarios 1, 4, 6; R = 20 each). Rule fixed before
% the pilot is run: keep the cubic basis (order 4, consistent with the
% manuscript and the approximation theory) unless the quadratic basis is
% better on BOTH criteria by a material margin:
%   (a) mean exact-recovery rate over the four supports and three scenarios
%       is higher by at least 0.03, and
%   (b) mean of the failure and censoring IMSE is lower by at least 10%.
% Either basis must have an acceptance rate of at least 0.95; otherwise the
% engine must be fixed before the formal campaign.
out_root = char(out_root);
summ = struct();
for order = [4 3]
    root = fullfile(out_root, sprintf("pilot_order%d", order));
    S = plaph_v2_aggregate(root, "Experiments", "spline_pilot");
    if ~isfield(S, "spline_pilot")
        error("PLAPH:PilotMissing", "No pilot results under %s.", root);
    end
    T = S.spline_pilot;
    summ.(sprintf("order%d",order)) = struct( ...
        'recovery', mean(mean([T.exact_beta T.exact_gamma T.exact_phiT T.exact_phiC],2)), ...
        'imse', mean(mean([T.imse_phiT T.imse_phiC],2)), ...
        'accept', mean(T.accept_rate), 'table', T);
end
c = summ.order4; qd = summ.order3;
better_recovery = qd.recovery >= c.recovery + 0.03;
better_imse = qd.imse <= 0.90*c.imse;
if better_recovery && better_imse
    chosen = 3;
else
    chosen = 4;
end
decision = struct('chosen_order', chosen, ...
    'cubic_recovery', c.recovery, 'quadratic_recovery', qd.recovery, ...
    'cubic_imse', c.imse, 'quadratic_imse', qd.imse, ...
    'cubic_accept', c.accept, 'quadratic_accept', qd.accept, ...
    'acceptance_ok', c.accept >= 0.95 && qd.accept >= 0.95);
disp(decision);
writetable(struct2table(rmfield(decision, {})), ...
    fullfile(out_root, "spline_decision.csv"));
if ~decision.acceptance_ok
    warning("PLAPH:PilotAcceptance", ...
        "Acceptance below 0.95 in the pilot: inspect the logs before the formal campaign.");
end
end
