function V = plaph_survival_uniforms(family,U)
%PLAPH_SURVIVAL_UNIFORMS Convert a copula draw to latent survival probabilities.
%
% The event-time construction uses T = H_0^{-1}(-log(V)/exp(eta)).
% Therefore a target survival copula C requires (V_T,V_C) to have copula C.
% For radially symmetric Frank and Gaussian copulas, 1-U has the same
% copula as U; retaining 1-U preserves the historical random realization.
% Clayton and Gumbel are not radially symmetric, so U itself must be used.

    arguments
        family {mustBeTextScalar}
        U (:,2) double
    end

    if any(~isfinite(U),"all") || any(U <= 0 | U >= 1,"all")
        error("plaph_survival_uniforms:InvalidDraw", ...
            "Copula draws must be finite and lie strictly between zero and one.");
    end

    family_l = lower(string(family));
    switch family_l
        case {"frank","gaussian"}
            V = 1-U;
        case {"clayton","gumbel"}
            V = U;
        otherwise
            error("plaph_survival_uniforms:UnsupportedFamily", ...
                "Unsupported copula family: %s",family_l);
    end
end
