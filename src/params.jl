struct Params

    log_polynomial_length :: Int # in HEAAN: N
    log_lo_modulus :: Int # in HEAAN: logQ
    log_hi_modulus :: Int # in HEAAN: == logQ
    gaussian_noise_stddev :: Float64 # in HEAAN: sigma
    secret_key_length :: Int # in HEAAN: h

    function Params(; log_polynomial_length::Int=16, log_lo_modulus::Int=300)
        # TODO: this all (including log_polynomial_length) should be derived
        # from the security parameter `lambda`.

        # From Cheon et al.,
        # "We need to set the ring dimension N that satisfies the security condition
        # `N >= (λ+110) / 7.2 * (log_hi_modulus + log_lo_modulus)` ...
        # it suffices to assume that P [`log_hi_modulus`] is
        # approximately equal to qL [`log_lo_modulus`]."

        gaussian_noise_stddev = 3.2
        secret_key_length = 64
        log_hi_modulus = log_lo_modulus

        new(
            log_polynomial_length,
            log_lo_modulus,
            log_hi_modulus,
            gaussian_noise_stddev,
            secret_key_length,
            )
    end
end

