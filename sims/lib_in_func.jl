
# ω -> sine frequency
# lvl -> value

function in_sin(t, ω, A,  A1, A2, T1, T2)
    # ω: Sine frequency
    # A: Sine amplitude
    # A1: First step value
    # A2: Second step value
    # T1: Time to first step
    # T2: Time for the second step

    if t > T1
        if t > T2
            sin(ω * 2 * π * t) * A + A2
        else
            A2
        end
    else
        A1
    end
end

function in_sin2(t, ω, A,  A1, A2, T1, T2)
    # ω: Sine frequency
    # A: Sine amplitude
    # A1: First step value
    # A2: Second step value
    # T1: Time to first step
    # T2: Time for the second step

    if t > T1
        if t > T2
            sin(ω * 2 * π * t) * A + A2
        else
            min(A2*(t- T1)^2*0.75, A2)
        end
    else
        A1
    end
end

function in_sin3(t, ω, A,  A1, A2, T1, T2)
    # ω: Sine frequency
    # A: Sine amplitude
    # A1: First step value
    # A2: Second step value
    # T1: Time to first step
    # T2: Time for the second step

    if t > T1
        if t > T2
            sin(ω * 2 * π * t) * A + A2
        else
            smooth_startup_10(t, A2)
        end
    else
        A1
    end
end

function smooth_startup_10(t, A)
    # return min((A + A/100) ./(1 .+ 15.0.*exp.(-t.*1.4 .+ 4)), A)
    # 1280 ./ (1 .+ exp.(-t .* 2 .+ 15))
    return min((A + A/100) ./ (1 .+ exp.(-t .* 2 .+ 15)), A)
end