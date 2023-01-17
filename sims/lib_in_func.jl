
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
