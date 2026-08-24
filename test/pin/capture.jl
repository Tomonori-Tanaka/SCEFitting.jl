# Regenerate the pins.  MANUAL: this is never run by the suite -- a pin that
# rewrites itself detects nothing.  See PIN.md for when recapture is allowed.
#
#   julia --project=test/pin -t 4 test/pin/capture.jl
include(joinpath(@__DIR__, "payload.jl"))
include(joinpath(@__DIR__, "fixtures.jl"))
mkpath(joinpath(@__DIR__, "pins"))

# `PIN_ONLY=id1,id2` regenerates just those fixtures.  Without it this script
# rewrites EVERY pin, including ones the change under way does not touch -- and a
# rewritten pin detects nothing, so recapturing a fixture is a deliberate act that
# needs its own entry in PIN.md (rule 5).
const _ONLY = let v = get(ENV, "PIN_ONLY", "")
    isempty(v) ? nothing : Set(strip.(split(v, ",")))
end
_wanted(id) = _ONLY === nothing || id in _ONLY

for (fx, payload) in vcat([(fx, pin_payload) for fx in PIN_FIXTURES],
                          [(fx, pin_moment_payload) for fx in PIN_MOMENT_FIXTURES])
    _wanted(fx.id) || continue
    d = payload(fx)
    d["meta"] = Dict{String,Any}(
        # Both are passed in, not sniffed: recapturing is a deliberate act, and
        # the entry in PIN.md that explains WHY is written at the same moment.
        "captured" => get(ENV, "PIN_DATE", "SET-ME"),
        "julia" => string(VERSION),
        "platform" => string(Sys.ARCH, "-", Sys.KERNEL),
        "threads" => Threads.nthreads(),
        "commit" => get(ENV, "PIN_COMMIT", "SET-ME"))
    open(joinpath(@__DIR__, "pins", fx.id * ".toml"), "w") do io
        TOML.print(io, d; sorted = true)
    end
    println("captured ", fx.id, "  n_salcs = ", d["L0"]["n_salcs"],
            "  folded = ", length(d["L1"]["folded"]))
end
