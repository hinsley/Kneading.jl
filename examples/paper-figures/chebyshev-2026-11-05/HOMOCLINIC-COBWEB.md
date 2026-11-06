# Figure 8: homoclinic and nearby chaotic cobwebs

Run `julia --project=. compute_cubic_homoclinic.jl` once, then
`julia --project=. render_cubic_critical_cobweb.jl` for cached-data style changes.
The data cache is `cubic_homoclinic_data.jld2`; numerical diagnostics are in
`cubic_homoclinic_validation.txt`. Use the PNG in Overleaf to preserve thin
red strokes. The paper includes it at single-column width.

At v=-1/2, solve f(f(u))=f(u) in u in (0.83,0.84), using 384-bit bisection.
Since f(-1/2)=u, x*=f(u) is the fixed point reached on iterate two.
The selected root is u=0.8365770164010115762693555543660572848012368816407477583...
The black marker is x*=0.0561629347546044306913047706746148116782598256047603...

The blue cobweb includes 80 backward iterates on the middle inverse branch,
x0=-1/2, x1=u, and the fixed forward tail. In particular x[-1]=1/2 exactly.
The inverse maps [-1/5,3/10] strictly into itself, with derivative magnitude
bounded by 1/((u+1/2)*(3/2-6*(3/10)^2)) < 1. Once the backward orbit enters
this interval, contraction proves convergence to the same fixed point.
The finite backward tail is within 1.8e-24 of the fixed point. All plotted
homoclinic transitions satisfy the map equation to better than 1e-110 before
conversion to plotting precision. The infinite fixed tail is represented
by the fixed point, not by an unstable floating-point forward integration.

The red cobweb shows 5000 iterates from -1/2+1e-6 at the same parameters.
Its finite-time Lyapunov exponent after 1000 transient iterates is about
0.3166, and the last 4001 samples occupy all 200 equal bins of [-1/2,u].
These diagnostics verify the intended numerical chaotic illustration;
they are not a proof of chaos for that particular rounded numerical seed.

The previous critical-only renderer is recoverable from repository history
(commit aaad283f9129f9a4f1607154c95bdf0ab5ba8d0a).
