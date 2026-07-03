"""Guidance laws and trajectory optimization for proximity operations.

Modules
-------
glideslope  Classical glideslope approach guidance (Hablani), core-compliant.
targeting   Multi-impulse CW/TH targeting wrappers around podium.core.cw.
convex      Layer-0 convex trajectory optimization: DPP-compiled LP/SOCP
            transcription on exact CW/YA/ROE discretizations, approach
            cone, hyperplane KOZ, plume, passive-safety scenarios, LCvx
            finite burns with validity audits, MIB quantization bridge.
scp         Layer-1 PTR/SCvx* successive convexification: true nonconvex
            keep-out constraints with trust regions, virtual buffers,
            penalty ramp, and exact-flow continuous-time cuts.
safety      Passive-abort / free-drift safety metrics (e/i separation,
            RN-plane minimum separation: bounded scan + analytic form).
tumbling    Tumbling-target terminal guidance (scoped study): known-
            tumble port capture stays CONVEX (rotating corridor as
            per-node cones, terminal port match as a boundary state).
arch        ARCH-COMP rendezvous benchmark model + hybrid-automaton
            export for the CI reachability gate.
"""
