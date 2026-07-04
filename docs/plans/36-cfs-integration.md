# 36 — cFS integration example: verified kernels on a software bus

GitHub issue: https://github.com/adi-oltean/podium/issues/36
(GitHub #35 was a mis-numbered first filing, closed and refiled as #36
so title == number; the #35 slot had been used by the viewer
frame-blending commit, which shipped without a tracking issue.)

## Fix (landed)

`podium.emit.cfsapp.generate(outdir, kernels_c)` writes a Core Flight
System relative-navigation app + a portable shim + README; the
generated tree is also checked in at `examples/cfs_nav_app/`.

- `podium_nav_app.c/.h`: the app in genuine cFS idiom — CFE_ES_RunLoop
  main loop, CFE_SB CreatePipe/Subscribe/ReceiveBuffer/TransmitMsg,
  CFE_EVS_SendEvent, CFE_MSG headers, message IDs. The ONLY arithmetic
  is the measurement handler, which calls the emitted
  podium_ekf_update_sequential then podium_ekf_predict over a
  maintained 6-state estimate and 6x6 covariance — so the flight
  numerics carry the full verification stack (EVA-proven,
  CompCert-checked, golden-vector bit-exact).
- `cfe_shim.c/.h`: a minimal, portable stand-in for the cFE SB/ES/EVS
  API backed by a stdin/stdout hex-float message transport (the same
  %a transport the golden vectors use). NOT cFE — it lets the same app
  source run standalone without OSAL/full cFE. Swap for real cfe.h in
  a mission build unchanged.
- README: maps every piece to real cFS AND to F´ (component, input/
  output ports, parameters, commands, telemetry).

## Receipts (tests/test_cfs_app.py)

- The compiled app (gcc -Wall -Werror clean) reproduces the Python EKF
  reference to 1e-11 (scaled) on a 60-step recorded nav stream — same
  init, same transport, only the emitted predict's BLAS-vs-naive
  matmul order differing (the documented test_cemit class).
- The app announces init on EVS and processes a stream to completion:
  the cFS run loop + SB plumbing work end to end.

## The tie-back worth noting

The app runs r_var at the EVA-proven envelope FLOOR (0.01 = (10 cm)^2).
A tighter sensor variance would sit below the region the sound-value
gate certified overflow-free in #30. So the verification constrains
the flight CONFIGURATION, not merely the code — the proof envelope is
an operational fact the app respects.

## Deferred

Real cFS mission-tree build in CI (needs the external cFE/OSAL);
F´ topology XML; command ingest (reset/reinit); multi-app bus demo.

## Push/merge instructions

Single commit on main: `36 — cFS integration example (#36)`; push;
close. Marks v0.6 COMPLETE.
