#!/bin/bash
# Tier-2 cross-architecture golden vectors. Runs INSIDE a debian
# container: installs the aarch64 cross toolchain + qemu user emulation,
# cross-compiles the emitted kernels + golden driver STATICALLY for
# aarch64, and replays each kernel's recorded input vectors under
# qemu-aarch64-static. IEEE-754 binary64 with -ffp-contract=off is
# architecture-independent, so the outputs must match the x86/CPython
# reference kernel-for-kernel (bit-exact for arithmetic+sqrt; the
# libm-trig tolerance is the same cross-libm one as on x86, and it
# vanishes in the correctly-rounded build).
#
# /work layout (written by tools/tier2_qemu.py):
#   kernels.c driver.c        default-mode translation unit + driver
#   kernels.txt               kernel names, one per line
#   <name>.in                 hex-float input vectors per kernel
#   coremath/                 (optional) CR-mode: sin.c cos.c coremath.h
#   kernels_cr.c driver_cr.c  (optional) CR-mode stm unit + driver
# Outputs: <name>.out and (optional) stm_cr.out
set -e
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq gcc-aarch64-linux-gnu qemu-user-static >/dev/null

CC=aarch64-linux-gnu-gcc
FLAGS="-std=c99 -O2 -ffp-contract=off -static"

$CC $FLAGS /work/kernels.c /work/driver.c -lm -o /work/kern
while read -r name; do
    qemu-aarch64-static /work/kern "$name" < "/work/$name.in" \
        > "/work/$name.out"
done < /work/kernels.txt

if [ -f /work/kernels_cr.c ]; then
    $CC $FLAGS -I /work/coremath /work/kernels_cr.c /work/driver_cr.c \
        /work/coremath/sin.c /work/coremath/cos.c -lm -o /work/kern_cr
    qemu-aarch64-static /work/kern_cr stm < /work/stm_cr.in \
        > /work/stm_cr.out
fi
echo TIER2_OK
