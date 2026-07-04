/* coremath.h — declarations for the vendored CORE-MATH correctly-
 * rounded binary64 functions. Emitted podium kernels built in
 * correctly-rounded mode call these in place of libm sin/cos, so their
 * results are the unique nearest double to the true value and agree
 * bit-for-bit with any correctly-rounded reference (see
 * podium.emit.croracle). CORE-MATH is MIT-licensed; see LICENSE. */
#ifndef PODIUM_COREMATH_H
#define PODIUM_COREMATH_H

double cr_sin(double x);
double cr_cos(double x);

#endif
