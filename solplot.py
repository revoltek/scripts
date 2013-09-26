import math
import numpy
import pylab
from solfetch import fetch

styles = ["%s%s" % (x, y) for y in ["-", ":"] for x in ["b", "g", "r", "c", "m",
    "y", "k"]]

def unwrap(phase, tol=0.25, delta_tol=0.25):
    """
    Unwrap phase by restricting phase[n] to fall within a range [-tol, tol]
    around phase[n - 1].

    If this is impossible, the closest phase (modulo 2*pi) is used and tol is
    increased by delta_tol (tol is capped at pi).
    """

    assert(tol < math.pi)

    # Allocate result.
    out = numpy.zeros(phase.shape)

    # Effective tolerance.
    eff_tol = tol

    ref = phase[0]
    for i in range(0, len(phase)):
        delta = math.fmod(phase[i] - ref, 2.0 * math.pi)

        if delta < -math.pi:
            delta += 2.0 * math.pi
        elif delta > math.pi:
            delta -= 2.0 * math.pi

        out[i] = ref + delta

        if abs(delta) <= eff_tol:
            # Update reference phase and reset effective tolerance.
            ref = out[i]
            eff_tol = tol
        elif eff_tol < math.pi:
            # Increase effective tolerance.
            eff_tol += delta_tol * tol
            if eff_tol > math.pi:
                eff_tol = math.pi

    return out

def unwrap_windowed(phase, window_size=5):
    """
    Unwrap phase by estimating the trend of the phase signal.
    """

    # Allocate result.
    out = numpy.zeros(phase.shape)

    windowl = numpy.array([math.fmod(phase[0], 2.0 * math.pi)] * window_size)

    delta = math.fmod(phase[1] - windowl[0], 2.0 * math.pi)
    if delta < -math.pi:
        delta += 2.0 * math.pi
    elif delta > math.pi:
        delta -= 2.0 * math.pi
    windowu = numpy.array([windowl[0] + delta] * window_size)

    out[0] = windowl[0]
    out[1] = windowu[0]

    meanl = windowl.mean()
    meanu = windowu.mean()
    slope = (meanu - meanl) / float(window_size)

    for i in range(2, len(phase)):
        ref = meanu + (1.0 + (float(window_size) - 1.0) / 2.0) * slope
        delta = math.fmod(phase[i] - ref, 2.0 * math.pi)

        if delta < -math.pi:
            delta += 2.0 * math.pi
        elif delta > math.pi:
            delta -= 2.0 * math.pi

        out[i] = ref + delta

        windowl[:-1] = windowl[1:]
        windowl[-1] = windowu[0]
        windowu[:-1] = windowu[1:]
        windowu[-1] = out[i]

        meanl = windowl.mean()
        meanu = windowu.mean()
        slope = (meanu - meanl) / float(window_size)

    return out

def phase_normalize(phase):
    """
    Normalize input phase to the range [-pi, pi].
    """

    # Convert to range [-2*pi, 2*pi].
    out = numpy.fmod(phase, 2.0 * numpy.pi)

    # Convert to range [-pi, pi]
    out[out < -numpy.pi] += 2.0 * numpy.pi
    out[out > numpy.pi] -= 2.0 * numpy.pi

    return out

def plot(sol, fig=None, sub=None, scatter=False, stack=False, sep=5.0,
    sep_abs=False, labels=None, show_legend=False):
    """
    Plot a list of signals.

    If 'fig' is equal to None, a new figure will be created. Otherwise, the
    specified figure number is used. The 'sub' argument can be used to create
    subplots. The 'sep' and 'stack' argument can be used to control placement
    of the plots in the list.

    The figure number of the figure used to plot in is returned.
    """
    global styles

    cf = pylab.figure(fig)

    if fig is not None and sub is None:
        pylab.clf()

    if sub is not None:
        pylab.subplot(sub)

    offset = 0.0
    for i in range(0,len(sol)):
        if labels is None:
            if scatter:
                pylab.scatter(range(0, len(sol[i])), sol[i] + offset,
                    edgecolors="None", c=styles[i % len(styles)][0], marker="o")
            else:
                pylab.plot(sol[i] + offset, styles[i % len(styles)])
        else:
            if scatter:
                pylab.scatter(range(0, len(sol[i])), sol[i] + offset,
                    edgecolors="None", c=styles[i % len(styles)][0], marker="o",
                    label=labels[i])
            else:
                pylab.plot(sol[i] + offset, styles[i % len(styles)],
                    label=labels[i])

        if stack:
            if sep_abs:
                offset += sep
            else:
                offset += sol[i].mean() + sep * sol[i].std()

    if not labels is None and show_legend:
        pylab.legend()

    return cf.number
