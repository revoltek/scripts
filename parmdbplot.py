#!/usr/bin/python
#
# Copyright (C) 2007
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id$

import sys
import lofar.parmdb as parmdb
import copy
import math
import numpy

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties

from PyQt4.QtCore import *
from PyQt4.QtGui import *

__styles = ["%s%s" % (x, y) for y in ["-", ":"] for x in ["b", "g", "r", "c",
    "m", "y", "k"]]

def contains(container, item):
    try:
        return container.index(item) >= 0
    except ValueError:
        return False

def common_domain(parms):
    if len(parms) == 0:
        return None

    domain = [-1e30, 1e30, -1e30, 1e30]
    for parm in parms:
        tmp = parm.domain()
        domain = [max(domain[0], tmp[0]), min(domain[1], tmp[1]), max(domain[2], tmp[2]), min(domain, tmp[3])]

    if domain[0] >= domain[1] or domain[2] >= domain[3]:
        return None

    return domain

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

def normalize(phase):
    """
    Normalize phase to the range [-pi, pi].
    """

    # Convert to range [-2*pi, 2*pi].
    out = numpy.fmod(phase, 2.0 * numpy.pi)

    # Convert to range [-pi, pi]
    out[out < -numpy.pi] += 2.0 * numpy.pi
    out[out > numpy.pi] -= 2.0 * numpy.pi

    return out

def plot(fig, y, x=None, clf=True, sub=None, scatter=False, stack=False,
    sep=5.0, sep_abs=False, labels=None, show_legend=False, title=None,
    xlabel=None, ylabel=None):
    """
    Plot a list of signals.

    If 'fig' is equal to None, a new figure will be created. Otherwise, the
    specified figure number is used. The 'sub' argument can be used to create
    subplots.

    The 'scatter' argument selects between scatter and line plots.

    The 'stack', 'sep', and 'sep_abs' arguments can be used to control placement
    of the plots in the list. If 'stack' is set to True, each plot will be
    offset by the mean plus sep times the standard deviation of the previous
    plot. If 'sep_abs' is set to True, 'sep' is used as is.

    The 'labels' argument can be set to a list of labels and 'show_legend' can
    be set to True to show a legend inside the plot.

    The figure number of the figure used to plot in is returned.
    """
    global __styles

    if clf:
        fig.clf()

    if sub is not None:
        fig.add_subplot(sub)

    axes = fig.gca()
    if not title is None:
        axes.set_title(title)
    if not xlabel is None:
        axes.set_xlabel(xlabel)
    if not ylabel is None:
        axes.set_ylabel(ylabel)

    if x is None:
        x = [range(len(yi)) for yi in y]

    offset = 0.0
    for i in range(0,len(y)):
        if labels is None:
            if scatter:
                axes.scatter(x[i], y[i] + offset, edgecolors="None",
                    c=__styles[i % len(__styles)][0], marker="o")
            else:
                axes.plot(x[i], y[i] + offset, __styles[i % len(__styles)])
        else:
            if scatter:
                axes.scatter(x[i], y[i] + offset, edgecolors="None",
                    c=__styles[i % len(__styles)][0], marker="o",
                    label=labels[i])
            else:
                axes.plot(x[i], y[i] + offset, __styles[i % len(__styles)],
                    label=labels[i])

        if stack:
            if sep_abs:
                offset += sep
            else:
                offset += y[i].mean() + sep * y[i].std()

    if not labels is None and show_legend:
        axes.legend(prop=FontProperties(size="x-small"), markerscale=0.5)

class Parm:
    def __init__(self, db, name, elements=None, isPolar=False):
        self._db = db
        self._name = name
        self._elements = elements
        self._isPolar = isPolar
        self._value = None
        self._value_domain = None
        self._value_resolution = None

        self._readDomain()

    def name(self):
        return self._name

    def isPolar(self):
        return self._isPolar

    def empty(self):
        return self._empty

    def domain(self):
        return self._domain

    def value(self, domain=None, resolution=None, asPolar=True, unwrap_phase=False, phase_reference=None):
        if self.empty():
            return (numpy.zeros((1,1)), numpy.zeros((1,1)))

        if self._value is None or self._value_domain != domain or self._value_resolution != resolution:
            self._readValue(domain, resolution)

            # Correct negative amplitude solutions by taking the absolute value
            # of the amplitude and rotating the phase by 180 deg.
            if self.isPolar():
                self._value[1][self._value[0] < 0.0] += numpy.pi
                self._value[0] = numpy.abs(self._value[0])

        if asPolar:
            if self.isPolar():
                ampl = self._value[0]
                phase = normalize(self._value[1])
            else:
                ampl = numpy.sqrt(numpy.power(self._value[0], 2) + numpy.power(self._value[1], 2))
                phase = numpy.arctan2(self._value[1], self._value[0])

            if not phase_reference is None:
                assert(phase_reference.shape == phase.shape)
                phase = normalize(phase - phase_reference)

            if unwrap_phase:
                for i in range(0, phase.shape[1]):
                    phase[:, i] = unwrap(phase[:, i])

            return [ampl, phase]

        if not self.isPolar():
            re = self._value[0]
            im = self._value[1]
        else:
            re = self._value[0] * numpy.cos(self._value[1])
            im = self._value[0] * numpy.sin(self._value[1])

        return [re, im]

    def _readDomain(self):
        if self._elements is None:
            self._domain = self._db.getRange(self.name())
        else:
            if self._elements[0] is None:
                self._domain = self._db.getRange(self._elements[1])
            elif self._elements[1] is None:
                self._domain = self._db.getRange(self._elements[0])
            else:
                domain_el0 = self._db.getRange(self._elements[0])
                domain_el1 = self._db.getRange(self._elements[1])
                self._domain = [max(domain_el0[0], domain_el1[0]), min(domain_el0[1], domain_el1[1]), max(domain_el0[2], domain_el1[2]), min(domain_el0[3], domain_el1[3])]

        self._empty = (self._domain[0] >= self._domain[1]) or (self._domain[2] >= self._domain[3])

    def _readValue(self, domain=None, resolution=None):
        if self._elements is None:
            value = numpy.array(self.__fetch_value(self.name(), domain, resolution))
            self._value = (value, numpy.zeros(value.shape))
        else:
            el0 = None
            if not self._elements[0] is None:
                el0 = numpy.array(self.__fetch_value(self._elements[0], domain, resolution))

            el1 = None
            if not self._elements[1] is None:
                el1 = numpy.array(self.__fetch_value(self._elements[1], domain, resolution))

            assert((not el0 is None) or (not el1 is None))

            if el0 is None:
                el0 = numpy.zeros(el1.shape)

            if el1 is None:
                el1 = numpy.zeros(el0.shape)

            self._value = [el0, el1]

        self._value_domain = domain
        self._value_resolution = resolution

    def __fetch_value(self, name, domain=None, resolution=None):
        if domain is None:
            tmp = self._db.getValuesGrid(name)[name]
        else:
            if resolution is None:
                tmp = self._db.getValues(name, domain[0], domain[1], domain[2], domain[3])[name]
            else:
                tmp = self._db.getValuesStep(name, domain[0], domain[1], resolution[0], domain[2], domain[3], resolution[1])[name]

        if type(tmp) is dict:
            return tmp["values"]

        # Old parmdb interface.
        return tmp

class PlotWindow(QFrame):
    def __init__(self, parms, selection, resolution=None, parent=None, title=None):
        QFrame.__init__(self, parent)

        if not title is None:
            self.setWindowTitle(title)

        self.parms = parms
        self.selected_parms = [self.parms[i] for i in selection]

        self.fig = Figure((5, 4), dpi=100)

        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        self.toolbar = NavigationToolbar(self.canvas, self)

        self.axis = 0
        self.index = 0
        self.axisSelector = QComboBox()
        self.axisSelector.addItem("Frequency")
        self.axisSelector.addItem("Time")
        self.connect(self.axisSelector, SIGNAL('activated(int)'), self.handle_axis)

        self.spinner = QSpinBox()
        self.connect(self.spinner, SIGNAL('valueChanged(int)'), self.handle_spinner)

#        self.slider = QSlider(Qt.Horizontal)
#        self.slider.setMinimum(0)
#        self.slider.setMaximum(159)
#        self.connect(self.slider, SIGNAL('sliderReleased()'), self.handle_slider)

        self.show_legend = False
        self.legendCheck = QCheckBox("Legend")
        self.connect(self.legendCheck, SIGNAL('stateChanged(int)'), self.handle_legend)

        self.polar = True
        self.polarCheck = QCheckBox("Polar")
        self.polarCheck.setChecked(True)
        self.connect(self.polarCheck, SIGNAL('stateChanged(int)'), self.handle_polar)

        self.unwrap_phase = False
        self.unwrapCheck = QCheckBox("Unwrap phase")
        self.connect(self.unwrapCheck, SIGNAL('stateChanged(int)'), self.handle_unwrap)

        self.reference = None
        self.reference_index = 0
        self.referenceSelector = QComboBox()
        self.referenceSelector.addItem("None")
        for parm in self.parms:
            if parm.isPolar():
                self.referenceSelector.addItem("%s (polar)" % parm.name())
            else:
                self.referenceSelector.addItem(parm.name())
        self.connect(self.referenceSelector, SIGNAL('activated(int)'), self.handle_reference)

        self.referenceLabel = QLabel("Phase reference")

        hbox = QHBoxLayout()
        hbox.addWidget(self.axisSelector)
        hbox.addWidget(self.spinner)
        hbox.addWidget(self.legendCheck)
        hbox.addWidget(self.polarCheck)
        hbox.addWidget(self.unwrapCheck)
        hbox.addWidget(self.referenceSelector)
        hbox.addWidget(self.referenceLabel)
        hbox.addStretch(1)

        layout = QVBoxLayout()
        layout.addLayout(hbox);
        layout.addWidget(self.canvas, 1)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)

        self.domain = common_domain(self.selected_parms)

        self.resolution = None
        if (not self.domain is None) and (not resolution is None):
            self.resolution = [min(max(resolution[0], 1.0), self.domain[1] - self.domain[0]),
                min(max(resolution[1], 1.0), self.domain[3] - self.domain[2])]

        self.shape = (1, 1)
        if not self.domain is None:
            self.shape = (self.selected_parms[0].value(self.domain, self.resolution)[0].shape)
            assert(len(self.shape) == 2)

        self.spinner.setRange(0, self.shape[1 - self.axis] - 1)
        self.plot()

    def plot(self):
        el0 = []
        el1 = []
        labels = []

        if not self.domain is None:
            for parm in self.selected_parms:
                value = parm.value(self.domain, self.resolution, self.polar, self.unwrap_phase, self.reference)
                if value[0].shape != self.shape or value[1].shape != self.shape:
                    print "warning: inconsistent shape; will skip parameter:", parm.name()
                    continue

                if self.axis == 0:
                    el0.append(value[0][:, self.index])
                    el1.append(value[1][:, self.index])
                else:
                    el0.append(value[0][self.index, :])
                    el1.append(value[1][self.index, :])
                labels.append(parm.name())

        legend = self.show_legend and len(labels) > 0
        xlabel = ["Time (sample)", "Freq (sample)"][self.axis]
        if self.polar:
            plot(self.fig, el0, sub="211", labels=labels, show_legend=legend, xlabel=xlabel, ylabel="Amplitude")
            plot(self.fig, el1, clf=False, sub="212", stack=True, scatter=True, labels=labels, show_legend=legend, xlabel=xlabel, ylabel="Phase (rad)")
        else:
            plot(self.fig, el0, sub="211", labels=labels, show_legend=legend, xlabel=xlabel, ylabel="Real")
            plot(self.fig, el1, clf=False, sub="212", labels=labels, show_legend=legend, xlabel=xlabel, ylabel="Imaginary")

        # Set x-axis scale in number of samples.
        for ax in self.fig.axes:
            ax.set_xlim(0, self.shape[self.axis] - 1)

        self.canvas.draw()

    def handle_spinner(self, index):
        self.index = index
        self.plot()

    def handle_axis(self, axis):
        if axis != self.axis:
            self.axis = axis
            self.spinner.setRange(0, self.shape[1 - self.axis] - 1)
            self.spinner.setValue(0)
            self.plot()

    def handle_legend(self, state):
        self.show_legend = (state == 2)
        self.plot()

    def handle_unwrap(self, state):
        self.unwrap_phase = (state == 2)
        self.plot()

    def handle_polar(self, state):
        self.polar = (state == 2)
        self.referenceSelector.setEnabled(self.polar)
        self.referenceLabel.setEnabled(self.polar)
        self.unwrapCheck.setEnabled(self.polar)
        self.plot()

    def handle_reference(self, index):
        if index != self.reference_index:
            if index == 0:
                reference = None
            else:
                parm = self.parms[index - 1]
                value = parm.value(self.domain, self.resolution)
                if value[1].shape != self.shape:
                    print "warning: inconsistent shape; will not change phase reference to parameter:", parm.name()
                    return
                reference = value[1]

            self.reference_index = index
            self.reference = reference
            self.plot()

class MainWindow(QFrame):
    def __init__(self, db):
        QFrame.__init__(self)
        self.db = db
        self.figures = []
        self.parms = []

        layout = QVBoxLayout()

        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.list, 1)

        self.useResolution = True
        checkResolution = QCheckBox("Use resolution")
        checkResolution.setChecked(True)
        self.connect(checkResolution, SIGNAL('stateChanged(int)'), self.handle_resolution)

        self.resolution = [QLineEdit(), QLineEdit()]
#        validator = QDoubleValidator(self.resolution[0])
#        validator.setRange(1.0, 2.0)
#        self.resolution[0].setValidator(validator)
        self.resolution[0].setAlignment(Qt.AlignRight)
        self.resolution[1].setAlignment(Qt.AlignRight)

        hbox = QHBoxLayout()
        hbox.addWidget(checkResolution)
        hbox.addWidget(self.resolution[0])
        hbox.addWidget(QLabel("Hz"))
        hbox.addWidget(self.resolution[1])
        hbox.addWidget(QLabel("s"))
        layout.addLayout(hbox)

        self.plot_button = QPushButton("Plot")
        self.connect(self.plot_button, SIGNAL('clicked()'), self.handle_plot)
        self.close_button = QPushButton("Close figures")
        self.connect(self.close_button, SIGNAL('clicked()'), self.handle_close)

        hbox = QHBoxLayout()
        hbox.addWidget(self.plot_button)
        hbox.addWidget(self.close_button)
        layout.addLayout(hbox)

        self.setLayout(layout)
        self.populate()

    def populate(self):
        names = self.db.getNames()
        while len(names) > 0:
            parm = names.pop()
            split = parm.split(":")

            if contains(split, "Real") or contains(split, "Imag"):
                if contains(split, "Real"):
                    idx = split.index("Real")
                    split[idx] = "Imag"
                    other = ":".join(split)
                    try:
                        names.pop(names.index(other))
                    except ValueError:
                        other = None
                    elements = [parm, other]
                else:
                    idx = split.index("Imag")
                    split[idx] = "Real"
                    other = ":".join(split)
                    try:
                        names.pop(names.index(other))
                    except ValueError:
                        other = None
                    elements = [other, parm]

                split.pop(idx)
                self.parms.append(Parm(self.db, ":".join(split), elements))

            elif contains(split, "Ampl") or contains(split, "Phase"):
                if contains(split, "Ampl"):
                    idx = split.index("Ampl")
                    split[idx] = "Phase"
                    other = ":".join(split)
                    try:
                        names.pop(names.index(other))
                    except ValueError:
                        other = None
                    elements = [parm, other]
                else:
                    idx = split.index("Phase")
                    split[idx] = "Ampl"
                    other = ":".join(split)
                    try:
                        names.pop(names.index(other))
                    except ValueError:
                        other = None
                    elements = [other, parm]

                split.pop(idx)
                self.parms.append(Parm(self.db, ":".join(split), elements, True))
            else:
                self.parms.append(Parm(self.db, parm))

        self.parms = [parm for parm in self.parms if not parm.empty()]
        self.parms.sort(cmp=lambda x, y: cmp(x.name(), y.name()))

        domain = common_domain(self.parms)
        if not domain is None:
            self.resolution[0].setText("%.6f" % ((domain[1] - domain[0]) / 100.0))
            self.resolution[1].setText("%.6f" % ((domain[3] - domain[2]) / 100.0))

        for parm in self.parms:
            name = parm.name()
            if parm.isPolar():
                name = "%s (polar)" % name

            QListWidgetItem(name, self.list)

    def close_all_figures(self):
        for figure in self.figures:
            figure.close()

        self.figures = []

    def handle_resolution(self, state):
        self.useResolution = (state == 2)

    def handle_plot(self):
        selection = [self.list.row(item) for item in self.list.selectedItems()]
        selection.sort()

        resolution = None
        if self.useResolution:
            resolution = [float(item.text()) for item in self.resolution]

        self.figures.append(PlotWindow(self.parms, selection, resolution, title="Figure %d" % (len(self.figures) + 1)))
        self.figures[-1].show()

    def handle_close(self):
        self.close_all_figures()

    def closeEvent(self, event):
        self.close_all_figures()
        event.accept()


if __name__ == "__main__":
    if len(sys.argv) <= 1 or sys.argv[1] == "--help":
        print "usage: parmdbplot.py <parmdb>"
        sys.exit(1)

    db = parmdb.parmdb(sys.argv[1])

    app = QApplication(sys.argv)
    window = MainWindow(db)
    window.show()

#    app.connect(app, SIGNAL('lastWindowClosed()'), app, SLOT('quit()'))
    app.exec_()
