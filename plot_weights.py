#!/usr/bin/env python
""" Provides acces to the weights of an MS file. """
from __future__ import division

import matplotlib.colors as colors
import matplotlib.patheffects as pe
import numpy as np

from casacore.tables import taql
from matplotlib import cm
from matplotlib.pyplot import subplots
from matplotlib import ticker
from matplotlib.pyplot import show

__author__ = 'Frits Sweijen'
__credits__= 'Francesco de Gasperin'
__version__ = '1.0.0'
__maintainer__ = 'Frits Sweijen'
__email__ = 'sweijen <at> strw.leidenuniv.nl'
__status__ = 'Development'

def normalize(x, nmin, nmax):
    normed = (x - nmin) / (nmax - nmin)
    return normed

class WeightPlotter:
    def __init__(self, msfile):
        if msfile.endswith('/'):
            self.msfile = msfile[:-1]
        else:
            self.msfile = msfile

        self.colors = ['C1', 'C2', 'C3', 'C4']
        # Open the table and ignore rows that are completely flagged.
        print 'Loading '+msfile+'...'
        mstable_init = taql('SELECT TIME, ANTENNA1, FIELD_ID, DATA, WEIGHT_SPECTRUM, FLAG FROM $msfile WHERE !ALL(FLAG)')
        self.mstable = taql('SELECT TIME, MEANS(GAGGR(DATA), 0) AS DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHTS, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV, FLAG FROM $mstable_init GROUPBY TIME')

        # Additional processing before we can use the data.
        self.flags = self.mstable.getcol('FLAG')
        self.time = self.mstable.getcol('TIME')
        self.elevation = self.mstable.getcol('ELEV')

        data_um = self.mstable.getcol('DATA')
        self.data = np.ma.MaskedArray(data=data_um, mask=self.flags)

        weights_um = self.mstable.getcol('WEIGHTS')
        self.weights = np.ma.MaskedArray(data=weights_um, mask=self.flags)
        print 'FLAG table applied to DATA and WEIGHT_SPECTRUM.'
        # Obtain polarization setup.
        temp = taql('SELECT CORR_TYPE from '+msfile+'/POLARIZATION')
        if temp.getcol('CORR_TYPE')[0] in np.asarray([5, 6, 7, 8]):
            # Circular polarization.
            self.polarization = ['RR', 'LL', 'RL', 'LR']
        elif temp.getcol('CORR_TYPE')[0] in np.asarray([9, 10, 11, 12]):
            self.polarization = ['XX', 'YY', 'XY', 'YX']
        print 'Polarzation setup:', ', '.join(self.polarization)

        # Obtain channel frequencies.
        chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
        # Select the first table, column CHAN_FREQ and convert to MHz.
        self.freq = chan_freq[0]['CHAN_FREQ'] * 1e-6

        print self.msfile+' loaded.'

    def plot_data_2D(self):
        print 'Plotting data matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Normalized log10(1 + |DATA|) '+self.msfile, fontweight='bold')
        axes = axes.ravel()

        for (pol, pol_label) in enumerate(self.polarization):
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            #print np.min(np.mean(weights[:, :, pol], axis=0)), np.max(np.mean(weights[:, :, pol], axis=0))
            #weights = np.ma.masked_where(self.weights > 1e5, self.weights)
            data = np.log10(1 + np.abs(self.data))
            if pol == 0:
                nmin = np.nanmin(data[:, :, pol])
                nmax = np.nanmax(data[:, :, pol])
            nweights = normalize(data[:, :, pol], nmin, nmax)
            #print nmin, nmax
            #print nweights
            #im = axes[pol].imshow(nweights, extent=[freq[0], freq[-1], time[0], time[-1]])
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nweights, vmin=0.0, vmax=1.0)
            axes[pol].set_aspect('auto')
            axes[pol].label_outer()
            fig.colorbar(im, ax=axes[pol])
        imgname = 'data_2D.png'
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_variance_2D(self, delta=3):
        print 'Plotting variance matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Inverse Variance '+self.msfile+', box =%dx%d'%(delta, delta), fontweight='bold')
        axes = axes.ravel()
        for (pol, pol_label) in enumerate(self.polarization):
            print 'Processing polarization '+pol_label
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            # Calculate variance in the visibilities.
            var = np.abs(self.visibility_variance2D(self.data, self.weights, pol, delta=delta))
            var = np.ma.MaskedArray(data=var, mask=self.flags[:, :, pol])
            if pol == 0:
                nmin = np.nanmin(var)#; print nmin
                nmax = np.nanmax(var)#; print nmax
                #print np.median(var)
                #print np.percentile(var, 10)
                nmax = np.percentile(var, 99)
            nvar = normalize(var, nmin, nmax)
            #im = axes[pol].imshow(nvar, extent=[freq[0], freq[-1], time[0], time[-1]])
            #im = axes[pol].imshow(nvar, extent=[freq[0], freq[-1], time[0], time[-1]], norm=colors.LogNorm(vmin=1e-4, vmax=1.0), origin='lower')
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nvar, norm=colors.LogNorm(vmin=1e-4, vmax=1.0))
            axes[pol].set_aspect('auto')
            fig.colorbar(im, ax=axes[pol])
            #break
        imgname = 'variance_2D_box_%dx%d.png' % (delta, delta)
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_2D(self):
        print 'Plotting weight matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Normalized Weights '+self.msfile, fontweight='bold')
        axes = axes.ravel()

        for (pol, pol_label) in enumerate(self.polarization):
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            weights = self.weights[:, :, pol]
            # Deal with outliers in the weights to avoid plotting issues/biases.
            Q1 = np.percentile(weights, 25)
            Q3 = np.percentile(weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_max = (wmedian + 5*IQR)
            outlier_thresh_min = (wmedian - 5*IQR)
            if np.max(weights) > outlier_thresh_max or np.min(weights) < outlier_thresh_min:
                fig.suptitle('Normalized Weights '+self.msfile+' (IQR filtered)', fontweight='bold')
                weights = np.ma.masked_where(weights > outlier_thresh_max, weights)
                weights = np.ma.masked_where(weights < outlier_thresh_min, weights)
            if pol == 0:
                nmin = np.nanmin(weights[:, :])
                nmax = np.nanmax(weights[:, :])
            nweights = normalize(weights[:, :], nmin, nmax)
            #print nmin, nmax
            #print nweights
            #im = axes[pol].imshow(nweights, extent=[freq[0], freq[-1], time[0], time[-1]])
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nweights, vmin=0.0, vmax=1.0)
            axes[pol].set_aspect('auto')
            axes[pol].label_outer()
            fig.colorbar(im, ax=axes[pol])
        imgname = 'weight_2D.png'
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_frequency(self, delta=10):
        print 'Plotting weights vs. frequency...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Time-mean Weights ' + self.msfile + ' $\\Delta=%d$'%(delta,), fontweight='bold')
        axes = axes.ravel()
        handles = []
        labels = []

        for (pol, pol_label) in enumerate(self.polarization):
            data_shift = np.roll(self.data, -1, axis=1)
            data_sub = self.data - data_shift
            # Calculate the variance in frequency space.
            fdata = np.nanmean(data_sub[:, :, pol], axis=0)
            fvar = np.zeros(shape=fdata.shape[0] // delta)
            f_axis = np.empty(shape=fdata.shape[0] // delta)
            for i in xrange(fdata.shape[0] // delta):
                vr = np.nanvar(fdata.real[delta * i:delta * (i+1)])
                vi = np.nanvar(fdata.imag[delta * i:delta * (i+1)])
                v = (vr + vi) / 2.
                if v != 0 and np.isfinite(v):
                    fvar[i] = 1. / v
                f_axis[i] = np.mean(self.freq[delta * i:delta * (i+1)])
            nmin = np.nanmin(fvar)
            nmax = np.nanmax(fvar)
            nvar = normalize(fvar, nmin, nmax)
            # Deal with outliers in the weights to avoid plotting issues/biases.
            Q1 = np.percentile(self.weights, 25)
            Q3 = np.percentile(self.weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights) > outlier_thresh_high or np.min(self.weights) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights > outlier_thresh_high, self.weights)
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                axes[pol].set_title(pol_label+' (IQR filtered)')
            else:
                weights = self.weights
                axes[pol].set_title(pol_label)
            fweights = np.mean(weights, axis=0)
            # Normalize w.r.t. XX or RR.
            nmin = np.nanmin(fweights[:, 0])
            nmax = np.nanmax(fweights[:, 0])
            nweights = normalize(fweights[:, pol], nmin, nmax)
            # Take mean in time to plot weights as a function of frequency.
            #tweights = np.mean(weights[:, :, :], axis=0)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel=self.polarization[0]+' Normalized Weight')
            axes[pol].label_outer()
            axes[pol].plot(self.freq, nweights, color=self.colors[pol], label=pol_label+' MS Weights')
            axes[pol].plot(f_axis, nvar, '--d', color=self.colors[pol], label=pol_label+' Variance', linewidth=2, path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

            lhandle, llabel = axes[pol].get_legend_handles_labels()
            handles.append(lhandle)
            labels.append(llabel)
        leg = axes[3].legend(np.ravel(handles), np.ravel(labels), loc='upper center', bbox_to_anchor=(-0.1, -0.3), ncol=2, borderaxespad=0.0)
        imgname = 'weight_frequency_mean_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Time-median Weights ' + self.msfile + ', $\\Delta=%d$'%(delta,), fontweight='bold')
        axes = axes.ravel()
        handles = []
        labels = []
        for (pol, pol_label) in enumerate(self.polarization):
            # Deal with outliers.
            Q1 = np.percentile(self.weights, 25)
            Q3 = np.percentile(self.weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights) > outlier_thresh_high or np.min(self.weights) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights > outlier_thresh_high, self.weights)
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                fweights = np.median(weights, axis=0)
                axes[pol].set_title(pol_label+' (IQR filtered)')
            else:
                fweights = np.median(self.weights, axis=0)
                axes[pol].set_title(pol_label)

            # Normalize w.r.t. XX or RR.
            nmin = np.nanmin(fweights[:, 0])
            nmax = np.nanmax(fweights[:, 0])
            nweights = normalize(fweights[:, pol], nmin, nmax)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel=self.polarization[0]+' Normalized Weight')
            axes[pol].label_outer()
            axes[pol].plot(self.freq, nweights, color=self.colors[pol], label=pol_label+' MS Weights')
            axes[pol].plot(f_axis, nvar, '--d', color=self.colors[pol], label=pol_label+' Variance', linewidth=2, path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
            lhandle, llabel = axes[pol].get_legend_handles_labels()
            handles.append(lhandle)
            labels.append(llabel)
        leg = axes[3].legend(np.ravel(handles), np.ravel(labels), loc='upper center', bbox_to_anchor=(-0.1, -0.3), ncol=2, borderaxespad=0.0)
        imgname = 'weight_frequency_median_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_time(self, mode='mean', delta=100):
        print 'Plotting weights vs. time...'
        fig, axes = subplots(nrows=1, ncols=1)
        fig.suptitle('Weights ('+mode+') for '+self.msfile+', $\\Delta=100$', fontweight='bold')
        axes.set(xlabel='Time', ylabel=self.polarization[0]+' Normalized Weight')
        axes.label_outer()
        axes_elev = axes.twinx()
        axes_elev.set(ylabel='Elevation [deg]')

        # Subtract adjacent channels to remove any signal. This should leave us with noise.
        data_shift = np.roll(self.data, -1, axis=1)
        data_sub = self.data - data_shift

        # Take the mean or median over frequency to marginalize into the time domain.
        if mode == 'median':
            tdata = np.median(data_sub, axis=1)
        elif mode == 'mean':
            tdata = np.mean(data_sub, axis=1)
        else:
            print 'Unknown mode, using median.'
            tdata = np.median(data_sub, axis=1)
        variance = np.zeros(shape=(tdata.shape[0] // delta, tdata.shape[1]))

        tdatar = tdata.real
        tdatai = tdata.imag
        tdataf = np.abs(tdata)
        t_axis = []
        for i in xrange(tdata.shape[0] // delta):
            t = np.mean(self.time[delta * i: delta * (i+1)])
            t_axis.append(t)
            vr = np.nanvar(tdatar[delta * i: delta * (i+1), :], axis=0)
            vi = np.nanvar(tdatai[delta * i: delta * (i+1), :], axis=0)
            v = (vr + vi) / 2.
            variance[i] = np.where(np.isfinite(1. / v), 1. / v, np.nan)

        for (pol, pol_label) in enumerate(self.polarization):
            Q1 = np.percentile(self.weights[:, :, pol], 25)
            Q3 = np.percentile(self.weights[:, :, pol], 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights[:, :, pol])
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights[:, :, pol]) > outlier_thresh_high or np.min(self.weights[:, :, pol]) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights[:, :, pol] > outlier_thresh_high, self.weights[:, :, pol])
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                fig.suptitle('Weights ('+mode+') for '+self.msfile+', $\\Delta=100$ (IQR filtered)', fontweight='bold')
            else:
                weights = self.weights[:, :, pol]

            # Take median or mean in frequency to plot weights as a function of time.
            if mode == 'median':
                tweights = np.median(weights, axis=1)
            elif mode == 'mean':
                tweights = np.mean(weights, axis=1)
            else:
                print 'Unknown mode, using median.'
                tweights = np.median(weights, axis=1)
            #if mask_threshold is not None:
                # Deal with some extreme outliers that screw with the median in this case.
            tweights = np.ma.masked_where(tweights > 1e3, tweights)
            # Normalize w.r.t. XX or RR.
            if pol == 0:
                nmin = np.nanmin(tweights); nmax = np.nanmax(tweights)
            nweights = normalize(tweights, nmin, nmax)
            nvar = normalize(variance[:, pol], np.nanmin(variance[:, 0]), np.nanmax(variance[:, 0]))

            axes.plot(self.time, nweights, label=pol_label, alpha=0.5, color=self.colors[pol])
            axes.plot(t_axis, nvar, '--d', label='Variance '+pol_label, color=self.colors[pol])
        axes_elev.plot(self.time, self.elevation * 180/np.pi, '-k', linewidth=2, label='Elevation')
        handles, labels = axes.get_legend_handles_labels()
        handles2, labels2 = axes_elev.get_legend_handles_labels()
        leg = axes.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, borderaxespad=0.0)
        imgname = 'weight_'+mode+'_time_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

    def visibility_variance2D(self, vis, msweights, pol=None, delta=3):
        if delta % 2 == 0:
            raise Exception('Box size must be odd!')
        s = delta // 2
        if pol is not None:
            vis = vis[:, :, pol]
            weight = np.copy(msweights[:, :, pol])
            msweights = msweights[:, :, pol]
        max_x = vis.shape[1] - 1
        max_y = vis.shape[0] - 1

        for (y, x), _ in np.ndenumerate(vis):
            #if x == 0:
            #    print 'At row', y
            # Horizontal bounds.
            if (x - s) <= 0:
                bound_left = 0
            else:
                bound_left = (x - s)
            if (x + s) >= max_x:
                bound_right = max_x
            else:
                bound_right = (x + s) + 1
            # Vertical bounds.
            if (y - s) <= 0:
                bound_low = 0
            else:
                bound_low = (y - s)
            if (y + s) >= max_y:
                bound_high = max_y
            else:
                bound_high = (y + s) + 1
            box_vis = vis[bound_low:bound_high, bound_left:bound_right]
            if np.any(box_vis) and not np.all(box_vis.mask):
                v = np.nanvar(box_vis)
                if v != 0.0:
                    w = 1. / v
                else:
                    w = 0.
            else:
                w = 0.
            if np.isfinite(w):
                weight[y, x] = w
        return weight

if __name__ == '__main__':
    import sys
    # Get the MS filename.
    msfile = sys.argv[1]
    wp = WeightPlotter(msfile)
    wp.plot_weight_time(mode='mean')
    wp.plot_weight_time(mode='median')
    wp.plot_weight_frequency(delta=10)
    wp.plot_weight_2D()
    wp.plot_variance_2D(delta=3)
    #wp.plot_variance_2D(delta=7)
    #wp.plot_variance_2D(delta=9)
    wp.plot_data_2D()
    #from matplotlib.pyplot import show
    #show()
