# plot corners for mcmc

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from typing import Union, List, Dict, Any, Iterable

def scientificNotationTitle(v, *lu, digit=3):
	### show value, lower, and upper uncertainty in scientific notation format
	e = np.floor(np.log10(np.abs(v)))
	### dont convert xx.xx to x.xxx*10^1
	if e == 1:
		e = 0
		mind = digit-2 #keep dd.d if digit=3
	else: mind = digit-1 #keep d.dd if digit=3
	#mind is the digits after dot
	vs = v/10**e
	lus = [k/10**e for k in lu]
	### show more digit if any error is zero
	for d in range(mind, 11):
		lustr = [f'%+.{d}f' % k for k in lus]
		if all([float(k) != 0.0 for k in lustr]): break
	vstr = f'%+.{d}f' % vs
	lustr = [f'%+.{d}f' % k for k in lus]
	### print latex form
	vlatex = f'%+.{d-int(e)}f' % (v)
	lulatex = [f'%+.{d-int(e)}f' % (k) for k in lu]
	latex = '& %s' % vlatex
	for i in range(0, len(lulatex), 2):
		latex += ' & $^{%s}_{%s}$' % (lulatex[i+1], lulatex[i])
	print(latex) #export latex table
	if e==0:
		return '$%s_{%s}^{%s}$' % (vstr, *lustr[:2])
	else:
		return '$%s_{%s}^{%s}$ x$10^{%i}$' % (vstr, *lustr[:2], e)


def corners(steps, probs, bins=20, ranges=0.0, digits=3, order=None,
	fig=None, ax = None,
	color=None, alpha=None, linewidth=None, value_color=None, label=None,
	figure_kws = dict(figsize=(9,9)),
	hist_kws = {},
	labels = None,
	values = None,
	values_limit=[0.16, 0.84],
	hist_value_kws = dict(marker='v'),
	hist2d_value_kws = dict(marker='+'),
	hist_value_limit_kws = dict(color='k'),
	hist2d = True,
	hist2d_kws = dict(cmap='Blues'),
	hexbin = False,
	hexbin_kws = dict(cmap='Blues'),
	scatter = False,
	scatter_kws = dict(s=0.2, cmap='RdYlBu', alpha=0.3),
	contour = True,
	contour_smooth_sigma = 1,
	contour_kws = dict(levels=[0.2, 0.4, 0.6, 0.8]),
	legend_kws = dict(), 
	output = 'test.out',
	):
	'''
	steps:
		numpy.array of shape (nsteps x ndim)
	probs:
		numpy.array of shape (nsteps)
		weights of steps
	bins:
		int or list
		bins on hist and hist2d
	ranges:
		float or list
		determine range with step[probs >quantile(probs, float)]
	digits:
		int or list
		significant digits in title
	fig, ax:
		plt.Figure, array of plt.Axes
		figure and axes to overlay, created by previous corners of the same ndim
	color, alpha, linewidth:
		keywords for both hist and contour
	label:
		str
		handle for legend
	labels:
		str or list
		label name for ndim
	values:
		float or list
		value for marker on axes
	values_limit:
		[float, float]
		vertical lines for lower and upper limit
	figure_kws:
		dict
		keywords for figure
	hist_kws:
		dict
		keywords for hist
	hist_value_kws:
		dict
		keywords for markers on hist
		pos = 'bartop' / 'top'
	hist2d_value_kws:
		dict
		keywords for markers on hist2d
	hist_value_limit_kws:
		dict
		keywords for two limit lines on hist
	hist2d, hist2d_kws:
		bool, dict
		whether to plot hist2d and the corresponding keywords
		set histtype='bar' if need filled hist
	scatter, scatter_kws:
		bool, dict
		whether to plot scatter and the corresponding keywords
	contour, contour_kws:
		bool, dict
		whether to plot contour and the corresponding keywords
	contour_smooth_sigma:
		float
		width of kernel to smooth contour
	legend_kws:
		dict
		keywords for legend
	'''

	#default settings
	default_color = color if color else 'k'
	default_alpha = alpha if alpha else 1.0
	default_linewidth = linewidth if linewidth else 1.0
	default_values_color = value_color if value_color else 'tab:red'

	hist_kws = {**dict(histtype='step', density=True, color=default_color, linewidth=default_linewidth, alpha=default_alpha, zorder=100), **(hist_kws or {})} # hist at [i,i]
	hist_value_kws = {**dict(pos='bartop', marker='v', color=default_values_color, zorder=200), **(hist_value_kws or {})} # marker in hist
	hist2d_value_kws = {**dict(marker='+', color=default_values_color, markersize=8, markeredgewidth=2, markeredgecolor=default_values_color, zorder=200), **(hist2d_value_kws or {})} # marker in hist2d
	hist_value_limit_kws = {**dict(linestyle='--', color='k', linewidth=1, alpha=0.5, zorder=200), **(hist_value_limit_kws or {})} # vertical lines of upper/lower error

	hist2d_kws = {**dict(cmap='Blues', zorder=0), **(hist2d_kws or {})}
	hexbin_kws = {**dict(cmap='Blues', zorder=0), **(hexbin_kws or {})}
	scatter_kws = {**dict(s=0.2, cmap='RdYlBu', alpha=0.3, zorder=10), **(scatter_kws or {})}
	contour_kws = {**dict(colors=default_color, levels=[0.2, 0.4, 0.6, 0.8], linewidths=default_linewidth, alpha=default_alpha, zorder=100), **(contour_kws or {})}

	legend_kws = {**dict(loc='upper right', frameon=False)}

	nstep, ndim = steps.shape

	### sort steps (high prob steps are put in front)
	steps = steps[np.argsort(probs)]
	probs = np.sort(probs)

	### get limit of axes
	if ranges is None: ranges = [None] * ndim
	elif isinstance(ranges, (int, float)):
		ranges = [ranges] * ndim
	elif isinstance(ranges, Iterable):
		if len(ranges) != ndim:
			raise ValueError("Length of 'ranges' (%i) must match 'steps' dimensions (%i)" % (len(ranges), ndim))
	else:
		raise TypeError("'ranges' must be an int/float or a list/tuple")
	### get exact limit
	axlim = []
	for i,r in enumerate(ranges):
		if r is None:
			ran = (np.nanmin(steps[:,i]), np.nanmax(steps[:,i]))
		elif isinstance(r, (int, float)):
			idx = probs >= np.nanquantile(probs, r)
			ran = (np.nanmin(steps[idx,i]), np.nanmax(steps[idx,i]))
		else:
			ran = r[:2]
		axlim.append(ran)

	### get bins of axes
	if bins is None: bins = [10] * ndim
	elif isinstance(bins, int):
		bins = [bins] * ndim
	elif isinstance(bins, Iterable):
		if len(bins) != ndim:
			raise ValueError("Length of 'bins' (%i) must match 'steps' dimensions (%i)" % (len(bins), ndim))
	else:
		raise TypeError("'bins' must be an int/float or a list/tuple")

	### get labels of axes title
	if labels is None: draw_labels = False
	else:
		draw_labels = True
		if isinstance(labels, str):
			labels = [labels] * ndim
		elif isinstance(labels, Iterable):
			if len(labels) != ndim:
				raise ValueError("Length of 'labels' (%i) must match 'steps' dimensions (%i)" % (len(labels), ndim))
		else:
			raise TypeError("'labels' must be an int/float or a list/tuple")

	### get value of axes
	if values is None: draw_values = False
	else:
		draw_values = True
		if isinstance(values, (int, float)):
			values = [values] * ndim
		elif isinstance(values, Iterable):
			if len(values) != ndim:
				raise ValueError("Length of 'values' (%i) must match 'steps' dimensions (%i)" % (len(values), ndim))
		else:
			raise TypeError("'values' must be an int/float or a list/tuple")

	### get digits of axes title
	if digits is None: digits = [3] * ndim
	elif isinstance(digits, int):
		digits = [digits] * ndim
	elif isinstance(digits, Iterable):
		if len(digits) != ndim:
			raise ValueError("Length of 'digits' (%i) must match 'steps' dimensions (%i)" % (len(digits), ndim))
	else:
		raise TypeError("'digits' must be an int/float or a list/tuple")

	if order:
		steps = steps[:, order]
		axlim = [axlim[o] for o in order] # ranges
		bins = [bins[o] for o in order]
		if draw_labels: labels = [labels[o] for o in order]
		if draw_values:
			values = [values[o] for o in order]
			digits = [digits[o] for o in order]

	### create figure and axes.
	newfig = (fig is None) & (ax is None)
	if newfig:
		fig, ax = plt.subplots(nrows=ndim, ncols=ndim, **figure_kws)
		plt.subplots_adjust(bottom=0.08, top=0.91, left=0.08, right=0.91, wspace=0.04, hspace=0.04)
		fig.align_labels(ax)
	else:
		for i in range(ndim):
			lim = ax[i,i].get_xlim()
			axlim[i] = [min(*lim, *axlim[i]), max(*lim, *axlim[i])]

	### decide where to put marker in hist, default is on top of hist bar
	valuePos = hist_value_kws.pop('pos', 'bartop')
	### contour levels in percent
	contourLevels = np.array(contour_kws.pop('levels', [0.2, 0.4, 0.6, 0.8]))

	### output values and limits
	legend_handles = []
	for i in range(ndim):
		### hist at [i,i]
		h,xe,_ = ax[i,i].hist(steps[:,i], bins=np.linspace(*axlim[i], bins[i]), label=label, **hist_kws)
		### draw marker/upperlim/lowerlim in hist
		if draw_values:
			# value markers
			if valuePos == 'bartop':
				upperY = h[np.searchsorted(xe, values[i], side='right')-1] * 1.1
			elif valuePos == 'top':
				upperY = np.max(h)
			ax[i,i].plot(values[i], upperY, **hist_value_kws)

			### upper/lower limits
			lim = np.quantile(steps[:,i], values_limit)
			ax[i,i].axvline(lim[0], **hist_value_limit_kws)
			ax[i,i].axvline(lim[1], **hist_value_limit_kws)

			titleValue = scientificNotationTitle(values[i], lim[0]-values[i], lim[1]-values[i], digit=digits[i])
		else:
			titleValue = ''

		### labels
		if draw_values:
			if labels: title = '%s=%s' % (labels[i], titleValue)
			else: title = titleValue
			#ax[i,i].text(0, 1.02, title, transform=ax[i,i].transAxes, ha='left', va='bottom', **label_kws)
			ax[i,i].set_title(title, loc='left')

		### hide x ticklabels except the bottom row
		if i<ndim-1: ax[i,i].xaxis.set_ticklabels([])
		elif draw_labels:
			### show x labels for the bottom riow
			#ax[i,i].xaxis.set_label_coords(0.5, 0)
			ax[i,i].set_xlabel(labels[i])
		### hide all y tickslables
		ax[i,i].set_yticklabels([])

		### ticks
		ax[i,i].set_autoscale_on(True)
		ax[i,i].minorticks_on()
		ax[i,i].tick_params(top=True, left=False, direction='in', labelrotation=45, length=5)
		ax[i,i].tick_params(which='minor', top=True, left=False, direction='in', length=2)
		ax[i,i].set_xlim(axlim[i])


		### plot hist2d
		for j in range(i):
			h, xe, ye = np.histogram2d(steps[:,j], steps[:,i], bins=(np.linspace(*axlim[j], bins[j]), np.linspace(*axlim[i], bins[i])))
			### hist2d
			if hist2d:
				ax[i,j].imshow(h.T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect='auto', **hist2d_kws)
			if hexbin:
				if bins[j] == bins[i]: gs = bins[i]
				else: gs = (bins[j], int(bins[i]*0.6))
				ax[i,j].hexbin(steps[:,j], steps[:,i], gridsize=gs, extent=list(axlim[j])+list(axlim[i]), **hexbin_kws)
			### scatter
			if scatter:
				ax[i,j].scatter(steps[:,j], steps[:,i], c=probs, **scatter_kws)
			### contour
			if contour:
				if contour_smooth_sigma is not None:
					h = gaussian_filter(h, sigma=contour_smooth_sigma)
				xc = (xe[1:]+xe[:-1])/2
				yc = (ye[1:]+ye[:-1])/2
				ax[i,j].contour(xc, yc, h.T, levels=contourLevels*np.max(h), **contour_kws)

			### draw marker in hist2d
			if draw_values:
				ax[i,j].plot(values[j], values[i], **hist2d_value_kws)
			
			if i<ndim-1: 
				### hide x ticklabels except the bottom row
				ax[i,j].xaxis.set_ticklabels([])
			elif draw_labels:
				### show x labels for the bottom row
				#ax[i,j].xaxis.set_label_coords(0.5, 0)
				ax[i,j].set_xlabel(labels[j])
			if j>0:
				### hide y ticklabels except the left column
				ax[i,j].yaxis.set_ticklabels([])
			elif draw_labels:
				### show x labels for the bottom row
				#ax[i,j].yaxis.set_label_coords(0, 0.5)
				ax[i,j].set_ylabel(labels[i])

			### ticks
			ax[i,i].set_autoscale_on(False)
			ax[i,j].minorticks_on()
			ax[i,j].tick_params(top=True, right=True, direction='in', labelrotation=45, length=5)
			ax[i,j].tick_params(which='minor', top=True, right=True, direction='in', length=2)
			ax[i,j].set_xlim(axlim[j])
			ax[i,j].set_ylim(axlim[i])

			### hide upper right
			ax[j,i].axes.set_axis_off()


			'''
			legend_handles.append = [
				Line2D([0], [0], color='red', lw=2, label='Red Line'),
				Line2D([0], [0], color='blue', lw=2, label='Blue Line'),
				Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='Green Dot')
			]
			'''
	if label:
		ax[0,ndim-1].plot([0,1], color=hist_kws.pop('color'), alpha=hist_kws.pop('alpha'), linewidth=hist_kws.pop('linewidth'), label=label)
		ax[0,ndim-1].set_xlim(2,3)
		ax[0,ndim-1].legend(**legend_kws)
	return fig, ax


if __name__ == '__main__':
	plt.rcParams['axes.titlesize'] = 13
	plt.rcParams['axes.linewidth'] = 1
	plt.rcParams['axes.labelpad'] = 0
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['axes.labelweight'] = 'bold'
	plt.rcParams['xtick.labelsize'] = 12
	plt.rcParams['ytick.labelsize'] = 12
	plt.rcParams['legend.fontsize'] = 13

	if 0:
		kws = dict(hist2d=0, bins=30, ranges=0.2, linewidth=1.2, digits=[3,3,2,2,3,3,3])
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_100pc_sqrtmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_100pc_sqrtmass.npy')
		fig, ax = corners(steps, probs, color='#005500', alpha=0.5, label='w=$\sqrt{m}$ H=100pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_150pc_sqrtmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_150pc_sqrtmass.npy')
		best = steps[np.argmax(probs)]
		#labs = ['$a_1$', '$R_{w1}$', '$b_{w1}$', '$\phi_{w1}$']
		labs = ['$a_1$', '$R_{w1}$', '$\phi_{w1}$', '$a_2$', '$R_{w2}$', '$b_{w2}$', '$\phi_{w2}$']
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='k', alpha=0.5, values=best, labels=labs, label='w=$\sqrt{m}$ H=150pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_300pc_sqrtmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_300pc_sqrtmass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='green', alpha=0.2, label='w=$\sqrt{m}$ H=300pc', **kws)
		
		
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_100pc_mass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_100pc_mass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='#550000', alpha=0.5, label='w=m H=100pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_150pc_mass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_150pc_mass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='red', alpha=0.5, label='w=m H=150pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_300pc_mass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_300pc_mass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='red', alpha=0.2, label='w=m H=300pc', **kws)

		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_100pc_lgmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_100pc_lgmass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='#000055', alpha=0.5, label='w=$lg$m H=100pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_150pc_lgmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_150pc_lgmass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='blue', alpha=0.5, label='w=$lg$m H=150pc', **kws)
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/steps_errD_300pc_lgmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/twoCompExc/probs_errD_300pc_lgmass.npy')
		fig, ax = corners(steps, probs, fig=fig, ax=ax, color='blue', alpha=0.2, label='w=$lg$m H=300pc', **kws)

		fig.savefig('/Users/shaobo/Work/colleague/ysun/vertical/fig/corner_comp2_all.pdf', bbox_inches='tight')
	else:
		steps = np.load('/Users/shaobo/Work/colleague/ysun/vertical/oneCompExc/steps_errD_150pc_sqrtmass.npy')
		probs = np.load('/Users/shaobo/Work/colleague/ysun/vertical/oneCompExc/probs_errD_150pc_sqrtmass.npy')
		best = steps[np.argmax(probs)]
		labs = ['$a_1$', '$R_{w1}$', '$b_{w1}$', '$\phi_{w1}$']
		#labs = ['$a_1$', '$R_{w1}$', '$\phi_{w1}$', '$a_2$', '$R_{w2}$', '$b_{w2}$', '$\phi_{w2}$']
		fig, ax = corners(steps, probs, color='k', alpha=0.5, values=best, labels=labs, bins=21, ranges=0.1, linewidth=1.2, hist2d=0, hexbin=1, contour_smooth_sigma=0.5, hist_kws={'histtype':'bar', 'color':'tab:blue', 'alpha':0.5})

	plt.show()