# nicepy

Analysis package for the NICE experiment.

## Example Usage

```python
import nicepy as npy

params = {'delay': -3}

ts1 = TofSet(filenames, params, norm=False, fluor=True)

ts1.get_raw_means()

masses = {'Be': 9, 'BeOH': 26}
ts1.get_peaks(masses, pk_range=(-100, 100))
ts1.get_peak_means()

ts1.show_means()

ts1.show_peaks(norm=False, fmt=True, box_out=False, capsize=4, grid=True)
fig = ts1.peaks_fig_ax.loc['fig']
ax = ts1.peaks_fig_ax.loc['ax']
npy.format_ax(ax, box_out=False)
ax.set_xlabel('Delay (s)')
ax.set_ylabel('Ion Signal (arb)')
```
