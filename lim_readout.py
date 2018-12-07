from Readout import Readout


grid = Readout()
grid.print_readout()
grid.centerline(0)
grid.avg_imp_vely(1)
grid.te_contour(2)
grid.deposition_contour(3, probe_width=0.015, rad_cutoff=0.05, side='ITF')
grid.deposition_contour(4, probe_width=0.015, rad_cutoff=0.05, side='OTF')
grid.ne_contour(5)
grid.show_fig()
