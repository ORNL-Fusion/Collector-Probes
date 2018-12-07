from ThomsonClass import ThomsonClass


def ts_2d(shot, system, times, ref_time, average_ts=5):
    ts = ThomsonClass(shot, system)
    ts.load_ts()
    ts.map_to_efit(times, ref_time, average_ts)
    ts.heatmap(ref_time, te_lims=(0, 50))
