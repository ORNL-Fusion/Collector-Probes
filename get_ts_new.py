from ThomsonClass import ThomsonClass


def ts_2d(shot, system, times, ref_time, average_ts=5, detach_front_te=5):
    ts = ThomsonClass(shot, system)
    ts.load_ts()
    ts.map_to_efit(times, ref_time, average_ts)
    ts.heatmap(ref_time, detach_front_te=detach_front_te)

# Note: average_ts averages the neighboring amount of points together, where
# each point is separated by about 20 ms.
