import glob

def find_last_modecar():
    max_val = -1
    print((glob.glob('./MODECAR*')))
    for modecar in glob.glob('./MODECAR*'):
        try:
            i = int(float(modecar[9:]))
        except ValueError:
            continue
        if i > max_val:
            max_val = i
    return max_val
