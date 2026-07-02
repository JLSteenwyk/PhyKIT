def str2bool(v):
    if isinstance(v, bool):
        return v
    value = v.lower()
    if value in ("true", "t", "1"):
        return True
    elif value in ("false", "f", "0"):
        return False
    else:
        import argparse

        raise argparse.ArgumentTypeError("Boolean value expected.")
