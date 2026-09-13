def str2bool(v):
    """Parse case-insensitive boolean strings, passing bool objects through.

    Accept true/t/1/yes/y as True and false/f/0/no/n as False.
    Invalid strings raise argparse.ArgumentTypeError.
    """
    if isinstance(v, bool):
        return v
    value = v.lower()
    if value in ("true", "t", "1", "yes", "y"):
        return True
    elif value in ("false", "f", "0", "no", "n"):
        return False
    else:
        import argparse

        raise argparse.ArgumentTypeError("Boolean value expected.")
