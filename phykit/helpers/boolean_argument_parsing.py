import argparse

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('true', 'True', 't', 'T', '1'):
        return True
    elif v.lower() in ('false', 'False', 'f', 'F', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')