import argparse
import sys

def str2bool(v):
    if v is None:
        return None
    if isinstance(v, bool):
        return v
    if str(v).lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif str(v).lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

parser = argparse.ArgumentParser()
parser.add_argument("--test", type=str2bool, nargs='?', const=True, default=None)

try:
    args = parser.parse_args()
    print(f"RES: {args.test}")
except SystemExit:
    print("EXIT")
except Exception as e:
    print(f"ERR: {e}")
