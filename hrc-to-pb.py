import argparse
import sys

import hrc_instance
from pb_model import PBModel
            
def main(lines, max_bp, quiet, flatzinc):
    instance = hrc_instance.Instance(lines, PBModel(flatzinc), max_bp)
    instance.write(quiet)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Translate a MIN BP HRC instance to .opb format")
    parser.add_argument("max_bp", type=int,
            help="The maximum permitted number of blocking pairs")
    parser.add_argument("--quiet", "-q", action="store_true", required=False,
            help="Suppress most comments in output")
    parser.add_argument("--flatzinc", "-f", action="store_true", required=False,
            help="Output FlatZinc")
    args = parser.parse_args()
    
    if args.flatzinc:
        print "Error: flatzinc output not yet implemented"
        sys.exit(1)

    main([line.strip() for line in sys.stdin.readlines() if line.strip()],
            args.max_bp, args.quiet, args.flatzinc)
