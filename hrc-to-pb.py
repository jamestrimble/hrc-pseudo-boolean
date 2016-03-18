import argparse
import sys

import hrc_instance
from pb_model import PBModel
            
def main(lines, max_bp, quiet, flatzinc):
    instance = hrc_instance.Instance(lines, PBModel(flatzinc), max_bp)
    instance.write(quiet)

def show_sol(lines, sol_filename):
    instance = hrc_instance.Instance(lines, PBModel(False), 0, False)
    instance.show_sol(sol_filename)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Translate a MIN BP HRC instance to .opb format")
    parser.add_argument("max_bp", type=int,
            help="The maximum permitted number of blocking pairs")
    parser.add_argument("--quiet", "-q", action="store_true", required=False,
            help="Suppress most comments in output")
    parser.add_argument("--flatzinc", "-f", action="store_true", required=False,
            help="Output FlatZinc")
    parser.add_argument("--show-sol", type=str, required=False,
            help="Show a solution from file")
    args = parser.parse_args()
    
    if args.show_sol:
        show_sol([line.strip() for line in sys.stdin.readlines() if line.strip()],
                args.show_sol)
    else:
        main([line.strip() for line in sys.stdin.readlines() if line.strip()],
                args.max_bp, args.quiet, args.flatzinc)
