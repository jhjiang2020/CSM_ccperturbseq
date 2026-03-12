import argparse
import json
import logging
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def fix_run_proc_defs(sublibrary_dirs):
    """
    Parses the run_proc_def.json from each sublibrary and harmonizes the 
    'sample' section so that all sublibraries contain all sample definitions.
    This fixes the split-pipe --mode comb 'list index out of range' error
    caused by non-overlapping sample names across sublibraries.
    """
    proc_defs = []
    all_samples_dict = {}
    
    # 1. Load and parse
    for sdir in sublibrary_dirs:
        pdir = Path(sdir) / "process"
        def_file = pdir / "run_proc_def.json"
        
        if not def_file.exists():
            logging.error(f"Cannot find run_proc_def.json in {pdir}")
            sys.exit(1)
            
        with open(def_file, 'r') as f:
            data = json.load(f)
            proc_defs.append((def_file, data))
            
        # Collect all unique samples based on their name
        if "sample" in data and "samples" in data["sample"]:
            for samp in data["sample"]["samples"]:
                sname = samp["name"]
                if sname not in all_samples_dict:
                    all_samples_dict[sname] = samp
                else:
                    # Prefer samples with more wells if duplicated across sublibs
                    if len(samp.get("well_list", [[]])[0]) > len(all_samples_dict[sname].get("well_list", [[]])[0]):
                         all_samples_dict[sname] = samp

    # We now have the universal superset of samples
    all_samples_list = list(all_samples_dict.values())
    num_samples = len(all_samples_list)
    
    logging.info(f"Gathered {num_samples} unique samples across sublibraries: {list(all_samples_dict.keys())}")

    # 2. Harmonize and write back
    for def_file, data in proc_defs:
        # Create a backup just in case
        backup_file = def_file.with_suffix('.json.bak')
        if not backup_file.exists():
            shutil.copy(def_file, backup_file)
            
        # Modify the data
        if "sample" not in data:
            data["sample"] = {}
        
        data["sample"]["number_samples"] = num_samples
        
        # When compiling the sample list for this specific sublibrary,
        # we still want to keep the dummy arrays length consistent.
        # But simply writing all_samples_list satisfies the combine mode check.
        # SpSamp.__init__ will check length across bc_rounds, and then look for the specific well list for this sublib
        # split-pipe spoutdir spoutdir_runproc just reads this file.
        # 
        # Actually split-pipe combine does the following:
        # for slib in self._sublibs:
        #     ssamp = slib.get_samples(one=name)
        #     if ssamp:
        #         ss_well_list = ssamp._all_well_list[0]
        #     else:
        #         ss_well_list = [] * 3   <-- BUG evaluates to [], causing [] instead of [[],[],[]]
        # 
        # We can bypass this by ensuring slib.get_samples(one=name) returns a valid dummy sample instead of None.
        
        sublib_sample_names = {s["name"] for s in data["sample"].get("samples", [])}
        
        # Add missing samples as dummy samples (empty well_list for 1 bc_round)
        # to ensure it returns an object and gets an empty well list but with CORRECT dimensions (e.g., [[]])
        new_samples = list(data["sample"].get("samples", []))
        
        for sname, sdef in all_samples_dict.items():
            if sname not in sublib_sample_names:
                dummy = sdef.copy()
                # Empty out the well lists so we don't count wells that don't belong to this sublib
                dummy["well_list"] = [[]]
                dummy["wspec"] = ""
                dummy["num_wells"] = [0]
                new_samples.append(dummy)
                
        data["sample"]["samples"] = new_samples
        data["sample"]["number_samples"] = len(new_samples)
        
        # Write it
        with open(def_file, 'w') as f:
            json.dump(data, f, indent=4)
            
        logging.info(f"Updated {def_file} with unified sample definitions.")

def main():
    parser = argparse.ArgumentParser(description="Fix split-pipe run_proc_def.json for combine mode")
    parser.add_argument("sublib_dirs", nargs='+', help="Output directories of sublibraries to combine")
    args = parser.parse_args()
    
    fix_run_proc_defs(args.sublib_dirs)

if __name__ == "__main__":
    main()
