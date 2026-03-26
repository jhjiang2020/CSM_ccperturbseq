import sys
import os
import runpy

# Intercept wrapper-specific arguments
wrapper_config = None
splitpipe_dir = None

new_argv = [sys.argv[0]]
i = 1
while i < len(sys.argv):
    arg = sys.argv[i]
    if arg == '--wrapper_config':
        wrapper_config = sys.argv[i+1]
        i += 2
    elif arg == '--wrapper_splitpipe_dir':
        splitpipe_dir = sys.argv[i+1]
        i += 2
    else:
        new_argv.append(arg)
        i += 1
sys.argv = new_argv

if not splitpipe_dir or not wrapper_config:
    print("ERROR: splitpipe_wrapper.py requires --wrapper_config and --wrapper_splitpipe_dir")
    sys.exit(1)

sys.path.insert(0, splitpipe_dir)

import json
import yaml
import argparse

# --- Pre-execution Parse Sublibrary Metadata Patching ---
# To skip re-running step1a, we just update the sample names in the existing run_proc_def.json
# files to match what's in config.yaml. This tricks combine into working seamlessly.

def patch_run_proc_defs(argv):
    # Quick parse to find sublibraries passed to split-pipe
    sublibs = []
    try:
        idx = argv.index('--sublibraries')
        for arg in argv[idx+1:]:
            if arg.startswith('-'): break
            sublibs.append(arg)
    except ValueError:
        pass
        
    if not os.path.exists(wrapper_config) or not sublibs:
        return
        
    with open(wrapper_config, 'r') as f:
        cfg = yaml.safe_load(f)
        
    lib_map = {lib['sublibrary_name']: lib['sample_name'] for lib in cfg.get('libraries', [])}
    
    for sdir in sublibs:
        # Extract C1, C2 etc from something like "output/gex_C1"
        sublib_name = os.path.basename(sdir).replace('gex_', '')
        
        def_file = os.path.join(sdir, "process", "run_proc_def.json")
        if os.path.exists(def_file) and sublib_name in lib_map:
            target_sample = lib_map[sublib_name]
            
            with open(def_file, 'r') as f:
                data = json.load(f)
                
            changed = False
            # Safely navigate sample definitions
            if "sample" in data and "samples" in data["sample"]:
                for samp in data["sample"]["samples"]:
                    old_name = samp.get("name")
                    # We want to change the sample name to the new merged one
                    if old_name and old_name != target_sample:
                        samp["name"] = target_sample
                        changed = True

            # ALWAYS verify the directory is named correctly, even if JSON was already patched
            new_dir = os.path.join(sdir, target_sample)
            if not os.path.exists(new_dir):
                renamed = False
                for item in os.listdir(sdir):
                    potential_dir = os.path.join(sdir, item)
                    if os.path.isdir(potential_dir) and item != target_sample:
                        if os.path.exists(os.path.join(potential_dir, "DGE_unfiltered")):
                            os.rename(potential_dir, new_dir)
                            print(f"Wrapper: Robustly renamed directory {potential_dir} -> {new_dir}")
                            renamed = True
                            break
                if not renamed:
                    print(f"Wrapper Warning: Could not find DGE folders to rename to {target_sample} in {sdir}")

            if changed:
                with open(def_file, 'w') as f:
                    json.dump(data, f, indent=4)
                print(f"Wrapper: Patched {def_file} to use sample '{target_sample}'")

# Run the patch function on the command line arguments
patch_run_proc_defs(sys.argv)

try:
    from splitpipe import spsamp
except ModuleNotFoundError:
    print("WARNING: Could not import splitpipe. Running in environment without Parse Pipeline?")
    spsamp = None

if spsamp:
    orig_init = spsamp.SpSamp.__init__

    def patched_init(self, *args, **kwargs):
        # Call the original __init__
        orig_init(self, *args, **kwargs)
        
        # In combine mode, SpSamp iterates through sublibraries:
        # for slib in self._sublibs:
        #     ssamp = slib.get_samples(one=name)
        #     if ssamp:
        #         ss_well_list = ssamp._all_well_list[0]
        #     else:
        #         ss_well_list = [] * 3   <-- BUG evaluates to [], causing [] instead of [[],[],[]]
        # 
        # Here we fix `_all_well_list` to maintain the expected dimensionality.
        
        if getattr(self, '_sublibs', None):
            # Find out the correct number of bc rounds by looking at the first valid sub-sample
            # or defaulting to 3 (which split-pipe originally intended with `[] * 3`).
            # actually, SpSamp init defines: self._wplates (list of WellPlates for each round)
            # but in combine mode it might not be fully populated until later.
            
            # The safest approach is to ensure every element in _all_well_list is a list of lists 
            # if it's currently empty.
            # When `ss_well_list` becomes `[]`, `_all_well_list.append(ss_well_list)` adds `[]`.
            for i in range(len(self._all_well_list)):
                if not self._all_well_list[i]:
                    # split-pipe expects `self._all_well_list[sublib][bc_round]`
                    # We give it enough empty arrays to avoid IndexError for up to 3 rounds.
                    self._all_well_list[i] = [[], [], []]

    # Apply monkeypatch
    spsamp.SpSamp.__init__ = patched_init

    # Monkeypatch comb_max_sublibs limit
    try:
        from splitpipe import spclass
        orig_get_par_val = spclass.SplitPipe.get_par_val

        def patched_get_par_val(self, par, *args, **kwargs):
            if par == "comb_max_sublibs":
                return 1000
            return orig_get_par_val(self, par, *args, **kwargs)

        spclass.SplitPipe.get_par_val = patched_get_par_val
    except ModuleNotFoundError:
        print("WARNING: Could not import spclass for monkeypatching.")

# Execute the original split-pipe but with the patched method in memory
splitpipe_bin = os.path.join(splitpipe_dir, "split-pipe")

# Replace this wrapper script name with the original split-pipe program name in args
if len(sys.argv) > 0:
    sys.argv[0] = splitpipe_bin

try:
    runpy.run_path(splitpipe_bin, run_name="__main__")
except SystemExit as e:
    # split-pipe usually throws SystemExit(0) on success.
    # However, because it relies on multiprocessing, it often leaves dangling worker pools
    # that prevent Python from exiting cleanly, causing the wrapper to hang indefinitely.
    # We force a hard exit here to terminate all threads.
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        code = int(str(e.code))
    except (TypeError, ValueError):
        code = 0
    os._exit(code)
except Exception as e:
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(1)

sys.stdout.flush()
sys.stderr.flush()
os._exit(0)
