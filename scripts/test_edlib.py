import sys
import yaml

# Try importing edlib
try:
    import edlib
except ImportError:
    print("Error: edlib is not installed in this environment.")
    sys.exit(1)

config_path = '/medpop/rgupta/Jiahao/TIMI/CSM_ccperturbseq/qc_config.yaml'
try:
    with open(config_path) as f:
        config = yaml.safe_load(f)
except Exception as e:
    print(f"Error loading config at {config_path}: {e}")
    sys.exit(1)

stem1 = config['q3_constants']['stem_1']
stem2 = config['q3_constants']['stem_2']

print(f"Reference Stem 1: {stem1}")
print(f"Reference Stem 2: {stem2}")

# Calculate exact mismatch positions between the two stems
diffs = [i for i, (a, b) in enumerate(zip(stem1, stem2)) if a != b]
print(f"Stems differ at {len(diffs)} positions: {diffs}")

print("\n--- Testing Edlib HW Alignment Behavior ---")
# Simulate a read that contains exactly stem 1, but with sequencing errors
# Let's say it has 2 errors that DO NOT overlap the distinguishing bases
test_read = "A" * 15 + "GTT" + "A" + "GAGAGCT" + "C" + "AGCAGGA" + "C" * 15
print(f"Test Read (Stem 1 with 2 errors): ...{test_read[15:15+19]}...")

for max_mm in [2, 5]:
    s1_res = edlib.align(stem1, test_read, mode="HW", k=max_mm, task="path")
    s2_res = edlib.align(stem2, test_read, mode="HW", k=max_mm, task="path")
    
    # editDistance is -1 if not found within k
    d1 = s1_res['editDistance']
    d2 = s2_res['editDistance']
    
    print(f"\nMax Mismatches Allowed (k) = {max_mm}")
    print("  Matching against Stem 1:")
    if d1 >= 0:
        print(f"    FOUND (editDistance = {d1})")
        print(f"    CIGAR: {s1_res['cigar']}")
    else:
        print("    NOT FOUND")
        
    print("  Matching against Stem 2:")
    if d2 >= 0:
        print(f"    FOUND (editDistance = {d2})")
        print(f"    CIGAR: {s2_res['cigar']}")
        if max_mm == 5 and d1 == d2:
            print("    --> TIE! This read will be marked as AMBIGUOUS by the pipeline.")
    else:
        print("    NOT FOUND")

print("\n--- Why is this happening? ---")
print("Edlib 'HW' mode finds the substring in the read with the minimal edit distance.")
print("If max_mm is 5, a sequence that is 2 edits away from Stem 1 will ALSO be ")
print("at most 2 + 3 = 5 edits away from Stem 2. Because of insertion/deletion ")
print("alignments, it might even find a matching path of EXACTLY 2 edits away ")
print("by skipping the mismatched bases. Let's look at the CIGAR string for the ")
print("Stem 2 match to see if it's counting them as mismatches (X) or indels (I/D).")
