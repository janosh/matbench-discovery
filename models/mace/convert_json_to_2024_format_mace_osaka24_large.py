import json
import gzip
import os

# --- Execution Settings ---
# Specify input filenames here (supports both .gz and .json)
ref_filename = '../../data/phonons/2024-11-09-kappas-phononDB-PBE-noNAC.json.gz'  # Reference file
model_name = 'mace-osaka24-large'
json_name = 'force-sets.json.gz'
target_filename = f'./{model_name}/2025-11-21-kappa-103-FIRE-dist=0.01-fmax=0.0001-symprec=1e-05/{json_name}' # Target file
output_filename = f'./{model_name}/2025-11-21-kappa-103-FIRE-dist=0.01-fmax=0.0001-symprec=1e-05/reformated_{json_name}'

def load_json_smart(filepath):
    """
    Function to load JSON files.
    Automatically handles gzip decompression if the file ends with .gz.
    """
    print(f"Loading: {filepath}")
    try:
        if filepath.endswith('.gz'):
            # Open with gzip module in text mode (rt) if extension is .gz
            with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                return json.load(f)
        else:
            # Open normally for standard .json files
            with open(filepath, 'r', encoding='utf-8') as f:
                return json.load(f)
    except FileNotFoundError:
        print(f"Error: File {filepath} not found.")
        return None

def align_sort_order_gz_input(ref_file, target_file, output_file):
    print("Starting process...")

    # 1. Load the reference file (2024 version)
    # We read this to establish the "correct sort order" of mp_ids.
    ref_data = load_json_smart(ref_file)
    if ref_data is None: return

    # 2. Load the target file (2025 version)
    target_data = load_json_smart(target_file)
    if target_data is None: return

    # --- Establish Sorting Criteria ---
    # Retrieve the mp_id dictionary from the reference file.
    ref_mp_ids_dict = ref_data.get('mp_id', {})
    
    # Create a list of mp_ids sorted by their index (0, 1, 2...) to ensure correct order.
    sorted_indices = sorted(ref_mp_ids_dict.keys(), key=lambda x: int(x))
    ordered_mp_ids = [ref_mp_ids_dict[idx] for idx in sorted_indices]
    
    print(f"Number of materials in reference file: {len(ordered_mp_ids)}")

    # --- Prepare Target Data Lookup ---
    # Create a reverse lookup dictionary for the target file: {material_id: index}
    target_mp_ids_dict = target_data.get('material_id', {})
    mpid_to_target_idx = {v: k for k, v in target_mp_ids_dict.items()}

    # --- Construct Output Data ---
    data_output = {}
    
    # List of keys to keep in the output file
    keys_to_keep = [
        'ph_freqs', 'q_points', 
        'kappa_tot_rta', 'mode_kappa_tot_rta', 'kappa_p_rta', 'kappa_c', 'mode_weights',
        'max_stress', 'reached_max_steps', 'broken_symmetry', 'has_imag_ph_modes'
    ]

    # Initialize dictionary for mp_id
    data_output['mp_id'] = {}
    # Initialize dictionaries for other keys if they exist in the target data
    for key in keys_to_keep:
        if key in target_data:
            data_output[key] = {}

    print("Aligning and sorting data...")
    missing_count = 0

    # Iterate through the reference order (0, 1, 2...) and populate the new data
    for new_idx_int, mpid in enumerate(ordered_mp_ids):
        new_idx_str = str(new_idx_int) # New index "0", "1"...

        # Check if this mp_id exists in the target file (2025)
        if mpid in mpid_to_target_idx:
            old_idx_str = mpid_to_target_idx[mpid] # Original index in target file

            # 1. Save mp_id
            data_output['mp_id'][new_idx_str] = mpid

            # 2. Copy other data fields
            for key in keys_to_keep:
                if key in target_data and old_idx_str in target_data[key]:
                    data_output[key][new_idx_str] = target_data[key][old_idx_str]
        else:
            missing_count += 1
            # Skip if ID is not found (missing data)

    # --- Report Results ---
    if missing_count > 0:
        print(f"Warning: {missing_count} materials were not found in the target file.")
    else:
        print("All IDs matched successfully.")

    # --- Save Output (.gz) ---
    if not output_file.endswith('.gz'):
        output_file += '.gz'

    print(f"Saving to: {output_file}")
    with gzip.open(output_file, 'wt', encoding='utf-8') as f:
        json.dump(data_output, f, indent=None, separators=(',', ':'))
    
    print("Completed. The file is now aligned with the 2024 version.")

# Execution block
if __name__ == "__main__":
    align_sort_order_gz_input(ref_filename, target_filename, output_filename)
    pass
