## Tutorial: Transfer Folders in Chunks

This script allows you to transfer folders from a source directory to a destination directory in chunks, giving you fine-grained control over the process.

### Prerequisites

- A Bash environment.

### Instructions

1. **Save the script:** Save the code above as a `.sh` file (e.g., `transfer_folders.sh`).

2. **Make it executable:**
   ```bash
   chmod +x transfer_folders.sh
   ```

3. **Run the script:**
   ```bash
   ./transfer_folders.sh <username> <start_index> <end_index>
   ```
   - Replace `<username>` with the actual username.
   - Replace `<start_index>` with the starting index of the chunk (0-based).
   - Replace `<end_index>` with the ending index of the chunk (inclusive).

### Example
   ```bash
   ./transfer_folders.sh mmm0666 0 9   # Transfers the first 10 folders
   ./transfer_folders.sh mmm0666 10 19  # Transfers the next 10 folders
   ```

### Explanation
- The script first generates a list of all folders in the source directory (`/old_lustre/home/<username>`) and stores it in a file (`/home/<username>/Scratch/folder_list_<username>.txt`). This list is generated only once.
- You then specify the chunk of folders to transfer using the `start_index` and `end_index` arguments.
- The script copies the specified folders from the source directory to the destination directory (`/home/<username>/Scratch`).

### Important Notes
- The indices are 0-based, so the first folder has an index of 0.
- You can modify the `source_path` and `dest_path` variables in the script if needed.
- The `cp -r` command copies folders recursively (including their contents). If you only want to copy the folders themselves, remove the `-r` option.
