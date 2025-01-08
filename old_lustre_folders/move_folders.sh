#!/bin/bash

# Check if the user provided arguments
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <username> <start_index> <end_index>"
  echo "  <start_index>: Starting index of the chunk (0-based)"
  echo "  <end_index>: Ending index of the chunk (inclusive)"
  exit 1
fi

# Set the username and chunk indices
username="$1"
start_index="$2"
end_index="$3"

# Set the source and destination paths
source_path="/old_lustre/home/$username"
dest_path="/home/$username/Scratch"

# Check if the source path exists
if [ ! -d "$source_path" ]; then
  echo "Error: Source path '$source_path' does not exist."
  exit 1
fi

# Generate the folder list only once and store it in a file
folder_list_file="/home/$username/Scratch/folder_list_$username.txt" 
if [ ! -f "$folder_list_file" ]; then
  find "$source_path" -maxdepth 1 -type d -print0 > "$folder_list_file"
fi

# Read the folder list from the file into an array
folders=()
{ 
  read -r -d '' && while IFS= read -r -d '' line; do
    folders+=("$line")
  done
} < "$folder_list_file"

# Check if the indices are valid
if (( start_index < 0 || start_index >= ${#folders[@]} || end_index < start_index || end_index >= ${#folders[@]} )); then
  echo "Error: Invalid chunk indices."
  exit 1
fi

# Loop through the specified chunk of folders
for i in $(seq "$start_index" "$end_index"); do
  folder="${folders[$i]}"
  folder_name=$(basename "$folder")

  echo "Transferring: $folder_name"
  # Create the destination folder if it doesn't exist
  #mkdir -p "$dest_path/$folder_name"

  # Copy the folder contents to the destination folder
  cp -r "$folder/." "$dest_path/$folder_name"

  # Echo the folder that was transferred
  echo "Transferred: $folder_name"
done

echo "Finished copying chunk of folders for user '$username'."
