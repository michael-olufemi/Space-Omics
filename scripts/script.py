import requests
import os

def download_fastq_files_for_study(study_id, version="1", out_dir="downloads"):
    """
    Downloads all .fastq.gz files for a given study ID from NASA OSDR.
    
    Parameters:
      study_id (str): The numeric part of the study ID (e.g. "667" for OSD-667).
      version (str): Study version (default "1" for version 1).
      out_dir (str): Base directory to save downloaded files.
    """
    # Construct the API URL for the study.
    # The API expects an accession like "667.1" for study OSD-667 version 1.
    api_url = f"https://osdr.nasa.gov/osdr/data/osd/files/{study_id}.{version}"
    print(f"Fetching metadata for study OSD-{study_id} (version {version}) ...")
    
    r = requests.get(api_url)
    if r.status_code != 200:
        print(f"Error: Could not fetch study OSD-{study_id} (HTTP {r.status_code})")
        return
    
    data = r.json()
    study_key = f"OSD-{study_id}"
    if study_key not in data.get("studies", {}):
        print(f"Study key {study_key} not found in API response.")
        return

    study_files = data["studies"][study_key].get("study_files", [])
    # Filter to only files with a .fastq.gz extension.
    fastq_files = [f for f in study_files if f["file_name"].lower().endswith(".fastq.gz")]
    if not fastq_files:
        print(f"No .fastq.gz files found for study OSD-{study_id}.")
        return

    # Create a directory for the study if it doesn't exist.
    study_dir = os.path.join(out_dir, study_key)
    os.makedirs(study_dir, exist_ok=True)

    base_download_url = "https://osdr.nasa.gov"
    for file_meta in fastq_files:
        file_name = file_meta["file_name"]
        remote_url = file_meta["remote_url"]
        download_url = base_download_url + remote_url
        print(f"Downloading {file_name} from {download_url} ...")
        # Stream download to avoid loading entire file in memory
        with requests.get(download_url, stream=True) as resp:
            if resp.status_code == 200:
                file_path = os.path.join(study_dir, file_name)
                with open(file_path, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                print(f"Downloaded: {file_path}")
            else:
                print(f"Failed to download {file_name}: HTTP {resp.status_code}")

if __name__ == "__main__":
    # List your study IDs here (without the "OSD-" prefix)
    study_ids = ["667", "666", "665", "532", "525", "515", "512"]  # Example: OSD-667, OSD-86, and OSD-87
    for sid in study_ids:
        download_fastq_files_for_study(sid)
