import requests
import time
from typing import List
import json

def read_uniprot_ids(filename: str) -> List[str]:
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def fetch_protein_info(uniprot_id: str) -> dict:
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # First, get the UID for the UniProt ID
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "protein",
        "term": f"{uniprot_id}[Accession]",
        "retmode": "json"
    }
    
    search_response = requests.get(search_url, params=search_params)
    search_data = search_response.json()
    
    if not search_data.get("esearchresult", {}).get("idlist"):
        return {"error": f"No results found for {uniprot_id}"}
    
    uid = search_data["esearchresult"]["idlist"][0]
    
    # Then fetch the protein information
    fetch_url = f"{base_url}efetch.fcgi"
    fetch_params = {
        "db": "protein",
        "id": uid,
        "retmode": "json"
    }
    
    fetch_response = requests.get(fetch_url, params=fetch_params)
    return fetch_response.json()

def main():
    # Read UniProt IDs
    uniprot_ids = read_uniprot_ids("member-accessions-A8PE21 (2).txt")
    
    # Create a dictionary to store results
    results = {}
    
    # Fetch information for each ID
    for uniprot_id in uniprot_ids:
        print(f"Fetching information for {uniprot_id}...")
        try:
            protein_info = fetch_protein_info(uniprot_id)
            results[uniprot_id] = protein_info
            
            # Add a small delay to avoid overwhelming the API
            time.sleep(0.34)  # NCBI allows up to 3 requests per second
            
        except Exception as e:
            print(f"Error fetching {uniprot_id}: {str(e)}")
            results[uniprot_id] = {"error": str(e)}
    
    # Save results to a JSON file
    with open("protein_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("Fetching complete. Results saved to protein_results.json")

if __name__ == "__main__":
    main() 