"""Upload geometry optimization prediction files from all models to a Figshare article.

This script creates/updates a Figshare article containing ML-relaxed structures from all
models in the Matbench Discovery benchmark.
"""

import os
import traceback

from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.figshare import (
    BASE_URL,
    create_article,
    make_request,
    upload_file_to_figshare,
)
from matbench_discovery.models import MODEL_METADATA

__author__ = "Assistant"
__date__ = "2024-03-21"


def get_existing_files(article_id: int) -> dict[str, int]:
    """Get a mapping of filenames to file IDs for files already in the article."""
    files = make_request("GET", f"{BASE_URL}/account/articles/{article_id}/files")
    return {file["name"]: file["id"] for file in files}


def main() -> int:
    """Main function to upload geometry optimization files to Figshare."""
    metadata = {
        "title": "Matbench Discovery ML-Relaxed Structures",
        "description": """
        This dataset contains ML-relaxed structures from various models evaluated on the
        Matbench Discovery benchmark. Each file contains structures that were relaxed
        using a specific ML model's predicted forces and energies.

        For more information about the benchmark and models, visit:
        https://matbench-discovery.materialsproject.org
        """.strip(),
        "defined_type": "dataset",
        "tags": [
            "materials-science",
            "machine-learning",
            "crystal-structure-prediction",
        ],
        "categories": [
            25162,  # Structure and dynamics of materials
            25144,  # Inorganic materials (incl. nanomaterials)
            25186,  # Cheminformatics and Quantitative Structure-Activity Relationships
        ],
    }

    try:
        # Try to find existing article first
        articles = make_request(
            "GET", f"{BASE_URL}/account/articles?search={metadata['title']}"
        )
        article_id = articles[0]["id"] if articles else create_article(metadata)
        article_url = f"https://figshare.com/articles/{article_id}"

        print(f"\nArticle: {metadata['title']}")
        print(f"URL: {article_url}")

        existing_files = get_existing_files(article_id)
        print(f"Found {len(existing_files)} existing files")

        uploaded_files: dict[str, str] = {}
        new_files: dict[str, str] = {}

        pbar = tqdm(MODEL_METADATA)
        for model_name in pbar:
            geo_opt_metrics = (
                MODEL_METADATA[model_name].get("metrics", {}).get("geo_opt", {})
            )
            if not isinstance(geo_opt_metrics, dict):
                continue
            pred_file = geo_opt_metrics.get("pred_file")
            if not pred_file:
                continue
            if not os.path.isfile(pred_file):
                print(
                    f"Warning: ML relaxed structures for {model_name} not "
                    f"found, expected at {pred_file}"
                )
                continue
            pbar.set_description(f"Uploading {model_name}")

            filename = os.path.basename(pred_file)
            if filename in existing_files:
                file_id = existing_files[filename]
                file_url = f"https://figshare.com/ndownloader/files/{file_id}"
                uploaded_files[filename] = file_url
                continue

            file_id = upload_file_to_figshare(article_id, f"{ROOT}/{pred_file}")
            file_url = f"https://figshare.com/ndownloader/files/{file_id}"
            uploaded_files[filename] = file_url
            new_files[filename] = file_url

        print(f"\nTotal files: {len(uploaded_files)}")
        print(f"Newly added: {len(new_files)}")

        print("\nAll uploaded files:")
        for filename, url in uploaded_files.items():
            print(f"{filename}: {url}")

        if new_files:
            print("\nNewly added files:")
            for filename, url in new_files.items():
                print(f"{filename}: {url}")
        else:
            print("\nNo new files were added.")

        # Save URLs to YAML file
        yaml_path = f"{ROOT}/matbench_discovery/geo-opt-files.yml"
        with open(yaml_path, "w") as file:
            round_trip_yaml.dump(uploaded_files, file)

        print(f"\nSaved file URLs to {yaml_path}")

    except Exception:  # prompt to delete article if something went wrong
        if article_id := int(locals().get("article_id", 0)):
            model_name = locals().get("model_name", "")
            answer = ""
            while answer.lower() not in ("y", "n"):
                answer = input(
                    f"Upload failed for {model_name}. Delete article? [y/n]\n"
                    f"\n{traceback.format_exc()}"
                )
            if answer.lower() == "y":
                make_request("DELETE", f"{BASE_URL}/account/articles/{article_id}")
        raise

    else:
        return 0


if __name__ == "__main__":
    main()
