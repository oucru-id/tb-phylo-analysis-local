#!/usr/bin/env python3

import argparse
import base64
import json
import os
import sys
from datetime import datetime, timezone

import requests


def get_headers(access_token):
    return {
        "Accept": "application/fhir+json",
        "Authorization": f"Bearer {access_token}",
    }


def load_access_token(token_file):
    with open(token_file, "r") as f:
        data = json.load(f)
    access_token = data.get("access_token")
    if not access_token:
        raise ValueError(f"No access_token found in {token_file}")
    expires_at = data.get("expires_at")
    if expires_at:
        try:
            exp = datetime.fromisoformat(expires_at)
            if datetime.now(timezone.utc) > exp:
                print(
                    f"WARNING: Access token expired at {expires_at}. "
                    "Run get_access_token.py to refresh.",
                    file=sys.stderr,
                )
        except ValueError:
            pass
    return access_token


def upload_binary(base_url, access_token, matrix_path):

    with open(matrix_path, "rb") as f:
        b64_data = base64.b64encode(f.read()).decode("ascii")

    binary_resource = {
        "resourceType": "Binary",
        "contentType": "text/tab-separated-values",
        "data": b64_data,
    }

    headers = get_headers(access_token)
    headers["Content-Type"] = "application/fhir+json"

    url = f"{base_url}/Binary"
    print(f"Uploading Binary resource to {url} ...")
    resp = requests.post(url, headers=headers, json=binary_resource)
    resp.raise_for_status()

    result = resp.json()
    binary_id = result.get("id")
    if not binary_id:
        raise ValueError(f"FHIR server did not return an id for Binary resource: {result}")

    binary_url = f"{base_url}/Binary/{binary_id}"
    print(f"Binary resource created: {binary_url}")
    return binary_url


def upload_document_reference(base_url, access_token, binary_url, title):
    doc_ref = {
        "resourceType": "DocumentReference",
        "meta": {
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "tb-distance-matrix",
                    "display": "Genomics Tuberculosis Distance Matrix",
                }
            ]
        },
        "status": "current",
        "docStatus": "final",
        "description": f"Distance matrix for TB Genomics ({title})",
        "content": [
            {
                "attachment": {
                    "contentType": "text/tab-separated-values",
                    "url": binary_url,
                    "title": title,
                }
            }
        ],
    }

    headers = get_headers(access_token)
    headers["Content-Type"] = "application/fhir+json"

    url = f"{base_url}/DocumentReference"
    print(f"Creating DocumentReference at {url} ...")
    resp = requests.post(url, headers=headers, json=doc_ref)
    resp.raise_for_status()

    result = resp.json()
    doc_id = result.get("id")
    print(f"DocumentReference created: {base_url}/DocumentReference/{doc_id}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Upload distance_matrix.tsv to FHIR server as DocumentReference"
    )
    parser.add_argument("--matrix", required=True, help="Path to distance_matrix.tsv")
    parser.add_argument("--url", required=True, help="FHIR server base URL")
    parser.add_argument("--token-file", required=True, help="Path to access_token.json")
    args = parser.parse_args()

    if not os.path.isfile(args.matrix):
        print(f"ERROR: Matrix file not found: {args.matrix}", file=sys.stderr)
        sys.exit(1)

    base_url = args.url.rstrip("/")
    title = os.path.basename(args.matrix)

    print(f"Loading access token from: {args.token_file}")
    access_token = load_access_token(args.token_file)

    binary_url = upload_binary(base_url, access_token, args.matrix)
    result = upload_document_reference(base_url, access_token, binary_url, title)

    with open("upload_result.json", "w") as f:
        json.dump(result, f, indent=2)

    print("Upload complete. Result saved to upload_result.json")


if __name__ == "__main__":
    main()
