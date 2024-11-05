"""Update auto-generated API docs viewable on the site's /api page."""

import json
import os
from glob import glob

from lazydocs import generate_docs

SITE = f"{os.path.dirname(__file__)}/../site"

with open(f"{SITE}/package.json") as file:
    pkg = json.load(file)  # get repo URL from package.json

out_path = f"{SITE}/src/routes/api"

for path in glob(f"{out_path}/*.md"):
    os.remove(path)

generate_docs(
    ["matbench_discovery"],
    output_path=out_path,
    watermark=False,
    src_base_url=f"{pkg['repository']}/blob/-",
)

# Tweak lazydocs's markdown output:
for path in glob(f"{out_path}/*.md"):
    with open(path) as file:
        text = file.read()

    # remove bold tags since they break inline code
    text = text.replace("<b>", "").replace("</b>", "")

    # make badges linking to GitHub source code blue with flat style, add alt text
    text = text.replace(
        'src="https://img.shields.io/badge/-source-cccccc?style=flat-square"',
        'src="https://img.shields.io/badge/source-blue?style=flat" alt="source link"',
    )
    with open(path, mode="w") as file:
        file.write(text)
