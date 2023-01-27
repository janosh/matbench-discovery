import os
from glob import glob

from lazydocs import generate_docs

from matbench_discovery import ROOT, URLs

# Update auto-generated API docs.

out_path = f"{ROOT}/site/src/routes/api"

for path in glob(f"{out_path}/*.md"):
    os.remove(path)

generate_docs(
    ["matbench_discovery"],
    output_path=out_path,
    watermark=False,
    src_base_url=f"{URLs['Repo']}/blob/-",
)

# Tweak lazydocs's markdown output:
# - remove bold tags since they break inline code
# - make badges linking to GitHub source code blue with flat style, add alt text
for path in glob(f"{out_path}/*.md"):
    text = open(path).read()

    text = text.replace("<b>", "").replace("</b>", "")

    text = text.replace(
        'src="https://img.shields.io/badge/-source-cccccc?style=flat-square"',
        'src="https://img.shields.io/badge/source-blue?style=flat" alt="source link"',
    )
    open(path, "w").write(text)
