#!/usr/bin/env bash
# docs/build-site.sh [serve]: builds the docs site (Zensical, classic variant) from the repo's
# Markdown: docs/*.md plus README.md as the home page, copied into a staging folder with their
# links rewritten for the site, so each page has one source. Output: site/ at the repo root (not
# committed). `build-site.sh serve` previews it locally at http://127.0.0.1:8000/gtotree/
#
# Needs zensical (`pip install zensical`).
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
stage="$repo/.site-src"
rm -rf "$stage"; mkdir -p "$stage/images" "$stage/stylesheets"
for f in "$repo"/docs/*.md; do
    sed -e 's#](\.\./README\.md)#](index.md)#g' -e 's#](\.\./CHANGELOG\.md)#](https://github.com/AstrobioMike/GToTree/blob/master/CHANGELOG.md)#g' \
        "$f" > "$stage/$(basename "$f")"
done
cp -r "$repo"/docs/images/. "$stage/images/"
cp -r "$repo"/docs/stylesheets/. "$stage/stylesheets/"
# README: links into docs/ become site pages, and the repo-root files it points at go to GitHub
sed -e 's#](docs/images/#](images/#g' -e 's#src="docs/images/#src="images/#g' -e 's#href="docs/images/#href="images/#g' \
    -e 's#](docs/\([^)]*\))#](\1)#g' \
    -e 's#](CHANGELOG\.md)#](https://github.com/AstrobioMike/GToTree/blob/master/CHANGELOG.md)#g' \
    "$repo/README.md" > "$stage/index.md"
cd "$repo"
if [ "${1:-}" = "serve" ]; then
    zensical serve -f mkdocs.yml
else
    zensical build -f mkdocs.yml
fi
