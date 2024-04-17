pandoc -V geometry:margin=0.5in -V geometry:paperwidth=13.5in  README.md -o README.pdf
zip -r domainator.zip README.md README.pdf pyproject.toml conda_env.yml test src -x "src/domainator.egg-info/*" "src/domainator/__pycache__/*" "test/__pycache__/*" "test/.ipynb_checkpoints/*" "test/data/.ipynb_checkpoints/*" "src/Bio/__pycache__/*"

# maybe create a directory in the zip?
# maybe instead of having all of those exclusions, make a temporary directory, copy in everything important, zip it, and delete it.

