[index](README.md)
# Developing Domainator

## testing

run all tests
```
pytest test
```

examine test coverage
```
coverage run -m pytest test
coverage report -m
coverage html
```
then open htmlcov/index.html in a browser to see coverage


## Converting the manual into a pdf
on Ubuntu

Install pandoc
```
sudo apt-get install pandoc texlive-latex-base texlive-fonts-recommended texlive-extra-utils texlive-latex-extra librsvg2-bin
```

run pandoc
```
pandoc -V geometry:margin=0.5in -V geometry:paperwidth=13.5in  README.md -o README.pdf
```

## Projects that domainator depends on

  - hmmer3
  - usearch
  - python3
  - biopython
  - pytest
  - pytest-datadir
  - pandas
  - seaborn
  - cd-hit
  - scipy
  - pyhmmer
  - umap-learn
  - diamond
  - coverage
