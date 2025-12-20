#!/bin/bash
sphinx-apidoc -f -o ./source ../yaf2q
make html
cp -r ./build/html/* ../docs
touch ../docs/.nojekyll
