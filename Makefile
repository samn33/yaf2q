test:
	cargo test; cd tests; pytest -s .

ruff:
	ruff check .

install:
	pip install .

uninstall:
	pip uninstall yaf2q
