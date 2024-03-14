clean:
	rm -rf dist/ build/ *.egg-info
build:
	python -m build

upload:
	twine upload dist/*
format:
	black nanomux/*.py