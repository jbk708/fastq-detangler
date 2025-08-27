.PHONY: help install install-dev test lint format check clean build

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package in development mode
	pip install -e .

install-dev:  ## Install development dependencies
	pip install -r requirements-dev.txt
	pip install -e .

test:  ## Run tests with pytest
	pytest tests/ -v --cov=fastq_detangler --cov-report=term-missing

test-fast:  ## Run tests quickly without coverage
	pytest tests/ -v

lint:  ## Run flake8 linting
	flake8 fastq_detangler tests --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 fastq_detangler tests --count --exit-zero --max-complexity=10 --max-line-length=88 --statistics

format:  ## Format code with black
	black fastq_detangler tests

format-check:  ## Check if code is formatted correctly
	black --check --diff fastq_detangler tests

sort-imports:  ## Sort imports with isort
	isort fastq_detangler tests

sort-imports-check:  ## Check if imports are sorted correctly
	isort --check-only --diff fastq_detangler tests

type-check:  ## Run mypy type checking
	mypy fastq_detangler --ignore-missing-imports

check: format-check sort-imports-check lint type-check test  ## Run all checks

clean:  ## Clean up build artifacts and cache
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build:  ## Build the package
	python -m build

build-check:  ## Check the built package
	twine check dist/*

install-pre-commit:  ## Install pre-commit hooks
	pre-commit install

run-pre-commit:  ## Run pre-commit on all files
	pre-commit run --all-files

ci: check build-check  ## Run CI checks locally 