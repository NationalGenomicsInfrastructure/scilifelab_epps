# scilifelab_epps

Collection of EPPs for the Scilifelab Stockholm node.

Table of contents

- [scilifelab_epps](#scilifelab_epps)
  - [Overview](#overview)
  - [Installation](#installation)
  - [Development](#development)
    - [Automated linting](#automated-linting)
      - [GitHub Actions](#github-actions)
      - [Pre-commit](#pre-commit)
      - [VS Code automation](#vs-code-automation)
    - [Git blame suppression](#git-blame-suppression)

## Overview

The repo is installed inside a designated conda environment on a server hosting Illumina Clarity LIMS.

After installation, the scripts in `scripts/` are addded to the `bin` of the Conda environment, from which they can be called by LIMS.

## Installation

Inside the repo, run `pip install .`

## Development

Run `pip install requirements-dev.txt` to install packages used for development and `pip install -e .` to make the installation editable.

### Automated linting

This repo is configured for automated linting. Linter parameters are defined in `pyproject.toml`.

As of now, we use:

- [ruff](https://docs.astral.sh/ruff/) to perform automated formatting and a variety of lint checks.
  - Run with `ruff check .` and `ruff format .`
- [mypy](https://mypy.readthedocs.io/en/stable/) for static type checking and to prevent contradictory type annotation.
  - Run with `mypy **/*.py`
- [pipreqs](https://github.com/bndr/pipreqs) to check that the requirement files are up-to-date with the code.
  - This is run with a custom Bash script in GitHub Actions which will only compare the list of package names.

    ```
    # Extract and sort package names
    awk '{print $1}' $1 | sort -u > "$1".compare
    awk -F'==' '{print $1}' $2 | sort -u > "$2".compare

    # Compare package lists
    if cmp -s "$1".compare "$2".compare
    then
      echo "Requirements are the same"
      exit 0
    else
      echo "Requirements are different"
      exit 1
    fi
    ```

- [prettier](https://prettier.io/) to format common languages.
  - Run with `prettier .`
- [editorconfig-checker](https://github.com/editorconfig-checker/editorconfig-checker) to enforce `.editorconfig` rules for all files not covered by the tools above.
  - Run with
    ```
    editorconfig-checker $(git ls-files | grep -v '.py\|.md\|.json\|.yml\|.yaml\|.html')
    ```

#### [GitHub Actions](https://docs.github.com/en/actions)

Configured in `.github/workflows/lint-code.yml`. Will test all commits in pushes or pull requests, but not change code or prevent merges.

#### [Pre-commit](https://pre-commit.com/)

Will prevent local commits that fail linting checks. Configured in `.pre-commit-config.yml`.

To set up pre-commit checking:

1. Run `pip install pre-commit`
2. Navigate to the repo root
3. Run `pre-commit install`

This can be disabled with `pre-commit uninstall`

#### VS Code automation

To enable automated linting in VS Code, go the the user `settings.json` and include the following lines:

```
"[python]": {
    "editor.defaultFormatter": "charliermarsh.ruff",
}
```

This will run the `ruff`-mediated linting with the same parameters as the `GitHub Actions` and `pre-commit` every time VS Code is used to format the code in the repository.

To run formatting on save, include the lines:

```
"[python]": {
    "editor.formatOnSave": true,
}
```

### Git blame suppression

When a non-invasive tool is used to tidy up a lot of code, it is useful to supress the Git blame for that particular commit, so the original author can still be traced.

To do this, add the hash of the commit containing the changes to `.git-blame-ignore-revs`, headed by an explanatory comment.
