# Contributing to featuretablefilter

Thank you for your interest in contributing to this package! This document provides guidelines and instructions for contributing.

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check existing issues. When creating a bug report, include as much detail as possible:

- **Use a clear and descriptive title**
- **Describe the steps to reproduce the bug**
- **Include what you expected to happen**
- **Include what actually happened**
- **Provide your R version and OS information**
- **Include a reproducible example if possible**

### Suggesting Enhancements

Enhancement suggestions are welcome. Please include:

- **A clear and descriptive title**
- **A detailed description of the suggested enhancement**
- **Examples of how the enhancement would be used**
- **Explanation of why this enhancement would be useful**

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following the code style guidelines
3. **Add tests** for new functionality
4. **Update documentation** including roxygen comments
5. **Ensure all tests pass** with `Rscript -e "testthat::test_dir('tests/testthat')"`
6. **Check package build** with `R CMD build .` and `R CMD check .`
7. **Submit a pull request**

## Code Style Guidelines

### R Code Style

- Use **spaces** around operators (`<-`, `=`, `+`, etc.)
- Use **4 spaces** for indentation (not tabs)
- Place opening braces `{` on the same line
- Place closing braces `}` on their own line
- Function names: `snake_case`
- Variable names: `snake_case`
- Maximum line length: 80 characters (100 for comments)

### Documentation

- All exported functions must have **roxygen2 documentation**
- Include `@param` for all arguments
- Include `@return` describing the return value
- Provide **examples** using `@examples`
- Use `\dontrun{}` for examples that require external data or take too long

### Testing

- Write **focused tests** - one test per file per functionality
- Use **testthat** framework (edition 3)
- Test both **success cases** and **error conditions**
- Keep tests **fast** and **reproducible**
- Name tests descriptively: `test_that("function does X", { ... })`

## Development Setup

```r
# Install dependencies
install.packages(c("devtools", "roxygen2", "testthat"))

# Install package dependencies
install.packages(c("ggplot2", "tidyr", "S4Vectors"))

# Optional dependencies for full functionality
install.packages(c("phyloseq", "vegan", "pheatmap"))
BiocManager::install(c("TreeSummarizedExperiment", "SingleCellExperiment"))

# Load package for development
devtools::load_all(".")

# Run tests
testthat::test_dir("tests/testthat")

# Document changes
roxygen2::roxygenise()

# Build package
devtools::build(".")

# Check package
devtools::check(".")
```

## Commit Messages

Follow [conventional commits](https://www.conventionalcommits.org/) format:

- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `style:` - Code style changes (formatting, etc.)
- `refactor:` - Code refactoring
- `test:` - Adding or updating tests
- `chore:` - Maintenance tasks

Example:
```
feat: add support for SingleCellExperiment objects

- Add from_SCE() conversion function
- Update run_filtering_pipeline() to handle SCE input
- Add tests for SCE conversion
- Update documentation
```

## License

By contributing, you agree that your contributions will be licensed under the
same GPL-3 license as this project.

## Questions?

If you have questions, please open an issue or contact
vojtech.barton@gmail.com.
