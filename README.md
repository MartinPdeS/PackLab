# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |     Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | --------: | --------: |
| PackLab/analytical.py               |      212 |      212 |       18 |        0 |     0.00% |     1-443 |
| PackLab/analytical/distributions.py |      130 |      130 |       32 |        0 |     0.00% |     1-319 |
| PackLab/analytical/domain.py        |      218 |      218 |       92 |        0 |     0.00% |     1-549 |
| PackLab/analytical/solver.py        |      168 |      168 |       18 |        0 |     0.00% |     1-342 |
| PackLab/monte\_carlo/results.py     |      124 |       58 |       26 |        5 |    48.67% |71, 75, 79, 97-123, 139-148, 172, 174, 184, 193-200, 204-206, 260-305 |
| PackLab/monte\_carlo/utils.py       |       30 |       23 |        6 |        0 |    19.44% |     26-55 |
| **TOTAL**                           |  **896** |  **809** |  **192** |    **5** | **8.64%** |           |

1 file skipped due to complete coverage.


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/MartinPdeS/PackLab/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/MartinPdeS/PackLab/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2FMartinPdeS%2FPackLab%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.