# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|-------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PackLab/analytical/grid.py      |       16 |        3 |        6 |        3 |     72.73% |36, 41, 43 |
| PackLab/monte\_carlo/results.py |      153 |       45 |       22 |        5 |     69.14% |118-119, 131, 143, 154, 166, 178, 196-232, 248-257, 267, 286, 288, 328-331, 336-339 |
| PackLab/scattering/data.py      |      156 |       46 |       54 |       21 |     66.19% |37-\>39, 97, 100, 106, 108, 110, 121-122, 124, 126, 128, 137, 142-143, 145, 147, 149, 151, 153, 159, 227, 265-276, 310-323, 347-351, 381, 383, 385, 389 |
| PackLab/scattering/model.py     |       32 |       10 |       10 |        3 |     69.05% |61-69, 102, 107 |
| PackLab/scattering/plottings.py |       74 |       68 |       24 |        0 |      6.12% |17-24, 71-120, 161-201 |
| **TOTAL**                       |  **434** |  **172** |  **116** |   **32** | **57.09%** |           |

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