# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PackLab/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| PackLab/analytical/grid.py          |       16 |        3 |        6 |        3 |     72.73% |36, 41, 43 |
| PackLab/monte\_carlo/diagnostics.py |       86 |        8 |       24 |       10 |     83.64% |43, 51, 56-\>60, 109, 111, 113, 149-\>151, 158, 172, 186 |
| PackLab/monte\_carlo/persistence.py |      135 |       40 |       56 |       22 |     64.40% |72, 86-89, 93, 99-100, 110, 118-122, 135, 137, 139, 141, 143, 145, 147, 149, 152-154, 194, 198, 223-224, 259, 263-264, 271, 275-278, 281, 283, 285, 287, 289 |
| PackLab/monte\_carlo/results.py     |      158 |       47 |       22 |        5 |     68.89% |62-64, 142-143, 155, 167, 178, 190, 202, 219-253, 271-280, 292, 318, 320, 360-363, 368-371 |
| PackLab/monte\_carlo/structure.py   |       74 |       18 |       32 |       15 |     66.98% |52-53, 57, 63, 69, 76, 78, 130, 136, 139, 141, 146, 157, 162, 164-166, 176, 180-\>189 |
| PackLab/scattering/data.py          |      156 |       46 |       54 |       21 |     66.19% |37-\>39, 97, 100, 106, 108, 112, 125-126, 128, 130, 132, 141, 146-147, 149, 151, 153, 155, 157, 163, 231, 269-280, 318-331, 357-363, 393, 395, 397, 404 |
| PackLab/scattering/model.py         |       33 |       10 |       10 |        3 |     69.77% |66-74, 106, 111 |
| PackLab/scattering/plottings.py     |       74 |       68 |       24 |        0 |      6.12% |18-25, 75-124, 165-205 |
| **TOTAL**                           |  **735** |  **240** |  **228** |   **79** | **62.72%** |           |

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