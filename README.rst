NICEgame
========

E. Vayena, A. Chiappino-Pepe, H. MohammadiPeyhani, Y. Francioli, N. Hadadi, M. Ataman, J. Hafner, S. Pavlou, & V. Hatzimanikatis, A workflow for annotating the knowledge gaps in metabolic reconstructions using known and hypothetical reactions, Proc. Natl. Acad. Sci. U.S.A. 119 (46) e2211197119, https://doi.org/10.1073/pnas.2211197119 (2022).


Requirements - MATLAB
------------

You will need to have `Git LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/NICEgame.git /path/to/NICEgame
    cd /path/to/NICEgame
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2017b, and CPLEX 12.7 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_.

This module requires `matTFA <https://github.com/EPFL-LCSB/mattfa/>`_

Generating reduced models
-------------------------
1. Place the thermodynamic data for the corresponding orgnanism into the `matTFA thermoDatabases <https://github.com/EPFL-LCSB/matTFA/thermoDatabases>`_ folder.
2. Place the corresponding curated GEM into the `GEMs <https://github.com/EPFL-LCSB/redgem/GEMs>`_ folder.
3. Place the *get* file into the `runFileExample <https://github.com/EPFL-LCSB/redgem/runFileExample>`_  folder.
4. Run the *get* file


License
=======
The software in this repository is put under an APACHE licensing scheme


Requirements - python
------------


Further the following pip-python packages are required (can be found in detail in requirements.txt)


- bokeh>=0.12.1
- cobra>0.13
- equilibrator-api
- equilibrator-cache
- ipdb
- lxml
- networkx
- openpyxl
- pymysql
- pytest
- python-libsbml==5.11.4
- scipy
- sqlalchemy
- tabulate
- tqdm
- sphinx
- sphinx-rtd-theme

Container-based install
-----------------------

You might want to use this program inside of a container. The
|docker|_
subfolder has all the necessary information and source files to set it
up.

.. |docker| replace:: ``docker/``
.. _docker: https://github.com/EPFL-LCSB/NICEgame/python/nicegamepy/docker

.. code:: bash

    cd NICEgame/python/nicegamepy/docker
    ./build.sh
    ./run.sh

Building the docker image takes approximately 5 mins.


Setup
=====
If container-based installation is not preferred you can also install this module from source using ``pip``:
*For Python 3, you might have to use* ``pip3`` *instead of* ``pip``

.. code:: bash

    git clone https://github.com/EPFL-LCSB/NICEgame.git /path/to/NICEgame/python
    pip3 install -e /path/to/NICEgame

The installation process should not exceed a minute if the requirements are installed. If they are not, it might take longer as the installer installs them first.


