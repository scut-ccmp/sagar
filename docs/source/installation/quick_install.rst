.. _quick_install:

=============
快速安装
=============

这一章将简述安装步骤，该包建议安装在Unix类型的系统上，在 ``OS X`` 和 ``Ubuntu`` 以及 ``Debian``,
``Archlinux`` 均测试安装, `Windows`系统也可以安装, 并且测试成功, 但是不建议开发者使用`Windows`安装维护。
该章节分为下列几个部分：

1. 依赖软件的安装
2. 安装 ``sagar`` 软件包
3. (option) jupyter notebook的安装


---------------
unix系统下的安装
---------------

依赖软件的安装
++++++++++++++++++
软件包的安装依赖于下列软件，因此在安装 ``sagar`` 之前需要首先安装下列依赖。

* `git`_ (To download the ``sagar`` package)
* `python >= 3.5`_ (The programming language used for sagar)
* `python-pip`_ (Python package manager)
* `virtualenv`_ (Software to create a virtual python environment to install sagar in)

.. _git: https://git-scm.com/downloads
.. _python >= 3.5: https://www.python.org/downloads
.. _python-pip: https://packaging.python.org/installing/#requirements-for-installing-packages
.. _virtualenv: https://packages.ubuntu.com/xenial/virtualenv

软件的安装方式取决于你所使用的系统.
对于 Ubuntu 和其他 Debian 系列的发行版，用户可以使用下列命令安装上述软件::

    $ sudo apt-get install git python3.6 python-pip virtualenv

对于 OS X 使用 `Homebrew`_ 作为包管理系统来安装所需软件::

    $ brew install git python

.. _Homebrew: http://brew.sh

开发者使用的开发环境为 ``Archlinux`` ，在Archlinux下的安装使用 `pacman`_ 软件包管理::

    $ sudo pacman -S python git virtualenv python-pip

.. _pacman: https://wiki.archlinux.org/index.php/pacman


安装 ``sagar`` 软件包
+++++++++++++++++++++++++++

当前仅支持从源码安装。创建安装地址所在的文件夹并且克隆软件所在仓库到本地，
运行下述命令::

    $ mkdir <your_directory>
    $ cd <your_directory>
    $ git clone https://github.com/unkcpz/sagar.git

为了避免安装该软件时同时安装的依赖包与系统正在使用的python包产生冲突， *强烈建议* 使用
virtual environment进行独立python环境的管理。
有关虚拟环境的管理请参考内容 :ref:`virtual environments <virtual_environment>`.
要建立一个独立的新的虚拟环境，运行下列命令::

    $ echo "export PYENV=$HOME/PYENV" >> ~/.bashrc
    $ source ~/.bashrc
    $ mkdir -p $PYENV
    $ virtualenv --python=/usr/bin/python3 $PYENV/sagar
    $ source $PYENV/sagar/bin/activate
    $ which python

以上命令将在家目录中建立文件夹 ``sagar``，并进入该环境. 此时你会看到你所使用的python
为 ``/home/<username>/PYENV/sagar/bin/python`` 。
在激活该环境后，你会看到在你的命令行提示符前出现 ``(sagar) $ <command>`` 这表示你正在当前
环境下工作，调用和安装的工具都会来自于以上目录中。

.. note:: 你可能需要更新你创建的新的虚拟环境中的安装工具 ``pip`` 和 ``setuptools`` ::

    (sagar) $ pip install -U setuptools pip

最后，在克隆的 ``sagar`` 项目所在的目录下执行安装::

    (sagar) $ pip install -e sagar

(option) jupyter notebook的安装
+++++++++++++++++++++++++++++++++++

 ``sagar`` 中提供了许多包装好的工具函数和类，提供给用户以组建新的需求。
 建议使用 ``jupyter notebook`` 来实现这些需求。jupyter的详细介绍请参考 `jupyter`_

在jupyter notebook中的使用案例，可以参考 :ref:`jupyter examples <examples>`

 .. _jupyter: http://jupyter.org/

安装jupyther,只需要在虚拟环境中执行::

    (sagar) $ pip install jupyter

(option) 运行nose2执行单元测试
+++++++++++++++++++++++++++++++++++++++++

进入sagar目录，执行::

    $ nose2 -v

测试通过则软件确保可用。

---------------
Windows系统下的安装
---------------

依赖软件的安装
+++++++++++++++++++++++++++

windows下的软件包安装同样依赖于上述软件，上述软件在anaconda软件中有着很好的管理，因此在安装 ``sagar`` 之前可以优先安装anaconda。

* `anaconda`_ (一个开源的Python发行版本，其包含了conda、Python等180多个科学包及其依赖项)

.. _anaconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/

首先下载最新版 ``anaconda`` 的windows32位或64, 安装时一直按确定或继续直到安装成功。

下载 ``sagar`` 软件包
++++++++++++++++++++++++++++

这里给出代码仓库的地址: `resource`_ ,

.. _resource: https://github.com/scut-ccmp/sagar

点击 ``Clone or download`` 中的 ``Download ZIP`` 选项,下载sagar软件包到本地.


安装Visual build tools
++++++++++++++++++++++++++++++++

该包由于依赖于`spglib`这个包, 需要C++的编译器进行编译, 但是一般windows上并没有
安装Visual Studio这样的编译器, 因此安装的时候会报
错: `visual studio C++ 14.0 is required`. 所以我们建议使
用 `visualcppbuildtoos.exe` 这个安装包来安装一系列所
需要的编译工具, 地址在这:  `vsbuildtools`_.

.. _vsbuildtools: https://github.com/ChangChunHe/Sundries/blob/master/visualcppbuildtools_full.exe

点击右上角的 ``Download`` 按钮.


安装 ``sagar`` 软件包
+++++++++++++++++++++++++++++++

首先在本地解压sagar-master.zip,进入sagar-master目录，打开 ``setup.py`` 文件，删除 ``"spglib==1.9.9.18"`` 部分。

然后在 ``cmd`` 下进入到sagar目录，执行::

	(根据自己的文件位置进行操作，例如)
	> G：
	> cd software/sagar-master/sagar-master

最后，执行::

	> pip install -e .

------------------------------
(option) 运行nose2执行单元测试
------------------------------

进入sagar目录，执行::

    $ nose2 -v

测试通过则软件确保可用。


现在你就可以开始使用软件包，或者继续查看后续手册 :ref:`get-started`.
