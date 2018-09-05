.. _virtual_environment:

===================
python虚拟环境
===================

为什么要使用python虚拟环境?
+++++++++++++++++++++++++++++++++++

python的软件通常会依赖于其他的第三方python包，而且很多情况下，不同的软件有不同的依赖需求。
如果pyyabc以全局的方式安装，则它可能会升级或者降级已有的第三方python库，这将可能导致其他的软件不可用。
相对的，如果在安装了pyyabc之后又安装了其他的python库或软件，可能打破现有的依赖，这就可能导致pyyabc不可用。

因此，我们 *强烈* 建议对几种合作使用的软件进行环境独立，以避免软件的不可用。

什么是python虚拟环境?
++++++++++++++++++++++++++++++

一个特定python的虚拟环境就是一个包含有环境独立的所有python软件和库的文件夹，文件夹内包含了下列内容：

* python 可执行文件(executable)
* python 标准包(standard packages)
* python包管理软件 如 ``pip``
* 一个激活和设置环境变量参数的脚本，用来设置 ``PYTHONPATH`` 和 ``PATH`` 变量

激活环境的脚本确保了新创建的，指定环境文件夹下的python可执行文件在 ``PATH`` 环境变量的开头，以及所有python的软件包都在虚拟环境的文件夹中。
这就使得可以通过建立多个文件夹来独立python的使用环境。

并且，虚拟环境所依赖的文件夹可以在用户的目录下，这就使得用户创建虚拟环境和安装软件时无须去的管理员权限。

Creating a virtual environment
++++++++++++++++++++++++++++++
There are different programs that can create and work with virtual environments.
An example for python virtual environments is called ``virtualenv`` and can be installed with for example ``pip`` by running::

    $ pip install --user -U virtualenv

As explained before, a virtual environment is in essence little more than a directory containing everything it needs.
In principle a virtual environment can thus be created anywhere where you can create a directory.
You could for example opt to create a directory for all your virtual environments in your home folder::

    $ mkdir ~/.virtualenvs

Using ``virtualenv`` you can then create a new virtual environment with python 2.7 by running::

    $ virtualenv --python=<path/to/python2.7> ~/.virtualenvs/my_env

This will create the environment ``my_env`` and automatically activate it for you.
If you open a new terminal, or you have deactivated the environment, you can reactivate it as follows::

    $ ~/.virtualenvs/my_env/bin/activate

If it is activated successfully, you should see that your prompt is prefixed with the name of the environment::

    (my_env) $

To leave or deactivate the environment and set all the settings back to default, simply run::

    (my_env) $ deactivate


.. _aiida_path_in_virtualenv:

Creating a ``.aiida`` folder in your virtualenvironment
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

When you run AiiDA in multiple virtual environments, it can be convenient to use a separate ``.aiida`` folder for each virtualenv. To do this, you can use the :ref:`AIIDA_PATH mechanism <directory_location>` as follows:

1. Create your virtualenv, as described above
2. Create a ``.aiida`` directory in your virtualenv directory::

    $ mkdir ~/.virtualenvs/my_env/.aiida
3. At the end of ``~/.virtualenvs/my_env/bin/activate``, add the following line::

    export AIIDA_PATH='~/.virtualenvs/my_env'
4. Deactivate and re-activate the virtualenv
5. You can test that everything is set up correctly if you can reproduce the following::

    (my_env)$ echo $AIIDA_PATH
    >>> ~/.virtualenvs/my_env

    (my_env)$ verdi profile list
    >>> Configuration folder: /home/my_username/.virtualenvs/my_env/.aiida
    >>> Stopping: No configuration file found
    >>> Note: if no configuration file was found, it means that you have not run
    >>> 'verdi setup' yet to configure at least one AiiDA profile.
6. Continue setting up AiiDA with ``verdi setup`` or ``verdi quicksetup``.
