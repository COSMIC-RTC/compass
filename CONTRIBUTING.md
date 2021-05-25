## General code guidelines
- [PEP8](https://www.python.org/dev/peps/pep-0008/) naming convention
- Explicit variable name
- All class attributes must be declared in the __init__ method of the class
- Force explicit use of keyword arguments
- Documentation is generated with doxygen (C++) and doxypypy (Python)
- Use [typing](https://docs.python.org/3/library/typing.html) module to give hint on the method/function signatures

## PEP 8 naming convention

Basicaly, it can be summarized like this:

  - class name = CamelCase
  - function and variable names = snake_case

For example:

```python
class MyBeautifulClass:
    """ Class description
    """
    def my_favorite_function(self, my_variable_with_explicit_name):
        """ snake_case for functions and variable name
        Explicit names are required : avoid single letter variable for example
        """
```

All names must be explicit.

For example, ```applyVoltGetSlopes``` instead of ```apply_volt_and_get_slopes```

This convention is also the baseline for the C++ API.

More info about [PEP8](https://www.python.org/dev/peps/pep-0008/) and [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)

## Explicit variable names
All the variables declared in the code must have an explicit name which is coherent with its assignation.
Then :
- Avoid single letter variable (exception can be made with "for" loop index)
- Only use an acronym if it is used widely

## Class attributes
All the attributes of a class must be declared in its __init__ method.
If it is assigned by an other method, declares it as None in the __init__ method.

## Force explicit use of keyword arguments
If a method/function use keyword arguments, force to use explicitly each keyword when a call is made. Python 3 allows to do it easily using ```*``` :
```python
def my_favorite_function(self, variable_1, *, kwarg_1=None, kwarg_2=None):
    ...
```
With this convention, calling ```my_favorite_function(0, 1)``` will result in an error as the code expects only one positional argument.
Then, the correct way to call this function is ```my_favorite_function(0, kwarg_1=1)```, and this is much more explicit and easy to read and understand.

Exception can be made for method/function using only one keyword argument without any positional argument :
```python
def my_favorite_function(self,kwarg_1=None):
    ...
```


## Code documentation with doxypypy and typing

All the class and method must be documented following a doxypypy compatible format.
In the class documentation, attributes of the class must be listed and described :

```python
class MyBeautifulClass(MyBeautifulParentClass):
    """ Brief description of the class

    Detailed description

    Attributes inherited from MyBeautifulParentClass:
        titi : (int) : titi description

    Attributes:
        toto: (dict) : toto description

        tata : (float) : tata description
    """

```

In the method documentation, follow the docstring template below:

```python
def my_favorite_function(self, file_path : str, wfs_index : int,
                         optional_argument : bool=False) -> int:
    """ Brief description of the function

    Detailed description

    Args:
        file_path : (str) : parameter description

        wfs_index : (int) : parameter description

    Kwargs:
        keyword_argument : (bool) : parameter description. Default value is False

    Returns:
        something : (int) : return description

    Raises:
        ZeroDivisionError, AssertionError, & ValueError.

    Examples:
        >>> my_favorite_function(file_path, wfs_index)
        0
    """
```

The [typing](https://docs.python.org/3/library/typing.html) python module must be used to also described the variable and return type in function signature, as in the example above.
Respect the docstring template described above, including attention to alinea, line break between parameters...

As COMPASS had become quite huge with time, code documentation will be a long term objective which will require contributions from everyone : do not hesitate to reformat, complete or modify docstrings that do not respect the above rules. We are counting on you all...

For Visual Studio code users, we recommends to use [Doxygen Documentation Generator](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen). It really useful, it will generate a template based on the function signature.

[more info](http://www.doxygen.nl/manual/docblocks.html)

## Merge request (test phase)

The develop branch is protected, which means that you cannot push on it directly. The typical process is as follows:

* git co -b myname/mybranch # create your branch
* ... # make your changes
* git ci -m "my message" # as much as you want
* git push --set-upstream origin myname/mybranch # as much as you want

When you are satisfied with your changes, in the message, it offers to make a merge request on develop of the style:

```bash
remote:
remote: To create a merge request for myname/mybranch, visit:
remote: https://gitlab.obspm.fr/compass/compass/-/merge_requests/new?merge_request%5Bsource_branch%5D=nono%2Ftest
remote:
```

If you have already committed to develop, you will not be able to push them on gitlab.
To avoid losing your changes in progress:

* git co -b myname/mybranch
* git push --set-upstream origin myname/mybranch # as much as you want

then make your merge request as defined above

## Citation
If you use COMPASS to publish, please cite one of those reference:
- [Ferreira, F. et al, “COMPASS: an efficient GPU-based simulation software for adaptive optics systems”, HPCS 2018](https://doi.org/10.1109/HPCS.2018.00043)
- [Ferreira, F. et al., “Real-time end-to-end AO simulations at ELT scale on multiple GPUs with the COMPASS platform”, SPIE 2018](https://doi.org/10.1117/12.2312593)
- [Gratadour, D. et al, “GPUs for adaptive optics: simulations and real-time control”, SPIE, 2012](https://doi.org/10.1117/12.925723)
