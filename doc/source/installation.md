# Installing GliderFlight

## Download

The software's repository is hosted on
[github](https://github.com/smerckel/gliderflight), from where you can
download the source using git or download it as a zip file.

## Installing

After having obtained the source code, you can build and install it
using [uv](https://docs.astral.sh/uv/) (or pip). From the root of the
repository:

```
$ uv build
$ uv pip install .
```

Alternatively, install directly with pip:

```
$ pip install .
```

## Dependencies

The gliderflight module depends on numpy and scipy.

The glidertrim script -- useful to check and adjust the ballast trim
of a glider during deployment -- additionally depends on matplotlib,
gsw, and dbdreader. All these packages are available from PyPI and will
be downloaded automatically when not present. However, you may prefer
to install the complex packages numpy, scipy and matplotlib using your
distribution's package manager (when on linux).

The dependencies are declared in `pyproject.toml`.

## PyPI

Glider flight is available from PyPI:

```
pip install gliderflight
```
