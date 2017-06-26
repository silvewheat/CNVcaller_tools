pip uninstall cnvcallertools
python setup.py sdist
python setup.py bdist_wheel
pip install dist/cnvcallertools-0.1.0-py3-none-any.whl