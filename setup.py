from setuptools import setup

setup(
    name='rnasleek',
    version='0.2.0',
    py_modules=['rnasleek'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'rnasleek = rnasleek:cli',
        ],
    },
)