from setuptools import setup

setup(name="gspn",
    version="0.1",
    description="Discrete stochastic processes in continuous time",
    long_description="""These classes implement a basic discrete
    event simulation in continuous time.
    """,
    classifiers=[
      "Development Status :: 3 - Alpha",
      "Intended Audience :: Science/Research",
      "Natural Language :: English",
      "License :: OSI Approved :: BSD License",
      "Programming Language :: Python :: 3.4",
      "Topic :: Scientific/Engineering :: Mathematics"
    ],
    keywords="stochastic dynamic simulation Markov Gillespie",
    url="http://github.com/adolgert/PyGSPN",
    author="Drew Dolgert",
    author_email="ajd27@cornell.edu",
    license="BSD 3-clause",
    packages=["gspn"],
    install_requires=[
        "numpy",
        "scipy"
    ],
    include_package_date=True,
    zip_safe=False)
